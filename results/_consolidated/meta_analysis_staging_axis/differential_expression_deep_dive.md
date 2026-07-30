# Differential Expression Methods: A Deep-Dive Review

*Companion note to the psoriasis meta-analysis white paper. Purpose: a self-contained
reference on the statistical models behind RNA-seq differential expression (DE), how the three
dominant tools differ, and how to defend each choice. Grounded where possible in measurements
on our own data (SRP035988, 178 samples, Gencode v26).*

---

## 0. The problem DE has to solve

We have, per gene, a handful-to-many replicate counts in each group and we ask: *is the
difference between groups larger than what replicate-to-replicate noise would produce?* Two
features make this hard and shape every method below:

1. **Counts, not continuous measurements.** The data are integer read counts with a
   mean-dependent variance — ordinary linear models (which assume constant variance) do not
   apply directly.
2. **Many genes, few replicates.** Tens of thousands of genes, each with few replicates, make
   per-gene variance estimates unstable. Every modern method's core trick is **borrowing
   statistical strength across genes** (empirical Bayes).

---

## 1. Poisson vs. Negative Binomial

### 1.1 Poisson — the count baseline

The Poisson distribution models independent counting events (reads landing on a gene). It has
one parameter (mean λ) and a defining constraint: **variance = mean.** It captures only
**technical** sampling noise — the randomness of drawing reads from a fixed pool.

### 1.2 Why Poisson fails: over-dispersion

Real experiments add **biological variability**: replicate samples (different individuals,
states) genuinely differ in a gene's true expression. This pushes **variance above the mean** —
**over-dispersion**. Note "over-dispersed" is defined *relative to the Poisson tie
variance = mean*, not "variance is large in absolute terms": binomial counts are
*under*-dispersed (variance < mean), so the direction is not automatic.

**Measured on our data (NN group, 83 samples, no disease effect):** the median
variance/mean ratio across genes was **24.2** (Poisson predicts 1), and **99.5%** of genes had
variance > 2× mean. RNA-seq is decisively over-dispersed.

### 1.3 Negative Binomial = Gamma-Poisson

The NB fixes this by letting the Poisson rate itself vary between samples according to a
**Gamma distribution**. Marginalizing a Poisson whose rate is Gamma-distributed yields exactly
the Negative Binomial — hence "Gamma-Poisson." It has two parameters (mean + dispersion) with
variance:

```
Var = mu + phi * mu^2
        |        |
     Poisson   biological
     (technical) over-dispersion
```

As phi -> 0 the NB collapses back to Poisson. This two-term variance is the statistical heart
of edgeR and DESeq2.

**Defense:** *Poisson assumes variance = mean (technical noise only); our replicates show
variance ~24x the mean. The NB adds a dispersion term so Var = mu + phi*mu^2, matching the
data. Using Poisson would understate variance and inflate false positives.*

---

## 2. The dispersion parameter

**Definition.** Dispersion **phi** quantifies biological variability beyond Poisson. Reported
as the **biological coefficient of variation, BCV = sqrt(phi)**.

**Measured on our data:** edgeR common dispersion **phi = 0.124**, i.e. **BCV = 35%** — a
typical value for human/clinical tissue (cell lines ~10%, human cohorts 30-40%). Concretely at
a count of 1000, Poisson expects variance 1000 but the NB expects ~125,000: ignoring dispersion
understates the noise ~100-fold.

**Why it is the crux.** Dispersion is effectively the denominator of the DE test — it sets how
big a fold-change must be to be called real. With few replicates, per-gene phi is noisy, so all
tools **shrink** it toward a global trend (empirical Bayes):

- **edgeR:** common + trended (mean-dependent) + EB-shrunk tagwise dispersion.
- **DESeq2:** fits a smooth dispersion-vs-mean trend, shrinks per-gene estimates toward it,
  flags outliers.
- **limma-voom:** does not estimate an NB dispersion at all — it converts the mean-variance
  trend into observation weights (Section 3).

---

## 3. The three tools: models and tests

Two philosophies: **model counts directly with an NB GLM** (edgeR, DESeq2) vs. **transform to
continuous and use a weighted normal linear model** (limma-voom).

| | **edgeR** | **DESeq2** | **limma-voom** (our choice) |
|---|---|---|---|
| Data modeled | raw counts | raw counts | log2-CPM (continuous) |
| Distribution | Negative Binomial | Negative Binomial | Normal (after voom) |
| Model | NB GLM | NB GLM | weighted linear model |
| Variance handling | EB-shrunk dispersion | dispersion trend + shrinkage | voom precision weights + EB variance moderation |
| Test | quasi-likelihood F-test (or LRT) | Wald test (default) or LRT | moderated t-test (F for multi-contrast) |
| Effect-size shrinkage | (optional) | `lfcShrink()` (apeglm/ashr) | via variance moderation; not LFC toward 0 by default |
| Multiple testing | BH FDR | BH FDR + independent filtering | BH FDR |

### 3.1 What voom does (our pipeline)

1. Fit an ordinary linear model on log2-CPM to get residuals.
2. Fit the **mean-variance trend** (how residual SD depends on expression level).
3. Convert that trend into a **precision weight per observation** — reliable high-count
   observations get large weights, noisy low-count ones small weights.
4. Run a **weighted linear model** + **moderated t-test**: each gene's variance is EB-shrunk
   toward a global trend (borrowing strength across genes).

So limma-voom respects the count mean-variance structure through weights rather than through an
NB likelihood, and solves the small-n problem via variance moderation rather than dispersion
shrinkage.

### 3.2 Why limma-voom here

At **n = 178** the three methods converge; limma-voom is fast, well-calibrated (does not
over-reject), and the most flexible for the complex designs we will need later (batch and
study covariates in the meta-analysis). edgeR/DESeq2 have an edge at very small n; that is not
our regime.

**Defense:** *edgeR and DESeq2 fit NB GLMs to counts (quasi-likelihood F-test; Wald test).
limma-voom transforms to log2-CPM, encodes the mean-variance relationship as precision weights,
and fits a weighted linear model with a moderated t-test. All three share information across
genes; at n=178 they agree and limma-voom is fast and flexible for later covariates.*

---

## 4. Shrunken (moderated) fold-change

**Problem.** Low-count genes (e.g. 2 vs 6 reads) can show huge *raw* fold-changes that are pure
noise, dominating a naive ranking.

**Shrinkage.** Pull unreliable (low-count / high-variance) log-fold-changes **toward zero**
while leaving well-measured genes untouched — the effect-size analogue of dispersion/variance
shrinkage. Big shrunken LFCs are both large *and* trustworthy.

**Per tool.**
- **DESeq2:** explicit — `lfcShrink()` with apeglm/ashr; MA-plots display shrunken LFCs.
- **edgeR:** `predFC` / shrinkage via the QL framework.
- **limma-voom (ours):** significance is protected by **variance moderation** in the moderated
  t-test, so a gene cannot be significant on a noisy fold-change alone. The reported `logFC` is
  not shrunk toward zero by default, but our **|log2FC| > 1 AND adjusted p < 0.05** thresholds
  enforce the same outcome: significant genes are both well-estimated and biologically large.
  `treat()` (test against a fold-change threshold) is available if we want the moderation folded
  directly into the hypothesis.

**Defense:** *Shrunken fold-changes stop noisy low-count ratios from topping the list. DESeq2
does this via lfcShrink; limma-voom protects inference through variance moderation, and our
|log2FC|>1 + adjusted-p<0.05 thresholds ensure significant genes are both large and reliable.*

---

## 5. Multiple testing — from raw p to FWER, FDR, and q-values

We test **24,528 genes at once**. At a raw threshold of alpha = 0.05, pure chance alone would
hand us ~1,226 "significant" genes even if nothing were truly different. So the raw p-value,
correct for a *single* gene, cannot be the reporting threshold across a whole transcriptome.
There are two philosophies for fixing this, and it helps to see them side by side.

### 5.1 The p-value itself

Each gene's p-value answers one narrow question: *if this gene were truly unchanged between PP
and NN, how often would noise alone produce a test statistic at least this extreme?* Small p =
the observed difference is hard to blame on chance. It is a **per-gene** statement; it says
nothing about the other 24,527 tests, which is exactly why testing many genes needs a
correction.

### 5.2 FWER and Bonferroni — "not even one false positive"

The **Family-Wise Error Rate (FWER)** is the probability of making **even a single** false
positive across the whole family of tests. **Bonferroni** controls it by dividing alpha by the
number of tests: a gene is significant only if p < 0.05 / 24,528 = **2.0e-6**.

This is the right tool when a single false positive is expensive — a confirmatory clinical
endpoint, or a GWAS hit you will chase with costly follow-up. For a **discovery-stage
transcriptome screen it is far too strict**: we *expect* thousands of genes to change in
lesional skin, and Bonferroni's "zero tolerance for any error" throws away most of that real
signal to avoid a handful of false ones. **We do not use it.** It answers the wrong question
for our purpose.

### 5.3 FDR and Benjamini-Hochberg — "a controlled fraction of my hits may be wrong"

The **False Discovery Rate (FDR)** asks a different, more sensible question for discovery: *of
the genes I call significant, what fraction do I expect to be false positives?* Controlling the
FDR at 5% means we accept that ~5% of our reported hits may be spurious in exchange for
capturing far more of the true signal.

The **Benjamini-Hochberg (BH)** procedure operationalises it:

1. Rank all N p-values ascending: p(1) <= p(2) <= ... <= p(N).
2. Find the largest rank k for which p(k) <= (k / N) x alpha.
3. Declare genes 1..k significant.

The **`adj.P.Val` column limma reports IS the BH-adjusted p-value.** Every threshold we quote
("adj.P < 0.05") is BH-FDR. This is the field standard for RNA-seq.

> **One-liner.** Bonferroni controls "any false positive at all" (too strict here); BH-FDR
> controls "the fraction of false positives in the hit list" (right for discovery). We use
> BH-FDR < 0.05, plus a |log2FC| > 1 effect-size filter.

### 5.4 Storey q-values and pi0 — estimating how many genes are truly null

BH makes one silent, conservative assumption: that **all** N genes *could* be null — formally,
that the proportion of true nulls **pi0 = 1**. In a tissue as remodelled as lesional psoriatic
skin that is obviously pessimistic; most of the transcriptome really has moved.

**Storey's q-value** improves on BH by **estimating pi0 from the data**. The intuition is in the
p-value histogram (Figure 7A): truly-null genes produce a *flat, uniform* carpet of p-values,
while truly-changed genes pile up as a spike near zero. The height of the flat shelf, relative
to a uniform distribution, estimates pi0 — the fraction of the histogram that is "just null."

For our data **pi0 = 0.133**: only ~13% of the 24,528 genes are estimated to be truly
unchanged, i.e. ~87% (~21,300 genes) carry some real signal. The **q-value** is then the FDR
analogue of the p-value — the minimum FDR at which a given gene is still called significant —
and it relates to BH by

> **q = pi0 x (BH adjusted p).**

Because pi0 = 0.133 < 1, every q-value is ~7.5x smaller than the corresponding BH value, so
Storey is uniformly **more powerful** (Figure 7B: all points on or below the 1:1 line). At the
p-value level this recovers thousands of extra genes (21,809 vs 18,773 at FDR < 0.05). **But
once our biologically-motivated |log2FC| > 1 filter is applied, BH and Storey converge to
essentially the same list (3,477 vs 3,479 genes)** — a reassuring sign that our headline result
is driven by strong, large-effect biology, not by the fine print of the FDR estimator. We report
both columns (`adj.P.Val` and `qvalue`) for transparency.

![Figure 7. Left (A): the p-value distribution — a tall spike near zero (real signal) sitting on a flat shelf of null genes; the dashed red line marks the estimated null level, pi0 = 0.133. Right (B): Storey q-value versus BH adjusted p; every point lies on or below the dotted 1:1 line because q = pi0 x BH, so the q-value is uniformly the more powerful of the two.]({{artifact:art_d5cd7635-7490-47a0-a3ec-dac25cfb089f}})

### 5.5 A subtlety worth stating: do not multiply by pi0 twice

A natural exercise: *given pi0 = 0.133, how many of the 18,773 BH-significant genes are expected
false positives at FDR < 0.05?* The answer is simply **0.05 x 18,773 ~ 939 genes** — the FDR
threshold already *is* the expected false-positive fraction.

It is tempting to then "correct" this by pi0 (0.133 x 939), but **that would be wrong**. Here is
why. pi0 is used **upstream**, inside Storey's machinery, to *convert p-values into q-values in
the first place* — it is the ingredient that makes the q-value less conservative than BH. Once a
gene already has a q-value (or you are thresholding at a stated FDR), that pi0 correction is
**already baked in**. Multiplying by pi0 a second time would **double-count** the same
adjustment and *understate* the true number of expected false positives. The rule: **pi0 is
applied once, when the q-value is built — never again afterward.** At a stated FDR of 5% on a
list of 18,773 genes, the expected false positives are 5% of 18,773, full stop.

### 5.6 Looking ahead: "combining p-values" in the meta-analysis

The p-value returns in a *different* role once we pool multiple psoriasis studies. There, each
study gives the **same gene** its own p-value, and meta-analysis fuses them. Two families exist:

- **Combine the p-values directly** — **Fisher's method** (sum of -ln p across studies follows a
  chi-square) or **Stouffer's Z** (weighted sum of per-study z-scores). Simple; needs only each
  study's p and direction; but discards effect-size magnitude.
- **Combine the effect sizes** — pool the per-study log2 fold-changes with **inverse-variance
  weights** (random-effects meta-analysis), which retains magnitude and models between-study
  heterogeneity. This is generally preferred.

So "combining on p-values" (a phrase you will meet in meta-analysis papers) means the *first*
family — one meta p-value per gene from the per-study p-values. Our plan is the effect-size
route as primary, with a p-combination method as a sensitivity cross-check. (Full treatment in a
later meta-analysis deep-dive.)

---

## 6. Testing coefficients and building null distributions

This section collects the statistical machinery that turns a fitted model into a p-value — and
the *non-parametric* alternative that pathway analysis (GSEA) will rely on. These are the tools
a reviewer is most likely to probe, so each is defined from first principles.

### 6.1 Linear models as the generalization of the t-test

A two-group comparison can be written as a linear model:

> y = beta0 + beta1 * group + error,  with group = 0 for NN, 1 for PP.

Then **beta1 is exactly the difference in group means** (our log2 fold-change), and **testing
beta1 = 0 is algebraically identical to the two-sample t-test** — same p-value. This is precisely
what limma fits per gene with the design `~ grp`; the `PPvsNN` coefficient *is* beta1. So we have
already "done" the t-test — as the one-binary-predictor special case of a linear model.

The reason to keep the linear-model framing rather than call `t.test()` is that it **generalizes
without changing method**, by editing the *design matrix*:

- **Covariate adjustment:** `~ group + age + sex + RIN` removes nuisance variation while testing group.
- **Multi-study meta-analysis (upcoming):** `~ group + study` tests the disease effect *adjusted
  for study of origin* — how we will pool datasets without batch confounding.
- **More than two groups / ANOVA:** extra columns for lesional / non-lesional / healthy, tested jointly.
- **Interactions:** `~ group * sex` asks whether the disease effect differs by sex.
- **Arbitrary contrasts:** with `~ 0 + group`, a **contrast matrix** (`limma::makeContrasts`)
  poses any comparison — "PP vs NN", "(PP+PN)/2 vs NN" — from a single fitted model.

The unifying idea: **the t-test is the special case; the linear model is the general engine, and
"complex comparisons" are just different design/contrast matrices fed to the same `lmFit`.** This
is the concrete reason we chose limma-voom over a battery of pairwise tests.

### 6.2 Testing a coefficient: t-test vs Wald test

Once the model is fit, each coefficient has an estimate `beta_hat` and a standard error
`SE(beta_hat)`. Both classical tests start from the **same ratio**:

> statistic = beta_hat / SE(beta_hat)

They differ only in **which reference distribution** that ratio is compared against — which in
turn reflects **how the uncertainty in the variance estimate is treated**:

| | **t-test** | **Wald test** |
|---|---|---|
| Ratio | beta_hat / SE | beta_hat / SE (same) |
| Compared to | **t-distribution** (df = residual d.f.) | **Normal / z** (or chi-square), asymptotic |
| Treats variance as | *estimated* from finite data | effectively *known* (large-n limit) |
| Small-n behaviour | heavier tails, conservative, correct | anti-conservative (too many false positives) |
| Used by | **limma** (moderated t), classic t-test | **DESeq2**, GLM Wald |

The conceptual split: the **t-test admits it estimated the variance** from limited replicates and
pays for that honesty with heavier tails; the **Wald test assumes the variance is essentially
known**, justified only asymptotically (large n), so it uses the thinner normal. Consequently
**small samples favour the t-test** (the Wald/z is anti-conservative when n is small — the classic
demonstration is Hauck & Donner 1977 [24]), while **large samples make them converge** (the
t-distribution approaches the normal as df grows).

What we actually use is better than either: limma's **moderated t** applies **empirical-Bayes
shrinkage** of each gene's variance toward a pooled mean-variance trend [15], borrowing strength
across all 24,528 genes. This stabilises the denominator (effectively adding degrees of freedom
beyond a single gene's handful of replicates) — a t-test whose variance has been regularised. At
n = 178 DESeq2's Wald agrees closely; the distinction matters most at small n.

A third form worth naming: the **likelihood-ratio test (LRT)** compares the fit of a full vs a
reduced model rather than a single coefficient's ratio, and is the natural choice for testing
*several coefficients at once* ("does *any* study term matter?"). edgeR and DESeq2 both offer it.
The LRT, Wald, and score tests are three asymptotically-equivalent views of the same likelihood;
Buse 1982 [25] is the classic one-page geometric explanation a reviewer will recognise.

### 6.3 Empirical p-values

Everything above is **parametric**: compute a statistic, look it up against an assumed theoretical
curve (t, normal, chi-square). An **empirical p-value** instead **builds the null distribution
from the data by resampling**:

1. Compute the real statistic on the actual data, S_obs.
2. Generate many "null" datasets in which the tested effect is, by construction, absent (§6.4).
3. Recompute the statistic on each, giving a distribution of S_null.
4. The empirical p-value is the fraction of null statistics at least as extreme as observed:

> p_emp = ( #{ S_null >= S_obs } + 1 ) / ( n_permutations + 1 )

The **+1 in both numerator and denominator** counts the observed data as one valid arrangement;
it keeps the p-value from ever being exactly zero and makes the test exact rather than optimistic
(Phipson & Smyth 2010 [22] — "permutation P-values should never be zero"). A practical
consequence: with 1,000 permutations the **smallest resolvable p-value is 1/1001 ~ 1e-3** — you
cannot report a permutation p finer than the permutation count supports, which is why permutation
FDRs look "chunky."

### 6.4 Permutation = non-parametric significance

**Non-parametric** means "no assumed form for the null distribution." **Permutation** manufactures
that null by **shuffling labels to deliberately break the association being tested**, then
measuring how large an effect chance alone produces. In gene-set / GSEA terms, two modes:

- **Phenotype (sample-label) permutation** — reshuffle the PP/NN labels across the 178 samples,
  re-rank all genes, recompute each pathway's enrichment score; repeat ~1,000x. This **preserves
  gene-gene correlation structure**, so it is the statistically preferred mode, and our sample
  sizes (83 / 95) are ample for it.
- **Gene-set permutation** — sample random gene sets of matched size as the null; used as a
  fallback when samples are too few to permute phenotypes. It **ignores gene-gene correlation**.

Why not a formula? Because the enrichment score is a path-dependent running-sum
(Kolmogorov-Smirnov-like) statistic with **no clean theoretical null** — there is no table to look
it up in, so the null is built empirically by permutation (Subramanian et al. 2005 [23]). This is
the exact machinery the next stage of the project runs on: **GSEA ranks all 24,528 genes by their
moderated-t statistic (6.1-6.2), then uses permutation (6.4) to assign each pathway an empirical
p-value (6.3).**

---

## 7. One-paragraph summary defense

> *RNA-seq counts are over-dispersed (our data: variance ~24x mean; BCV 35%), so Poisson is
> inadequate and the Negative Binomial (Gamma-Poisson) is the natural count model, adding a
> dispersion term Var = mu + phi*mu^2. edgeR and DESeq2 fit NB GLMs and test with a
> quasi-likelihood F-test and a Wald test respectively; limma-voom instead transforms to
> log2-CPM, turns the mean-variance trend into precision weights, and uses a moderated t-test.
> All three borrow strength across genes — shrinking dispersion, variance, or fold-change — to
> cope with few replicates per gene. We use limma-voom (fast, well-calibrated, flexible at
> n=178), threshold on adjusted p < 0.05 and |log2FC| > 1, and rely on variance moderation plus
> those thresholds to keep noisy low-count genes out of the results.*

---

## References

[3] Chen Y, Lun ATL, Smyth GK. *From reads to genes to pathways: differential expression
analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.*
F1000Research. 2016;5:1438. doi:10.12688/f1000research.8987.2. *(edgeR QL workflow.)*

[4] Robinson MD, Oshlack A. *A scaling normalization method for differential expression
analysis of RNA-seq data.* Genome Biology. 2010;11:R25. doi:10.1186/gb-2010-11-3-r25. *(TMM.)*

[8] Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq
data with DESeq2.* Genome Biology. 2014;15:550. doi:10.1186/s13059-014-0550-8. *(DESeq2 model,
Wald test, lfcShrink.)*

[11] Robinson MD, McCarthy DJ, Smyth GK. *edgeR: a Bioconductor package for differential
expression analysis of digital gene expression data.* Bioinformatics. 2010;26(1):139-140.
doi:10.1093/bioinformatics/btp616. *(edgeR; NB model.)*

[12] McCarthy DJ, Chen Y, Smyth GK. *Differential expression analysis of multifactor RNA-Seq
experiments with respect to biological variation.* Nucleic Acids Research. 2012;40(10):4288-4297.
doi:10.1093/nar/gks042. *(NB dispersion estimation, BCV.)*

[13] Law CW, Chen Y, Shi W, Smyth GK. *voom: precision weights unlock linear model analysis
tools for RNA-seq read counts.* Genome Biology. 2014;15:R29. doi:10.1186/gb-2014-15-2-r29.
*(voom precision weights.)*

[14] Ritchie ME, Phipson B, Wu D, et al. *limma powers differential expression analyses for
RNA-sequencing and microarray studies.* Nucleic Acids Research. 2015;43(7):e47.
doi:10.1093/nar/gkv007. *(limma; moderated t-test, empirical Bayes.)*

[15] Smyth GK. *Linear models and empirical Bayes methods for assessing differential expression
in microarray experiments.* Statistical Applications in Genetics and Molecular Biology.
2004;3:Article 3. doi:10.2202/1544-6115.1027. *(Original empirical-Bayes moderated t-statistic.)*

[16] Anders S, Huber W. *Differential expression analysis for sequence count data.* Genome
Biology. 2010;11:R106. doi:10.1186/gb-2010-11-10-r106. *(DESeq NB model, dispersion-mean
trend.)*

[17] Soneson C, Delorenzi M. *A comparison of methods for differential expression analysis of
RNA-seq data.* BMC Bioinformatics. 2013;14:91. doi:10.1186/1471-2105-14-91. *(Benchmark
comparison of DE methods.)*

[18] Benjamini Y, Hochberg Y. *Controlling the false discovery rate: a practical and powerful
approach to multiple testing.* J R Stat Soc B. 1995;57(1):289-300. *(FDR / BH procedure.)*

[19] Storey JD, Tibshirani R. *Statistical significance for genomewide studies.* Proc Natl Acad
Sci USA. 2003;100(16):9440-9445. doi:10.1073/pnas.1530509100. *(q-value, pi0 estimation.)*

[20] Fisher RA. *Statistical Methods for Research Workers.* Oliver & Boyd, Edinburgh; 1925.
*(Fisher's method for combining p-values.)*

[21] Whitlock MC. *Combining probability from independent tests: the weighted Z-method is
superior to Fisher's approach.* J Evol Biol. 2005;18(5):1368-1373.
doi:10.1111/j.1420-9101.2005.00917.x. *(Stouffer/weighted-Z p-value combination.)*

[22] Phipson B, Smyth GK. *Permutation P-values should never be zero: calculating exact P-values
when permutations are randomly drawn.* Stat Appl Genet Mol Biol. 2010;9(1):Article 39.
doi:10.2202/1544-6115.1585. *(The +1 correction; exact empirical p-values.)*

[23] Subramanian A, Tamayo P, Mootha VK, et al. *Gene set enrichment analysis: a knowledge-based
approach for interpreting genome-wide expression profiles.* Proc Natl Acad Sci USA.
2005;102(43):15545-15550. doi:10.1073/pnas.0506580102. *(GSEA; permutation null for the
enrichment score.)*

[24] Hauck WW, Donner A. *Wald's test as applied to hypotheses in logit analysis.* J Am Stat
Assoc. 1977;72(360):851-853. doi:10.1080/01621459.1977.10479969. *(Wald test anti-conservative /
aberrant at small samples.)*

[25] Buse A. *The likelihood ratio, Wald, and Lagrange multiplier tests: an expository note.* Am
Stat. 1982;36(3a):153-157. doi:10.1080/00031305.1982.10482817. *(Geometric relationship of the
three classical tests.)*
