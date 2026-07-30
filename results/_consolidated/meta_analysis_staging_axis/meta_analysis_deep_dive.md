# Deep dive: multi-study meta-analysis of the psoriasis transcriptome

*Companion to section 11 of the white paper. This document develops the statistical machinery
behind cross-study pooling, one idea at a time, and records the design decisions and honest
caveats that the main text summarises. Independent reference numbering.*

---

## 1. Why pool at all? The single-study trap

A single RNA-seq study, however large, confounds two things we would like to separate: the
biology of the disease, and the idiosyncrasies of that one experiment - its patient population,
its biopsy and library protocol, its sequencing batch, its analyst's choices. When our anchor
study (SRP035988, 178 samples) reports STAT3 up with p = 5.8e-52, that number answers the
question "given *these* samples, how surely is STAT3 up *here*?" It does not answer "is STAT3 up
in psoriasis?" - the question a clinician actually cares about.

Meta-analysis is the formal machinery for the second question. It takes an effect measured
independently in several studies and asks whether the *ensemble* supports a consistent
conclusion. The payoff is not usually a smaller p-value (as we will see, it can be larger); the
payoff is **generalisability**. A gene that pools to significance across independent cohorts is a
property of the disease. A gene significant in one study but not the pool was, at least partly, a
property of that study.

There are two broad strategies, and we deliberately chose the more conservative one.

- **Mega-analysis (merge then model).** Concatenate all samples into one matrix, add a "study"
  covariate (or a batch correction such as ComBat), and run one DE model. Maximum nominal power,
  but it assumes the studies are exchangeable after a single additive correction - an assumption
  that fails when protocols differ in ways that interact with biology, and one that is hard to
  audit gene-by-gene.
- **Meta-analysis (model then merge).** Run DE *separately* in each study, then combine the
  per-study effect sizes statistically. Slightly less nominal power, but every study keeps its
  own normalization and dispersion estimate, heterogeneity becomes an explicit, inspectable
  quantity, and no single study's batch structure can silently dominate. This is what we did.

The choice matters most precisely when studies disagree - and ours do (section 5 below).

---

## 2. The effect size and its standard error

Meta-analysis combines **effect sizes**, not p-values-turned-back-into-effects. For each gene in
each study we take the two numbers `limma` already gives us:

- the log2 fold-change, `yi` (lesional minus healthy, on the log2-CPM scale), and
- its standard error, `sei` (from the moderated model: `sqrt(s2.post) * stdev.unscaled`).

Two properties make log2FC an ideal effect measure for pooling. It is **on a common scale** across
studies (log2-CPM is dimensionless and TMM-normalised), and its **standard error is directly
available and interpretable** - a gene measured in a big, clean study has a small SE and will be
weighted heavily; a gene from a 4-vs-4 study has a large SE and is weighted down automatically.
The whole edifice below is just a principled way of taking a **precision-weighted average** of
the `yi`, where "precision" is `1/SE^2`.

---

## 3. Fixed effects vs random effects

Suppose k studies report effects `y1 ... yk` with variances `v1 ... vk` for one gene. Two models:

**Fixed-effect model.** Assume all studies estimate *one* true effect theta, and the only reason
they differ is sampling noise:

```
    yi = theta + ei,     ei ~ N(0, vi)
```

The optimal (minimum-variance) estimate is the inverse-variance weighted mean:

```
    wi = 1/vi
    theta_hat = sum(wi * yi) / sum(wi)
    SE(theta_hat) = sqrt(1 / sum(wi))
```

This is correct only if the studies really share one true effect. If they don't, it is
overconfident - the CI is too narrow because it ignores genuine between-study variation.

**Random-effects model.** Allow the true effect itself to vary across studies:

```
    yi = theta_i + ei,     theta_i ~ N(mu, tau^2),     ei ~ N(0, vi)
```

Now there are two variance components: the within-study `vi` (sampling noise) and the
between-study `tau^2` (real heterogeneity in the underlying effect - different patient mixes,
protocols, disease severity). The pooled estimate becomes an inverse-variance mean using the
*inflated* weights:

```
    wi* = 1 / (vi + tau^2)
    mu_hat = sum(wi* * yi) / sum(wi*)
    SE(mu_hat) = sqrt(1 / sum(wi*))
```

When `tau^2 = 0` this collapses back to the fixed-effect model. When `tau^2 > 0`, every weight
shrinks toward equality (a huge study no longer dominates a small one as completely), and the CI
widens to reflect that we are now generalising to a *population of studies*, not describing one.
For a transcriptome meta-analysis across independent cohorts, random effects is the honest
default, and it is what we used throughout section 11 [1].

---

## 4. Estimating tau^2: the DerSimonian-Laird estimator

The one remaining quantity is `tau^2`. We used the classic **DerSimonian-Laird (DL)** moment
estimator [1], which is closed-form (no iteration) and therefore trivially vectorisable across
24,000 genes at once. The recipe:

```
    # 1. fixed-effect weights and mean
    wi   = 1/vi
    ybar = sum(wi*yi)/sum(wi)

    # 2. Cochran's Q: observed weighted heterogeneity
    Q    = sum(wi * (yi - ybar)^2)

    # 3. method-of-moments tau^2 (floored at 0)
    C     = sum(wi) - sum(wi^2)/sum(wi)
    tau2  = max(0, (Q - (k - 1)) / C)
```

`Q` compares the observed spread of the `yi` to what sampling noise alone predicts (its
expectation under homogeneity is `k - 1`). If `Q <= k - 1`, there is no evidence of excess
heterogeneity and `tau^2` is floored at 0 (the model becomes fixed-effect for that gene). If `Q`
is large, `tau^2` is positive and the weights inflate.

We validated our vectorised DL implementation against `metafor::rma(..., method = "DL")` [2] on
STAT3 and several other genes: the pooled estimate, SE, and I-squared agreed to four decimal
places. The vectorised version exists only for speed (one pass over a genes-by-studies matrix
instead of 24,000 separate `rma` calls); it is numerically the same estimator.

---

## 5. Reading heterogeneity: Q, I-squared, tau-squared

Three related numbers describe how much the studies disagree, and they answer different questions:

- **tau^2** is the estimated between-study variance of the true effect, in the *squared units of
  the effect* (here, squared log2FC). It is the raw ingredient but hard to interpret directly.
- **Cochran's Q** is the weighted sum of squared deviations; its p-value (against a chi-squared
  with k-1 df) tests the null "all studies share one true effect." With few studies Q has low
  power, so a non-significant Q does not prove homogeneity.
- **I-squared = max(0, (Q - (k-1)) / Q) x 100%** is the share of total variation that is due to
  real heterogeneity rather than sampling noise. It is **scale-free and comparable across genes**,
  which is why we report it. Rough convention: < 25% low, 25-50% moderate, 50-75% substantial,
  > 75% considerable heterogeneity.

The key interpretive point, repeated in the white paper because it is so easy to get wrong:
**high I-squared does not mean the studies disagree about direction.** For STAT3, I-squared = 91%
- but all three studies agree STAT3 is *up*; they disagree about *how much* (+1.16, +1.63, +0.42).
Random effects handles this exactly right: it keeps the positive pooled estimate (+1.25) but rests
it on a wide base ([0.81, 1.69]) to reflect the magnitude uncertainty. A forest plot (Figure 25)
is the honest way to show this - the reader sees the three squares, their spread, and the diamond,
and can judge for themselves.

One refinement to the interval, not the estimate. The DL confidence interval above uses a normal
quantile and treats tau-squared as known. At k = 3 with I-squared = 91% that under-covers, because
tau-squared is itself very uncertain [9]. The **Hartung-Knapp-Sidik-Jonkman (HKSJ)** modification
replaces the normal quantile with a t-distribution on k-1 df and propagates the tau-squared
uncertainty [9]; it is the recommended primary method in this regime, with DL reported as a
sensitivity analysis [10]. Re-fitting STAT3 (Figure 31A): the estimate barely moves (+1.21) but the
interval roughly triples and crosses zero - REML+HKSJ +1.21 [-0.11, 2.53], p = 0.059, versus DL
+1.25 [0.81, 1.69], p = 3e-08. The width is driven by the uninformative n = 8 cohort (SE approx
0.48); the two well-powered cohorts each place STAT3 clearly up. This tempers the *magnitude* claim,
not the *direction* claim: leave-one-out re-pooling (Figure 31B) stays positive in every case
(+0.95 to +1.39), and with the large anchor cohort removed STAT3 is still up at +1.12 - so the
up-regulation is genuinely cross-cohort, not anchor-driven. The defensible summary: STAT3 is
reproducibly up across cohorts, with a pooled magnitude of approximately +1.2 that few-study
inference correctly reports as imprecise. See white paper section 11.5a.

Across all meta-significant PP-vs-NN genes, the I-squared distribution is **bimodal** (Figure
24B): a spike near 0 (524 genes with near-identical effects across cohorts - the truly
reproducible core) and a large mass above 75% (1,279 genes where direction agrees but magnitude
varies). That shape is itself informative: psoriasis has a hard reproducible core plus a
larger penumbra of context-dependent-magnitude genes.

---

## 6. Combining p-values: the Stouffer weighted-Z method

For the pathway and TF-activity meta-analyses (sections 11.5) we did not have a clean common
effect size per study - CAMERA and decoupleR return per-study statistics on different internal
scales. Here the appropriate tool is **p-value combination**, and we used **Stouffer's
weighted-Z** method:

```
    zi   = qnorm(1 - pi/2) * sign(effect_i)     # signed per-study z
    Z    = sum(sqrt(ni) * zi) / sqrt(sum(ni))   # sample-size-weighted combination
    p    = 2 * (1 - pnorm(abs(Z)))
```

Two design choices worth recording:

- **Weighting by sqrt(sample size).** An unweighted Stouffer treats a 4-sample study equally with
  a 178-sample study, which throws away information. Weighting by `sqrt(ni)` is the standard
  optimal weight when the per-study test statistics are on a common null scale [3].
- **Stouffer over Fisher.** Fisher's method (`-2 sum log(pi)`) combines p-values but is *unsigned*
  - it cannot distinguish "all studies up" from "two up, one down," which for a directional
  biological question is a serious defect. Stouffer's signed-Z preserves direction, so a pathway
  up in two studies and down in one correctly partly cancels rather than piling up as spurious
  significance. Whitlock [3] shows the weighted-Z is also generally more powerful than Fisher for
  the combined-effect question.

A numerical gotcha we hit and fixed: when a per-study p-value underflows to 0 (as CAMERA's can for
a strongly enriched set), `qnorm(1)` returns `Inf`, and a single `Inf` makes the whole Stouffer
sum `Inf`, destroying the ranking. The fix is to **cap each per-study z at a finite maximum**
(we used 38, corresponding to p ~ 1e-300) before summing. This is a display/ranking convenience,
not a change to any conclusion - every affected pathway was already unambiguously significant.

---

## 7. What the meta-analysis actually bought (and what it did not)

It is tempting to expect "more studies -> more significant genes." That expectation is wrong here,
and understanding why is instructive.

**Raw detection count went down, not up.** The anchor alone calls 18,771 genes at FDR < 0.05; the
random-effects pool calls 12,036. This is not a bug. Our anchor is unusually large (n = 178) and
clean, so it already has near-ceiling power. Adding smaller, noisier, or partially-disagreeing
cohorts and then *widening* the standard errors to absorb heterogeneity is a **regularisation** -
it pulls borderline single-study calls back toward the null. That is the conservative, correct
behaviour of random effects, and it is exactly what protects against single-study false positives.

**What was genuinely gained is replication, in three concrete forms:**

1. **264 newly-significant genes** at the matched strict threshold (FDR < 0.05 and |pooled logFC|
   > 1) that were not individually significant in the anchor - genes whose signal was too weak in
   any one study but consistent enough across cohorts to survive pooling (CCR1, C1QB, PTPN7, RHOH,
   MMP25, and others). These are the meta-analysis's positive discoveries.
2. **A reproducibility label on every gene** via I-squared - we can now separate the low-I-squared
   reproducible core from the high-I-squared context-dependent genes, which no single study can do.
3. **The retirement of a false lead** - the STAT3-beta isoform switch (section 8 below), which
   looked significant in the anchor and evaporated on replication.

The honest one-line summary: meta-analysis here did not make the answer *bigger*, it made the
answer *trustworthy*.

---

## 8. Case study in honesty: the STAT3 isoform switch that did not replicate

The single-study spliceosome analysis found lesional skin shifting toward the truncated,
dominant-negative STAT3-beta isoform (PSI-beta 4.97% -> 6.14%, Wilcoxon p = 0.017). This was an
appealing mechanistic story - a dominant-negative brake on STAT3 signalling being partially
released in disease. The meta-analysis tested it directly.

We extracted the identical junction signal (shared acceptor chr17:42,317,181; alpha donor
42,316,902; beta donor 42,316,852; two-base tolerance for coordinate conventions) from every
study's junction-level RSE, computed PSI-beta = beta / (alpha + beta) per sample, and pooled the
lesional-vs-healthy shift restricting to samples with junction depth >= 20 (below which the ratio
is too noisy to trust).

Only three studies cleared the depth floor for PP-vs-NN (the shallow SRP076982 and tiny SRP016583
did not). The result:

| Study | PP-NN shift | Wilcoxon p | note |
|---|---|---|---|
| SRP035988 (anchor) | +1.17 pp | 0.017 | reproduces section 9 exactly |
| SRP165679 (deep, well-powered) | -0.07 pp | 0.70 | no shift |
| SRP126422 (4 vs 4) | +1.01 pp | 0.86 | underpowered |
| **RE pooled** | **+0.78 pp [-0.05, +1.60]** | **0.065** | **not significant** |

The interpretation must be surgical, because it is easy to overstate in either direction:

- It does **not** overturn STAT3 the gene or STAT3 the pathway - those are confirmed and
  strengthened by the same meta-analysis (section 11.5). STAT3 protein is being made *more* in
  lesional skin; that is solid.
- It **does** retire the specific claim that the disease shifts the *isoform ratio* toward beta.
  The anchor's p = 0.017 is most parsimoniously read as a study-specific effect - real in that
  cohort, or a depth/batch artifact - not a reproducible feature of psoriasis.

Why keep a negative result in the paper at all? Because the alternative - reporting only the
anchor's p = 0.017 - would be exactly the kind of single-study over-claim that the replication
crisis is made of. The value of the whole multi-study exercise is concentrated in results like
this one: it is the mechanism that distinguishes a robust target from an over-interpreted detail.
If deeper junction-level data (targeted long-read or high-depth cohorts) later revisit the isoform
question, this is the honest prior they should start from.

---

## 9. Cross-study co-expression: what clust adds over single-study clustering

Section 10 clustered samples within one study and deferred **clust** [4] to the meta-analysis.
The reason is that clust's distinctive capability is finding modules co-expressed *consistently
across multiple datasets* - it optimises a consensus over per-dataset co-expression, so a module
it returns is one whose members move together in *every* input study, not just on average. That is
a far stronger, and far more disease-attributable, statement than a module from any single matrix.

Practical notes from running it:

- **clust is conservative by construction.** Requiring consensus across all five studies (including
  the 6- and 12-sample cohorts, whose co-expression estimates are unreliable) returned **zero**
  modules. Restricting to the three well-powered studies and loosening the tightness parameter
  (`-t 0`) returned exactly **one** high-confidence module of 65 genes. This is not a failure to
  tune - it is clust honestly reporting that only one coherent programme replicates across cohorts
  at that stringency.
- **The module is biologically legible.** All 65 genes are up in lesional skin (median logFC
  +3.43), fusing an IL-17 antimicrobial/keratinocyte programme (S100A7/8/9, DEFB103B, PI3, IL36A/G,
  CCL20, SERPINB3/4, SPRR2 and LCE3 families, KYNU, TGM1) with a mitotic cassette (BIRC5, CCNB1,
  CEP55, PLK1, PTTG1, RRM2). These are the two hallmarks of the psoriatic plaque - immune-driven
  antimicrobial response and keratinocyte hyperproliferation - recovered here as a single
  co-regulated unit purely from cross-study co-expression, with no disease labels used.
- **The gradient is the discovery.** The module's per-sample score orders healthy < peri-lesional
  < lesional, and in the one cohort with all three arms (SRP165679) peri-lesional sits cleanly
  between (-0.73, -0.33, +1.31). This is direct molecular evidence that clinically uninvolved skin
  in a psoriasis patient is a partially-activated intermediate, not simply healthy skin.

A tooling caveat for reproducibility: clust 1.18 predates NumPy 2.0 and calls the removed
`numpy.in1d`; it must be run in an environment pinned to NumPy < 2 (we used a dedicated conda env
with numpy 1.26). Its runtime also scales poorly with gene count, so we restricted the input to
the 2,380 meta-significant genes rather than the full transcriptome - which is also the more
interpretable choice, since we want disease-relevant modules.

---

## 10. Design decisions, recorded

For future extension of this meta-analysis, the choices that shaped it:

1. **Whole-skin bulk only.** Sorted-cell, blood, and explant studies were excluded because
   between-compartment variance dwarfs the lesional-vs-healthy signal.
2. **Gencode v26 throughout.** recount3's uniform annotation means the gene space is harmonised by
   construction - no cross-platform probe mapping, the historic bane of expression meta-analysis.
3. **Three contrasts (PP-NN, PN-NN, PN-PP)** to dissect the healthy -> peri-lesional -> lesional
   gradient rather than only the endpoint difference.
4. **Random effects, DL, validated against metafor.** Conservative and honest about heterogeneity.
5. **Depth floors** (junction depth >= 20 for PSI-beta; `filterByExpr` per study for genes) so
   shallow data down-weights itself rather than injecting noise.
6. **Two-stage scope.** Treatment/timepoint studies (ERP110816, SRP065812; ~230 further samples)
   were noted but deferred - their baseline/pre-treatment arms could extend the healthy-vs-lesional
   pool in a second pass without disturbing the clean cross-sectional design used here.

---

## References

[1] DerSimonian R, Laird N. **Meta-analysis in clinical trials.** *Control Clin Trials.*
1986;7(3):177-188. doi:10.1016/0197-2456(86)90046-2.

[2] Viechtbauer W. **Conducting meta-analyses in R with the metafor package.** *J Stat Softw.*
2010;36(3):1-48. doi:10.18637/jss.v036.i03.

[3] Whitlock MC. **Combining probability from independent tests: the weighted Z-method is superior
to Fisher's approach.** *J Evol Biol.* 2005;18(5):1368-1373. doi:10.1111/j.1420-9101.2005.00917.x.

[4] Abu-Jamous B, Kelly S. **Clust: automatic extraction of optimal co-expressed gene clusters from
gene expression data.** *Genome Biology.* 2018;19:172. doi:10.1186/s13059-018-1536-8.

[5] Higgins JPT, Thompson SG. **Quantifying heterogeneity in a meta-analysis.** *Stat Med.*
2002;21(11):1539-1558. doi:10.1002/sim.1186. *(I-squared statistic.)*

[6] Wu D, Smyth GK. **Camera: a competitive gene set test accounting for inter-gene correlation.**
*Nucleic Acids Res.* 2012;40(17):e133. doi:10.1093/nar/gks461.

[7] Badia-i-Mompel P, Vélez Santiago J, Braunger J, et al. **decoupleR: ensemble of computational
methods to infer biological activities from omics data.** *Bioinformatics Advances.* 2022;2(1):vbac016.
doi:10.1093/bioadv/vbac016.

[8] Wilks S, Charrad M, et al. **recount3: summaries and queries for large-scale RNA-seq expression
and splicing.** *Genome Biology.* 2021;22:323. doi:10.1186/s13059-021-02533-6.

[9] IntHout J, Ioannidis JPA, Borm GF. **The Hartung-Knapp-Sidik-Jonkman method for random effects
meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird
method.** *BMC Med Res Methodol.* 2014;14:25. doi:10.1186/1471-2288-14-25.

[10] Jackson D, Law M, Rücker G, Schwarzer G. **The Hartung-Knapp modification for random-effects
meta-analysis: A useful refinement but are there any residual concerns?** *Stat Med.*
2017;36(25):3923-3934. doi:10.1002/sim.7411.
