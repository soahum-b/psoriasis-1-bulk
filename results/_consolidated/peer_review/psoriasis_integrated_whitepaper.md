# The Psoriasis Transcriptome: A STAT3-Centred Meta-Analysis and a Molecular Staging Axis

*Integrated white paper — living document. Written alongside the analysis so that every
methodological choice is explained in plain language and can be defended in a talk. This
document integrates two analyses of the same uniformly reprocessed recount3 data:
(i) a five-cohort meta-analysis asking which genes, pathways, and regulators are
reproducibly deranged in lesional psoriasis and whether STAT3 sits at the network centre
(§0–§11), and (ii) a molecular staging axis that reframes the same data as an ordered
healthy → peri-lesional → lesional trajectory and asks in what order the disease programme
assembles (§12). A unified synthesis closes the paper (§13).*

---

## Abstract

Psoriasis carries a strong, reproducible molecular signature, which makes it an unusually
tractable disease for transcriptomics. We analyse five uniformly reprocessed recount3 cohorts
(548 samples: 341 lesional, 82 peri-lesional, 125 healthy; Gencode v26) to answer two linked
questions from one harmonised dataset.

**Part I — a STAT3-centred meta-analysis (§0–§11).** Across independent cohorts, the lesional
programme converges on an interferon / IL6-JAK-STAT3 / IL-17 inflammatory core. STAT3 is
**up-regulated in every cohort**, and the direction of that effect is robust to leave-one-out,
including removal of the large anchor cohort; its pooled magnitude is approximately +1.2 log2FC
but is imprecisely estimated at three studies, which few-study-appropriate inference
(Hartung-Knapp-Sidik-Jonkman) makes explicit (§11.5a). STAT3 also ranks in the top 3% of ~730
transcription factors by inferred regulon activity in both well-powered cohorts (pooled
Stouffer Z = 10.1), and IL6-JAK-STAT3 is a reproducibly enriched Hallmark programme. The
single-cohort STAT3-β isoform switch does **not** replicate (pooled p = 0.065), narrowing the
STAT3 hypothesis to expression/activity rather than isoform choice.

**Part II — a molecular staging axis (§12).** Using the deepest three-group cohort (SRP165679;
NN=38, PN=27, PP=28), we show that the psoriasis programme is overwhelmingly monotonic along
the healthy → peri-lesional → lesional axis: 85% of lesional DE genes are strictly monotonic,
99% show a clear monotonic rank trend, and peri-lesional skin sits a median 16% of the way to
lesional — a genuine early molecular stage. Ordering genes by how much of their lesional change
has already occurred at the peri-lesional stage separates an **early inflammatory/interferon
programme** from a **late keratinocyte-proliferation programme**: inflammation precedes
proliferation. The canonical STAT3/IL-17 hubs validate the axis, and STAT3's inferred *activity*
reaches ~30% of its lesional level at peri-lesional while its *transcript* reaches only ~10% —
consistent with activation-before-transcription, placing STAT3 activation early on the disease
trajectory.

Together the two analyses describe a reproducible, STAT3-associated inflammatory core that both
replicates across cohorts and initiates the disease sequence, with keratinocyte proliferation as
a downstream, lesion-stage consequence. Every number and figure is reproducible from public data
via a single pipeline.

---

## 0. Reading guide

This document grows one section per analysis step. Each section states, in order:
**(i)** the plain-language *goal*, **(ii)** the *decision* we made and the alternatives we
rejected, **(iii)** the *why* (with a citation where a choice is a known best practice),
and **(iv)** the *result* on our data. Terms are defined the first time they appear.

---

## 1. Motivation and question

**Clinical setting.** Psoriasis is a common, immune-mediated inflammatory skin disease.
Its plaques carry a strong, reproducible molecular signature, which makes it an unusually
good disease for transcriptomics — the study of which genes are switched on or off by
measuring their RNA.

**The question.** *Which genes and biological pathways are consistently deranged in
psoriatic lesional skin compared with healthy skin, and does the transcription factor
STAT3 — including its alternatively spliced isoforms — sit at the centre of that
signature?*

- **Transcription factor:** a protein that turns other genes on or off.
- **STAT3:** a transcription factor downstream of the IL-23/IL-17 inflammatory axis that
  drives psoriasis; a mechanistically motivated candidate.
- **Isoform / alternative splicing:** one gene can be assembled into more than one protein
  version ("isoform") by including or excluding pieces of its RNA ("splicing"). STAT3 has
  two well-known isoforms — a full-length **STAT3α** (activating) and a truncated
  **STAT3β** — that can have opposing effects, so *which* isoform dominates matters.

**Design choice (this project).** We compare **lesional psoriatic skin vs. healthy-control
skin**. A third common category, *uninvolved* skin (normal-appearing skin from a psoriasis
patient), is deliberately set aside for now because it shares the patient's genetic
background and is a different biological contrast.

---

## 2. Data source: why recount3

**Goal.** Obtain gene expression measurements that are comparable across independent
studies, so that a *meta-analysis* (combining several studies for robustness) is possible.

**Decision.** Use **recount3**, a public resource in which thousands of human RNA-seq
studies were reprocessed through a *single, uniform* pipeline. Uniform processing is what
makes cross-study comparison legitimate — otherwise differences between studies could be
artifacts of different software rather than biology.

**Key structures (defined once):**

- **Count matrix:** a table with **genes in rows, samples in columns**; each cell is a
  measure of how much that gene was expressed in that sample.
- **Metadata / sample sheet:** a companion table labelling each sample (e.g. lesional vs.
  normal). Counts tell us *how much*; metadata tells the statistical test *which group each
  sample belongs to*. Both are required.
- **RSE (`RangedSummarizedExperiment`):** a single container bundling counts + metadata +
  gene annotation so they cannot fall out of alignment.
- **Coverage → counts:** recount3 does not store tidy read counts; it stores base-level
  *coverage*. We convert coverage to estimated read counts with `transform_counts()`
  before any analysis. (Gene- and exon-level assays need this; junction-level assays are
  already raw counts and must **not** be re-transformed.)

**Why this enables the STAT3 isoform aim specifically.** recount3 provides not only
gene-level counts but also **exon-level** and **exon–exon junction** counts from the same
samples. Junctions are the direct fingerprint of splicing, so they are what let us ask
which STAT3 isoform is favoured — an analysis impossible from gene-level counts alone.

---

## 3. The starting dataset

**Study.** SRP035988 (a psoriasis RNA-seq study in recount3), loaded at the gene level.

**Annotation.** Quantified against **Gencode v26** (recount3's "**G026**" gene sums; genome
build **GRCh38/hg38**) — verified directly from the loaded object: `metadata(rse)$annotation`
returns `"gencode_v26"`, gene IDs are versioned Ensembl/Gencode identifiers (e.g.
`ENSG00000278704.1`), and the counts file is `human.gene_sums.G026.gtf.gz`. This is recount3's
uniform human annotation — the *same* gene model used for its GTEx and TCGA data — which is
what will let us pool studies cleanly in the later meta-analysis (every study is counted
against one identical gene definition).

**Dimensions as loaded:** **63,856 genes × 178 samples**.

**Group assignment.** Sample labels in the metadata are two clean strings — `normal_skin`
and `Psoriasis_skin` (`lesional psoriatic skin`) — which we map to:

| Group | Meaning | n |
|-------|---------|---|
| **NN** | normal / healthy-control skin | **83** |
| **PP** | lesional psoriatic skin | **95** |
| unassigned (`NA`) | — | **0** |

**Verification (a habit, not an afterthought).** Whenever groups are assigned by matching
text in the metadata, the first check is a count of how many samples landed in each group
and how many landed *nowhere*. Here every sample was labelled and nothing was silently
dropped, and this study contains no *uninvolved* samples — so the lesional-vs-healthy
contrast maps cleanly onto it. When we scale to additional studies, this same table is the
first thing we will inspect, because other datasets use messier labels.

---

## 4. Gene filtering: keep the genes worth testing

**Goal.** Of 63,856 measured genes, most are near-silent in skin (essentially zero counts
everywhere). Testing them adds noise and worsens the *multiple-testing burden* (the more
genes we test, the harder we must correct p-values, which costs power to detect real
signal). We therefore remove low-expression genes *before* testing.

**The intuition we started from.** Keep a gene only if it is expressed in at least *some
minimum number of samples* — not just one (a single-sample spike is usually a technical
artifact) — and tie that minimum to the size of the **smaller group**, so that a gene
switched on in only one condition is still retained.

**The important correction from the literature.** That intuition is right about *how many*
samples, but wrong if the count is taken **within a specific group**. The requirement for
valid filtering is that the filter be **independent of the statistic used to test for
differential expression**; a two-stage "filter then test" approach increases power *only*
when the filter criterion does not use the group labels [1]. A filter such as "expressed
in ≥50% of the *psoriasis* samples" uses the group labels and is therefore *not*
independent — it can inflate the false-discovery rate. edgeR's author, Gordon Smyth, makes
this point explicitly: group-oriented filtering is statistically invalid, whether applied
before or after the differential-expression analysis, because it inflates the error rate
[2, forum guidance consistent with 1].

**The professional rule (what we adopt).** Use edgeR's **`filterByExpr()`** [3]. It blends
both ideas correctly:

- it decides *how many* samples a gene must be expressed in from the **minimum group
  size** (so single-condition genes are protected), but
- it counts those samples **across the whole experiment**, never within one group (so the
  filter stays independent of the test).

Its defaults are `min.count = 10`, `min.total.count = 15`, `large.n = 10`,
`min.prop = 0.7`. This replaces the earlier hard-coded "expressed in ≥ 60 samples" rule,
which was neither tied to the design nor justified.

- **CPM (counts per million):** expression rescaled so that samples sequenced to different
  depths become comparable; "≥ 1 CPM" is a common floor for "meaningfully expressed."

**Takeaway sentence for a talk.** *"50% of the smaller group" is the right instinct for the
sample count, but you must apply it experiment-wide, never gated to one condition —
otherwise you bias which genes are even allowed into the test.*

**Result on our data.** Starting from 56,937 uniquely-named genes (after collapsing
duplicate gene symbols):

| Rule | Genes kept |
|------|-----------|
| `filterByExpr()` (adopted) | **24,528** |
| old hard-coded ≥ 1 CPM in ≥ 60 samples | 19,784 |
| "≥ 1 CPM in ≥ 42 samples" (50%-of-min-group intuition, applied experiment-wide) | 20,350 |

The decisive comparison: **every gene retained by the old ≥ 60 rule is also retained by
`filterByExpr` (0 lost), while `filterByExpr` keeps 4,744 additional genes.** Those extra
genes are legitimately testable transcripts that the arbitrary flat cutoff discarded.
`filterByExpr` is more permissive precisely where it is statistically safe to be, because
it derives the sample-count threshold from the library sizes and the ~83/95 group
structure rather than from a hand-picked number. **We carry 24,528 genes forward.**

---

## 5. Normalization (TMM): making samples fairly comparable

**Goal.** Even after filtering, raw counts cannot be compared *between* samples. RNA-seq
counts are a **fixed budget**: the counts in a sample sum to its library size, so if a few
genes are enormously expressed in one condition, they consume a large share of the reads
and every *other* gene's count is mechanically deflated — even genes whose true expression
did not change.

**Why this bites psoriasis specifically.** Psoriatic plaques massively over-express a
handful of skin genes (keratins, S100 proteins). Without correction, an unchanged
"housekeeping" gene would appear **artificially lower** in the lesional samples purely
because those few genes ate the read budget — a composition artifact, not biology.

**Decision.** Apply **TMM normalization** (Trimmed Mean of M-values; `calcNormFactors(...,
method = "TMM")`) [4]. TMM estimates, for each sample, a single scaling factor that removes
this compositional distortion, on the assumption that most genes are *not* differentially
expressed.

**How TMM works (plain mechanism).** Pick a reference sample. For every gene, compute the
**M-value** — the log2 ratio of its expression in the test sample vs. the reference,
`M = log2(test / reference)`. If there were no composition bias, *most* genes should be
unchanged, i.e. M ≈ 0. But genuinely up- and down-regulated genes have large positive and
negative M-values that would corrupt a naive average. TMM therefore:

1. **Trims the extremes** — discards the top and bottom 30% of genes by M-value (the likely
   differentially-expressed ones) and the most extreme 5% by overall expression level, then
2. takes a **weighted mean of the surviving M-values** (weighting by precision, so
   well-measured genes count more).

That trimmed mean *is* the sample's normalization factor. The load-bearing assumption is
that **the majority of genes are not differentially expressed**, so the "typical" log-ratio
should be zero — any systematic departure is composition bias, and TMM scales it away.
Because it trims rather than averages blindly, TMM is robust even when a few genes (the
keratins/S100s here) dominate the read budget.

**Why TMM and not the alternatives.** The choice of normalization must match the *comparison*
being made. Differential expression is a **between-sample** comparison of the *same gene*
across samples, which rules out the within-sample abundance units and the depth-only scalings:

| Method | Corrects depth | Corrects gene length | Corrects **composition** | Right for our DE? |
|---|:---:|:---:|:---:|---|
| RPKM / FPKM | ✓ | ✓ | ✗ (worsens) | ✗ — within-sample unit; RPKM even sums differently per sample [5, 6] |
| TPM | ✓ | ✓ | ✗ (fixed per-sample sum *builds in* the bias) | ✗ — within-sample abundance unit [5, 9] |
| CPM / RPM (total-count) | ✓ | ✗ | ✗ (total is inflated by the dominant genes) | ✗ for inference; ✓ only for visualization [6] |
| **TMM** (adopted) | ✓ | ✗ | **✓** | **✓ — composition-aware, native to edgeR/limma-voom** [4] |
| DESeq2 median-of-ratios | ✓ | ✗ | **✓** | ✓ — validated *sibling*; agrees closely with TMM [7, 8, 10] |

The reasoning in one line each: **gene-length correction cancels in DE** (the same gene has
the same length in both groups), so length-normalized units (TPM/RPKM) add noise and, worse,
TPM's constraint that every sample sum to 1,000,000 mechanically deflates all other genes when
the psoriatic keratins/S100s dominate — reintroducing the exact bias we must remove.
**Total-count CPM corrects only depth**, and the total is itself inflated by those dominant
genes, so it too fails on composition (we use log2-CPM for visualization only). **TMM and
DESeq2's median-of-ratios** are the two validated, composition-aware, between-sample methods;
they rest on the same "most genes unchanged" assumption and are benchmarked to agree [10]. We
use **TMM** because it is the native normalization of our edgeR/limma-voom framework. *(A
standalone deep-dive with worked definitions and the full reference set is in
`normalization_deep_dive.md`.)*

**Result on our data.**

- TMM scaling factors ranged **0.565 – 1.279** across the 178 samples. A factor of 0.565
  means that sample's effective library is scaled down to ~57% (a few genes had inflated
  its apparent depth); a factor of 1.28 scales another sample up. That spread is precisely
  the between-sample distortion TMM removes.
- **Figure 1** (`fig1_tmm_before_after.png`): per-sample log2-CPM distributions before and
  after TMM. Before, sample medians wobble; after, they align into a flat line across all
  178 samples — the visual signature of successful normalization.

![Figure 1. Per-sample log2-CPM distributions before and after TMM normalization.]({{artifact:art_fadf8e20-4342-4609-99d3-34e5fd47c132}})

- **Figure 2** (`fig2_pca.png`): PCA (Principal Component Analysis — a method that finds the
  axes of greatest variation in the data) of the normalized log2-CPM. **PC1 captures 32.3%**
  of all variation and cleanly separates lesional (PP) from normal (NN) skin.
  - **Caveat noted for follow-up:** one NN sample sits inside the PP cloud (and one PP among
    the NN) — a possible mislabel or biological outlier, to be revisited when we scale up.

![Figure 2. PCA of normalized log2-CPM; PC1 (32.3%) separates lesional from normal skin.]({{artifact:art_26deffce-f432-4f5f-951e-f4cfebb1d486}})


### 5.1 PCA diagnostics: does PC1 correspond to a real variable, and what are the eigenvectors?

PCA is worth interrogating rather than eyeballing. Two questions: (a) does any principal
component *correspond to the disease variable*, and (b) what do the **eigenvectors**
(loadings) tell us about which genes build each axis?

- **What PCA is, precisely.** PCA re-expresses the samples along new axes (principal
  components) that are ranked by how much variance they capture. Each PC has an
  **eigenvalue** (the amount of variance it carries — plotted in the scree plot) and an
  **eigenvector** (the "loadings": a weight per gene saying how much each gene contributes
  to that axis).

- **(a) PC1 is the disease axis — quantified.** We tested each PC's association with the
  lesional/normal label:

  | PC | Variance | Correlation with group | t-test p |
  |----|---------|------------------------|----------|
  | PC1 | 32.3% | **0.95** | 3×10⁻⁷⁴ |
  | PC2 | 11.7% | 0.05 | 0.48 (n.s.) |
  | PC3 | 4.7% | −0.14 | 0.06 (n.s.) |
  | PC4 | 4.2% | 0.07 | 0.38 (n.s.) |

  A linear model shows **90% of the variation along PC1 is explained by the group label**.
  So PC1 is not merely "some structure" — it is almost entirely the psoriasis-vs-normal
  contrast. Critically, **PC2–PC6 do *not* associate with disease** (all non-significant),
  which means the disease effect is concentrated in one dominant axis and the remaining
  variation reflects other, unlabelled sources (the outlier, inter-individual differences,
  or possible batch structure) — candidates to model as covariates later.

- **Scree plot (Figure 3, `fig3_scree.png`).** PC1's bar towers over the rest (32.3% vs.
  11.7% for PC2, then a long tail), the classic "elbow" showing a single dominant axis.

![Figure 3. Scree plot: variance explained per principal component.]({{artifact:art_7314744e-1f35-4742-86bb-e38d357693a4}})


- **(b) Eigenvector interpretation (Figure 4, `fig4_pc1_loadings.png`) — the payoff.** The
  genes with the largest PC1 loadings *are* the canonical psoriasis signature, recovered
  with **no knowledge of the labels**:
  - **Positive end (up in lesional):** DEFB4A/DEFB4B (β-defensins), S100A7A/S100A8/S100A9/
    S100A12 (S100 antimicrobial proteins), SERPINB4, PI3 (elafin), SPRR2A/2B/2F and LCE3A
    (cornified-envelope genes), IL36A — the antimicrobial-peptide and epidermal-barrier
    programs that define the psoriatic plaque.
  - **Negative end (up in normal):** AADACL3, DGAT2L6, AWAT2, ELOVL3, THRSP — **lipid and
    sebaceous-metabolism** genes, the normal-skin programs suppressed in lesions.

![Figure 4. PC1 gene loadings; top positive and negative contributors to the disease axis.]({{artifact:art_433e0c98-c526-4224-a53c-dfc8d927e65a}})

  **Interpretation for a talk:** the *unsupervised* structure of the data reconstructs the
  disease biology before any differential-expression test is run. This is strong evidence
  that the signal is real and dominant, not an artifact of the statistical model — the data
  "already knows" which samples are psoriatic.

**Should we evaluate eigenvectors as a matter of routine? Yes.** Inspecting PC loadings is a
cheap, powerful sanity check: if PC1 had instead been driven by, say, ribosomal or
mitochondrial genes, or by a single sample, that would flag a technical artifact rather than
biology. Here the loadings confirm biology, which justifies proceeding to formal DE testing.

> **Technical footnote — the linear algebra of our PCA (covariance, SVD, rank, and the
> pseudoinverse).** PCA is the eigendecomposition of the gene–gene **covariance matrix**
> **C**. We used the *covariance* form (`prcomp(scale.=FALSE)`), not the *correlation* form
> (`scale.=TRUE`): because all genes are already on a common log2-CPM scale, standardizing to
> unit variance would inflate the weight of near-noise, low-variance genes and let them
> compete with the biologically meaningful high-variance ones. Because **C** is symmetric and
> positive semi-definite, the spectral theorem guarantees real, non-negative eigenvalues (the
> per-PC variances) and **orthogonal** eigenvectors — so PCA is a rigid *orthogonal
> transformation* (a rotation) **Y = X V** with **VᵀV = I**, rotating from 24,528 correlated
> gene axes to 177 mutually **uncorrelated** PC axes. This orthogonality is what lets us test
> each PC's disease-association independently.
>
> In practice `prcomp` computes this via the **singular value decomposition (SVD)** of the
> centered data matrix, **X = U D Vᵀ**, rather than by ever forming the 24,528 × 24,528
> covariance matrix: **V** holds the gene loadings, **U** the sample scores, and each singular
> value obeys `eigenvalue = (singular value)² / (n − 1)` — verified on our data to a precision
> of 10⁻¹³. SVD is preferred because it never materializes the 600-million-entry covariance
> matrix and avoids the numerical error of squaring the condition number (the reason
> SVD-based `prcomp` is favored over covariance-based `princomp`).
>
> With **p ≫ n** (24,528 genes, 178 samples), centering leaves **C** with **rank n − 1 =
> 177**: the first 177 PCs carry 100% of the variance and the remaining directions form a
> **nullspace of dimension 24,528 − 177 = 24,351** — gene-combinations along which the samples
> show no variation at all and whose eigenvectors are arbitrary (hence meaningless to
> interpret). This degeneracy is why only the leading handful of PCs are examined, and it is
> the concrete face of the "small-n, large-p" problem in genomics. It also makes **C**
> **singular** (smallest singular value ≈ 3 × 10⁻¹²): any method that would invert it must use
> the **Moore-Penrose pseudoinverse** **X⁺ = V D⁺ Uᵀ**, where **D⁺** inverts the non-zero
> singular values and leaves the nullspace zeros untouched, yielding the least-squares,
> minimum-norm solution. We do not invert **C** for PCA itself (a forward decomposition), but
> the same rank-deficiency is exactly why the downstream differential-expression step avoids a
> joint 24,528-gene covariance inversion and instead fits one gene at a time with
> empirical-Bayes shrinkage (limma) — a statistical resolution of the same degeneracy the
> pseudoinverse resolves algebraically.

- **CPM / log2-CPM:** counts per million, log-transformed; the standard scale for
  visualizing and clustering expression because it is depth-adjusted and compresses the
  huge dynamic range of counts.

---

## 6. Differential expression: the statistical model

**Goal.** Move from the unsupervised evidence (PCA) to a *formal, per-gene* test: for each
gene, is the lesional-vs-normal difference larger than replicate noise, and by how much?

**Why not Poisson — the counts are over-dispersed.** The Poisson distribution (the baseline
count model) forces **variance = mean**, capturing only technical sampling noise. Real
replicates add **biological** variability, pushing variance above the mean — *over-dispersion*
(defined relative to the Poisson tie, not "variance is large": binomial counts are
*under*-dispersed). Measured on our 83 normal samples, the **median variance/mean ratio is
24.2** (Poisson predicts 1) and **99.5% of genes exceed 2× their mean** (Figure 5). The fix is
the **Negative Binomial = Gamma-Poisson** model, whose variance `Var = μ + φμ²` adds a
**dispersion** term φ for biological spread. We measured edgeR's common dispersion at
**φ = 0.124 (BCV = 35%)**, typical for clinical tissue.

![Figure 5. Mean-variance relationship in normal samples; genes track the negative-binomial curve, far above Poisson.]({{artifact:art_d4ff8c4e-dbba-44bb-b15a-1aeb7d0ff09c}})

**The three tools and why we use limma-voom.** edgeR and DESeq2 fit Negative-Binomial GLMs to
raw counts (testing with a quasi-likelihood F-test and a Wald test respectively). **limma-voom**
(our choice) instead transforms counts to log2-CPM, fits the mean-variance trend and turns it
into a **precision weight per observation**, then runs a **weighted linear model with a
moderated t-test** whose per-gene variance is empirical-Bayes shrunk across genes. At our sample
size (**n = 178**) the three converge; limma-voom is fast, well-calibrated, and the most
flexible for the batch/study covariates the later meta-analysis will need.

| Tool | Data | Distribution | Test |
|---|---|---|---|
| edgeR | counts | Negative Binomial GLM | quasi-likelihood F |
| DESeq2 | counts | Negative Binomial GLM | Wald |
| **limma-voom** (ours) | log2-CPM | Normal + precision weights | moderated t |

**Shrunken fold-changes.** Low-count genes can show huge but meaningless raw fold-changes.
DESeq2 shrinks these explicitly (`lfcShrink`); limma-voom protects the *inference* through
variance moderation, and our thresholds — **adjusted p (FDR) < 0.05 AND |log2FC| > 1** —
enforce that reported genes are both statistically solid and at least 2-fold. All tools control
multiple testing with **Benjamini-Hochberg FDR** across ~24,500 genes, so we always threshold on
*adjusted* p. *(Full model derivations, dispersion detail, and the reference set are in
`differential_expression_deep_dive.md`.)*

### 6.1 Results: the lesional-vs-normal signature

Running limma-voom on the 24,528 filtered genes yields **3,477 significant genes** at
**adjusted p (BH-FDR) < 0.05 AND |log2FC| > 1** — **1,511 up** and **1,966 down** in lesional
skin. The effect-size filter is doing real work: 18,773 genes pass adjusted-p alone, but only
3,477 also clear the 2-fold bar — the visual argument for never thresholding on p alone.

**The biology validates the whole pipeline.** The strongest up-regulated genes are the textbook
psoriasis signature, recovered without any prompting: **DEFB4A** (+13.2 log2FC, ~9,000-fold — an
antimicrobial beta-defensin), **S100A7A**, **SERPINB4**, **IL36A**, **SPRR2A**, **PI3**,
**S100A8/9/12**, **LCE3A** — the antimicrobial / cornified-envelope program of the psoriatic
plaque. The most down-regulated genes are barrier keratins and lipid/sebaceous programs
(**KRT77**, **BTC**, **RORC**, **IL34**).

![Figure 6. Volcano plot of lesional (PP) vs normal (NN) skin. Red = significantly up in psoriasis, blue = down, grey = not significant; dashed lines mark the adjusted-p < 0.05 and |log2FC| > 1 thresholds.]({{artifact:art_a7e7fe83-44dc-40a1-a14c-498d31566fc4}})

**Our target lit up.** STAT3 is significantly up (**+1.16 log2FC, ~2.2-fold, adj.p ~ 6e-52**),
and its entire signalling neighbourhood moves coherently with it — the activators upstream, the
canonical target gene downstream, and the IL-17 effector axis it drives:

| Gene | log2FC | adj.p (BH) | q (Storey) | Role |
|---|---|---|---|---|
| IL17A | +5.71 | 6.5e-52 | 8.7e-53 | Th17 effector cytokine |
| IL17F | +4.72 | ~0 | ~0 | Th17 effector cytokine |
| IL36G | +5.60 | ~0 | ~0 | amplification loop |
| IL23A | +2.13 | ~0 | ~0 | Th17-driving cytokine |
| SOCS3 | +2.11 | 1.2e-40 | 1.6e-41 | **direct STAT3 target gene** |
| IL6   | +2.02 | 7.0e-17 | 9.3e-18 | STAT3 activator |
| **STAT3** | **+1.16** | **5.9e-52** | **7.8e-53** | **our hub** |
| NFKB1 | +0.68 | ~0 | ~0 | NF-κB axis |
| EFTUD2 | +0.51 | ~0 | ~0 | spliceosome (isoform hypothesis) |

STAT3 up, its activators (IL6, IL23) up, its target (SOCS3) up, the IL-17 axis up: the
hypothesized circuit is transcriptionally active. One caveat to carry forward — STAT3's biology
is driven largely by *phosphorylation* and *isoform switching* (STAT3α vs the truncated
STAT3β), which this total-mRNA measurement cannot see. The +1.16 is the total-transcript signal;
the exon/junction isoform analysis is a later stage.

**Multiple-testing robustness (BH vs Storey).** We report both the BH-adjusted p (`adj.P.Val`)
and the **Storey q-value**, which estimates the true-null proportion **pi0 = 0.133** — i.e. only
~13% of the transcriptome is estimated to be truly unchanged, biologically sensible for a
massively remodelled tissue. Because q = pi0 × BH, the q-values are uniformly smaller, but after
the |log2FC| > 1 filter both methods converge (3,477 vs 3,479 genes): the headline list is
**robust to the multiple-testing method**, driven by strong, large-effect biology rather than
the FDR estimator. (Full treatment: `differential_expression_deep_dive.md` §5.)

![Figure 7. Multiple-testing diagnostic. (A) p-value distribution: a spike near zero (real signal) on a flat shelf of null genes; dashed line marks the estimated null level pi0 = 0.133. (B) Storey q-value vs BH adjusted p; points lie on/below the 1:1 line because q = pi0 x BH.]({{artifact:art_d5cd7635-7490-47a0-a3ec-dac25cfb089f}})

---

## 7. Pathway analysis: from a gene list to a coordinated program

Differential expression asks a **per-gene** question — "is STAT3 individually up?" (yes, +1.16).
But genes act in **teams** (pathways), and a pathway can be powerfully activated even when no
single member clears a fold-change threshold, if many members each shift **modestly and
consistently** in the same direction. Pathway analysis is how we recover that team-level signal.

**Gene Set Enrichment Analysis (GSEA)** is our primary method. It is **threshold-free**: rather
than first splitting genes into "significant" and "not," it ranks **all 24,528 genes** by the
limma moderated t-statistic (effect size ÷ its standard error, so a gene must be both large and
reproducible to rank high), then, for a given gene set, walks down the ranked list accumulating a
**running enrichment score** — stepping up at every set member (weighted by the statistic) and
down at every non-member. If the set's genes cluster near the top of the ranking, the score
climbs to a positive peak — the **Enrichment Score (ES)** — which is size-normalized to a
**Normalized Enrichment Score (NES)**. A **positive NES means the set is coordinately up in
lesional skin; negative means down.** Significance comes from permutation: `fgsea`'s multilevel
Monte-Carlo scheme generates a null distribution of scores, yielding a p-value that is then
BH-adjusted across sets [20,21]. We tested four reference collections — **MSigDB Hallmark** [22],
**Reactome** [23], **KEGG** [24], and **GO Biological Process** [25] — for 6,505 gene sets in
total.

**The STAT3 hypothesis holds at the pathway level.** Every pathway in the hypothesized circuit is
significantly and strongly up. Across the 50 Hallmark signatures, 31 are significant; the top of
that list is a textbook psoriasis portrait — E2F/MYC/G2M targets (hyperproliferating
keratinocytes), interferon-γ response, inflammatory response, TNFα-via-NF-κB, and
**IL6/JAK/STAT3 signaling** (NES +2.49, adj.p 6.4e-11). Reactome and KEGG resolve the specific
mechanistic sets — STAT3 nuclear events, IL-17 signaling, JAK-STAT signaling — and GO:BP supplies
the biological output: keratinization and cornified-envelope formation, the histology of a plaque.

![Figure 8. GSEA running enrichment score for the Hallmark IL6/JAK/STAT3 signaling set across all 24,528 ranked genes. Red ticks (bottom rug) mark the 81 pathway genes; because they concentrate at the up-in-psoriasis (left) end, the running score climbs to a peak of +0.635 (NES +2.49, adj.p 6.4e-11). The leading edge - the genes before the peak that drive the signal - includes STAT3 itself, alongside STAT1, SOCS3, OSMR and the IL/IFN receptor chains.]({{artifact:art_9d5174b2-a18c-43f2-aed3-0943c2495955}})

![Figure 9. GSEA pathway landscape. The top STAT3-axis and overall enriched sets across all four collections, ranked by NES. Every set points right (positive NES = up in lesional skin); dot size = gene-set size, colour = significance. The IL6/JAK/STAT3 axis, its upstream/parallel inflammatory drivers (TNFa/NF-kB, IFN-gamma, IL-17), and the downstream keratinocyte program all sit among the strongest, most significant sets.]({{artifact:art_0d2f191c-02a8-43ed-9ef0-2bce07b67691}})

**GSEA and over-representation analysis (ORA) can disagree — informatively.** ORA is the
threshold-based alternative: it fixes a DE cutoff (our 3,477 genes), then asks by **Fisher's
exact test** on a 2×2 table whether a pathway is over-represented among the hits. For **GO
Keratinization**, both methods agree emphatically (46/79 genes DE, 4.1× enriched, Fisher p =
1e-19; GSEA NES +2.54). But for **Reactome IL-17 signaling**, GSEA calls it significant (NES
+2.02, adj.p 0.006) while ORA calls it *depleted* (7/70 DE, p = 0.88). The reason is mechanistic:
only the cytokines and receptors (IL17A/F/C, IL17RE) clear the strict |log2FC| > 1 bar; the ~60
intracellular signaling genes (IRAKs, TRAFs, MAP kinases) shift coordinately but modestly and are
**discarded at the threshold by ORA but counted by GSEA**. The disagreement itself reveals *how*
the pathway is regulated. This is why we run both — agreement gives confidence, disagreement gives
insight. (Full treatment, with the permutation-scheme details and the worked 2×2 tables:
`pathway_analysis_deep_dive.md`.)

**The disagreement is systematic, not anecdotal.** To confirm the two hand-picked pathways
generalize, we ran ORA across **all** 6,505 gene sets (clusterProfiler's `enricher`, same universe
of 24,528 tested genes) and merged it with GSEA verdict-by-verdict. Restricting to the 3,584
up-enriched sets both methods evaluate, the two agree on **77%** — but the disagreement is almost
entirely one-directional: **1,204 sets are GSEA-significant yet ORA-null, versus only 64 the
reverse** (176 significant by both). ORA is essentially a conservative *subset* of GSEA: nearly
everything ORA flags, GSEA also flags, plus a large tier of coordinated-moderate programs ORA's
threshold discards. Per collection the gap is consistent — Hallmark 31 GSEA vs 14 ORA, Reactome
545 vs 42, GO:BP 1,419 vs 272, KEGG 67 vs 17.

![Figure 10. GSEA vs ORA for every up-enriched gene set tested by both methods. Each point is one gene set; axes are -log10 adjusted-p under GSEA (vertical) and ORA/Fisher (horizontal); dashed lines mark padj 0.05. The dense vertical band on the left is the 1,204 GSEA-significant / ORA-null sets - coordinated moderate signal below ORA's threshold. Keratinization sits top-right (both methods strongly agree); IL-17 signaling sits in the left band (GSEA-only), the teaching case in its systematic context.]({{artifact:art_025c8b87-fa15-45e3-8491-0fdc01198ab3}})

We also ran a **directional ORA** (restricting the hit list to up-regulated genes only) as the
sharpest possible test of the STAT3 axis. Fold-enrichment *rises* across every axis pathway —
keratinization 3.65×→5.14×, IL6/JAK/STAT3 2.07×→3.25×, keratinocyte differentiation 2.49×→3.63×,
TNFα/NF-κB 1.44×→2.12× — confirming these are specifically **up**-programs, not bidirectional
noise. The lone exception is IL-17 signaling, which stays null even directionally (1.10×, adj.p =
1.0): the teaching case survives the strictest ORA we can construct, precisely because its signal
lives in sub-threshold intracellular machinery that only a rank-based method can see.

The pathway layer thus corroborates the STAT3 hypothesis at a **third level of evidence**: the
single gene (§6.1), its immediate transcriptional neighbourhood (§6.1 table), and now the
coordinated IL6→JAK→STAT3→IL-17/NF-κB→keratinocyte **program** activated as a unit.

### 7.1 Robustness: does the STAT3 signal survive stricter gene-set tests?

Competitive GSEA (§7) has a known weakness — it treats the genes in a set as independent votes, but
pathway genes are **co-regulated**, so it can be over-confident [35,36]. Before we build a hypothesis
on the STAT3 regulon, we stress-tested it with three additional methods that plain GSEA does not
cover. (Concepts and worked statistics: `pathway_analysis_deep_dive.md` §9.)

**CAMERA and ROAST — two harder questions, same answer.** GSEA asks a *competitive* question ("more
changed than the rest of the genome?"). **CAMERA** [36] asks the same question but first estimates
the inter-gene correlation and discounts the redundant evidence; **ROAST** [37] asks a different,
*self-contained* question ("is this set engaged **at all**?") using rotation rather than permutation,
which is the natural test for a named hypothesis set. On the **STAT3 regulon** (371 measured CollecTRI
targets) all three agree the set is **up in lesional skin**, but the p-values are instructive: GSEA
7.1×10⁻¹⁴, **CAMERA 4.6×10⁻⁴**, ROAST 5.0×10⁻⁵. The ten-order-of-magnitude gap between GSEA and
CAMERA is exactly the correlation inflation the method corrects — so **we quote CAMERA's more
conservative p as the honest competitive statement**. Across the 50 Hallmark sets, GSEA and CAMERA
agree on direction and significance for 27, leaving the core immune/proliferation program intact.

![Figure 18. Enrichment robustness across three gene-set tests. (A) Hallmark GSEA NES (x) vs CAMERA -log10 FDR (y); 27/50 sets are significant in both, same direction (red), with the immune and proliferation programs - E2F, interferon-alpha/gamma, G2M, MYC, IL6/JAK/STAT3, inflammatory response - in the strong-and-certain corner. (B) The STAT3 regulon (371 targets) tested three ways: all agree it is up in lesional skin, but CAMERA (correlation-corrected, p 4.6e-4) is far less extreme than GSEA (p 7.1e-14), and ROAST (self-contained, p 5.0e-5) answers a different question. The CAMERA value is the one we report.]({{artifact:art_a1ba3a6e-b345-4945-a8e8-7e712d1ccfb7}})

**Per-sample pathway scores (ssGSEA) — from one p-value to 178 numbers.** All the tests above collapse
the cohort into a single verdict per set. **GSVA/ssGSEA** [38,39] instead scores each pathway in each
sample, converting the gene × sample matrix into a pathway × sample matrix — the object we can
correlate, cluster, and (crucially) pool across studies in the meta-analysis. Because the
Bioconductor `GSVA` package would not install in our sandbox, we implemented ssGSEA in ~20 lines of
base R (a fully transparent alternative to the black-box call) and validated it against the biology:
the per-sample IL6/JAK/STAT3 score cleanly separates the groups (NN mean 0.505 vs PP 0.636, Wilcoxon
p = 5.9×10⁻²⁸), reproducing the cohort-level GSEA verdict from an independent per-sample computation.
Two payoffs follow. First, the score **tracks STAT3 expression continuously** (Spearman ρ = 0.81
across all 178 samples) — pathway activity is a graded axis, not an on/off switch, and STAT3 sits on
it. Second, correlating the score against **PSI_β** reproduces the **Simpson's-paradox** signature we
found with EFTUD2 (§9): overall ρ = 0.09 (n.s.) but **within lesional skin ρ = −0.18** — higher
pathway drive goes with *less* β (more activating α), consistent with the β brake being regulated
separately from bulk pathway activity.

![Figure 19. Per-sample IL6/JAK/STAT3 ssGSEA score. (A) The score tracks STAT3 expression across all 178 samples (Spearman rho 0.81). (B) Score vs PSI_beta (STAT3-beta fraction): overall flat, but within lesional skin the slope is negative (rho -0.18) - the same Simpson's-paradox pattern as EFTUD2, arguing the beta brake is not simply driven by pathway activity.]({{artifact:art_bb597911-3cd9-406e-bf82-6b56f4ad5b52}})

![Figure 20. Per-sample ssGSEA pathway-score heatmap (row z-scored), samples ordered normal-to-lesional. The two-block structure is visible sample by sample: metabolic/lipid programs (myogenesis, adipogenesis, fatty-acid metabolism, OXPHOS) high in normal skin; immune and proliferation programs (IL6/JAK/STAT3, inflammatory response, interferon-alpha/gamma, allograft rejection, E2F/G2M/MYC) high in lesional. The row dendrogram splits along exactly that biology.]({{artifact:art_b69fcadb-0ca3-4a03-84d5-38033dac0022}})

**Leading-edge redundancy — six "immune pathways" are largely one signal.** GSEA reports many
significant immune sets, but their **leading edges** (the genes actually driving each score) overlap
heavily, because psoriasis has essentially one dominant immune program that many curated sets each
capture a slice of. An **UpSet plot** [40] makes this concrete: the leading edges of the six top
immune Hallmark pathways sum to 537 gene-memberships but only **373 unique genes (redundancy 1.44×)**,
and a shared core — IL6, IRF1, CXCL9/10/11, CCL2/5, IL1B, TLR2 — recurs in four or more of the six.
The lesson for interpretation: these are **not** six independent discoveries, and the STAT3/IL-17 core
is what they hold in common.

![Figure 21. Leading-edge overlap across the six top immune Hallmark pathways (UpSet). Bars are the number of genes in each specific intersection; the shared inflammatory core (IL6, IRF1, CXCL9/10/11, CCL2/5, IL1B, TLR2) recurs across most sets. Only 40-63 genes are unique to any single pathway - the six co-significant sets are largely one redundant signal, not six independent findings.]({{artifact:art_428678be-3cf6-4e27-a142-ba86639f7a9e}})

Taken together, §7.1 does not overturn anything in §7 — it **hardens** it. The STAT3 program is up
in lesional skin under the correlation-corrected competitive test (CAMERA), under the self-contained
test (ROAST), and as a graded per-sample score that tracks STAT3 itself; and the apparent multiplicity
of enriched immune pathways resolves to one shared inflammatory core.

---

## 8. Network analysis: from a coordinated program to its upstream drivers

Pathway analysis (§7) answered *which gene sets moved*. It cannot, by construction, answer the next
question — *which regulator drove the movement* — because a gene set is an unordered bag of genes
with no arrows. The largest-magnitude signature in the psoriasis transcriptome is **proliferation**
(cell-cycle/E2F/MYC sets, §7), but that is the phenotypic *output* of a plaque, not its cause. To
find the cause we move from a **list** to a **network**: nodes (genes) connected by directed,
signed edges (who regulates whom) drawn from prior biological knowledge, and we ask which node sits
*upstream* of the observed changes. (Full conceptual treatment — network types, algorithms, and
their trade-offs — in the companion `network_analysis_deep_dive.md`.)

**Transcription-factor (TF) activity inference** is the most direct network method for our question.
A TF's own mRNA level is a poor proxy for its activity — activity is set post-transcriptionally, by
phosphorylation, nuclear translocation, and dimerization. But an active TF leaves a **fingerprint**:
the genes it activates go up and the genes it represses go down. A **regulon** is the curated set of
a TF's targets, each signed +1 (activated) or −1 (repressed). We used **CollecTRI** [26], a
signed regulon covering 1,201 TFs, and **decoupleR** [27] with its recommended univariate linear
model (**ULM**): for each TF, the genome-wide limma t-statistics are regressed on that TF's signed
regulon membership, and the slope's t-value becomes the **activity score** (positive = active in
lesional skin). This infers *activity from targets*, sidestepping the mRNA-proxy problem entirely.

**STAT3 is a significantly active regulator, sitting in the predicted hub cluster.** Of 732 TFs with
a testable regulon, STAT3 scores **7.68 (adj.p = 6.3e-13), ranking 19th — the 97th percentile.** But
the more striking result is its *company*: the most-active TFs, in order, are the entire **NF-κB
family** (RELA, RELB, REL, NFKB1/2), then **STAT1**, **MYC**, **E2F1**, **STAT3**, **IRF1**. Every
hub we had inferred only *indirectly* from the enriched pathways in §7 — TNFα/NF-κB signalling,
interferon response, proliferation, IL6/JAK/STAT3 — now appears *directly* as a top-ranked upstream
regulator. This is independent corroboration from an orthogonal data model: gene-set overlap (§7)
and regulon-target coherence (§8) are different statistics on different inputs, and they name the
same circuit.

![Figure 12. TF-activity landscape. 732 transcription factors ranked by decoupleR ULM activity score on the CollecTRI regulon; red marks the STAT3 / NF-kB / interferon axis, which dominates the most-active extreme. STAT3 sits at the 97th percentile among all regulators.]({{artifact:art_28fbdc3e-f3e9-4acd-8452-39f40480d8eb}})

**The activity score is not a black box — it is grounded in target behaviour.** Figure 13 opens up
STAT3's score. Its full CollecTRI regulon has **435 targets** (≈396 activating, 39 repressing); of
these, **371 were present in our 24,528-gene tested set** and so entered the scoring. Among those
371 measured targets, the ~340 STAT3 *activates* shift coherently upward in lesional skin (mean t =
+3.84), which is exactly the fingerprint of an active transcription factor. Individually notable
targets moving up include SOCS3 (the classic STAT3 feedback gene), VEGFA, ICAM1, and CDKN1A.

![Figure 13. STAT3 target coherence. Each point is one STAT3 target gene, split by whether STAT3 activates (+) or represses (-) it; the y-axis is the observed limma t-statistic in lesional skin. The activated targets shift up as a group - the mechanistic basis of STAT3's high activity score.]({{artifact:art_e0a88d4a-8f71-444b-99db-b88e86fda3ed}})

**An honest limitation.** The signal is carried by the **activation** arm. CollecTRI annotates far
more activating than repressing edges for STAT3 (≈396 vs 39), and the repressed targets in our data
sit near zero (mean t = +1.47) rather than clearly negative. The activity call is robust — it rests
on hundreds of concordant activating targets — but we do not over-interpret the repression side, and
a regulon is only ever as good as its curation [26]. Two further caveats generalize to all TF-activity
work: activity is inferred, not measured (no protein or phospho-STAT3 data here), and CollecTRI is
not psoriasis- or keratinocyte-specific, so tissue-specific regulatory edges may be missing.

**Where §8 leaves the STAT3 hypothesis.** We now have **four independent lines of evidence**, each
from a different method and a different view of the data: (1) the single gene is up (§6.1); (2) its
transcriptional neighbourhood is up (§6.1); (3) its pathway is coordinately activated (§7); and (4)
it is inferred as an *active upstream regulator*, in a hub cluster with NF-κB and STAT1 (§8). The
driver-versus-output logic is now explicit: proliferation is the loudest signature, but the
regulators driving it are the inflammatory hubs — STAT3 among them. This is exactly the structure a
"STAT3 is central" hypothesis predicts, and it motivates the next rung — resolving *which STAT3
protein* (activating α vs truncated β) is being produced, from exon- and junction-level data.

The design is deliberately **modular for the meta-analysis**: TF activity is computed from a single
per-study t-statistic vector, so each study in the multi-cohort analysis produces its own STAT3
activity score against the identical, checkpointed CollecTRI regulon, and those scores pool across
studies exactly as the gene-level statistics do (§ forthcoming).

---

## 9. STAT3 isoform splicing: which STAT3 protein is lesional skin making?

Every layer so far — the gene (§6), the neighbourhood (§6), the pathway (§7), the inferred
regulator activity (§8) — treats *STAT3* as one entity. But the gene encodes **two proteins with
opposite functions**, produced by alternative splicing at the 3′ (C-terminal) end of the transcript.
This section asks the most mechanistically specific question in the whole analysis: *when lesional
skin turns STAT3 up, which STAT3 is it making?*

**The two isoforms.** STAT3 mRNA is spliced into two dominant forms that share the entire N-terminal
and DNA-binding region and diverge only at the C-terminus [31,32]:
- **STAT3α** (full length, 770 aa, ~92 kDa): the canonical activator. It retains the C-terminal
  **transactivation domain (TAD)** including the Ser727 regulatory site. This is the
  pro-inflammatory, pro-proliferative STAT3 that drives IL-17, SOCS3, BCL2, and MYC — the isoform
  the rest of this paper is implicitly about.
- **STAT3β** (truncated, 722 aa, ~83 kDa): an alternative splice removes the TAD, replacing the last
  55 residues of α with 7 unique ones. It was first described as a **dominant-negative** regulator —
  it still dimerizes and binds DNA but cannot transactivate, so it competes α off its targets
  [31,32] — though later work showed it also has specific transcriptional functions of its own [33].

The distinction matters because a rise in total STAT3 signal could mean *more activator* (α) or a
compensatory brake (β), and only isoform-resolved data can tell them apart.

**Reading the splice choice directly from junction reads.** recount3 provides **exon–exon junction
counts** — the direct quantitative readout of splicing — so we do not have to infer isoforms from
gene-level expression. We established the α-versus-β distinguishing junction rigorously rather than
by transcript-ID lookup: from Ensembl/GENCODE transcript models we confirmed the two protein
products (α = 770 aa ending `…CATSPM`; β = 722 aa diverging at **amino acid 716** and ending in the
unique heptapeptide **`FIDAVWK`** — the textbook β signature), and traced them to a single
**alternative-donor splice at a shared acceptor** (chr17:42,317,181, minus strand): α uses donor
42,316,902, β uses donor 42,316,852. Both junctions are well-supported in this study (α 37,144
reads; β 1,991 total), and **all 178 samples have ≥ 22 junction reads** at the locus, so the ratio
is real signal, not sparse-count noise. We summarise the splice choice per sample as
**PSI_β = β / (α + β)** — the percent-spliced-in of the β isoform.

**Result 1 — the β fraction rises significantly in lesional skin.** PSI_β increases from a mean of
**4.97 %** in normal skin to **6.14 %** in lesional skin (medians 5.0 % → 5.6 %; Wilcoxon
**p = 0.017**; Figure 14A). In absolute terms *both* isoforms go up — the activating α junction more
than doubles (mean 137 → 272 reads) and β rises in parallel (6.7 → 15.1; Figure 14B) — but β rises
*proportionally faster*, shifting the ratio toward the truncated isoform.

![Figure 14. STAT3 isoform splicing. (A) Percent-spliced-in of the truncated STAT3-beta isoform, PSI beta, per sample by group; beta inclusion rises in lesional (PP) versus normal (NN) skin (Wilcoxon p = 0.017). (B) Absolute alpha- and beta-defining junction read counts (log10); both isoforms are elevated in lesional skin with alpha dominant throughout.]({{artifact:art_9a97623e-72b4-41f8-9d06-8f29d5b4a044}})

**This is a genuinely two-sided result, and we read it honestly.** The absolute doubling of α is
fully consistent with the STAT3-activation story built in §§6–8: the activator isoform is
unambiguously elevated. But the *ratio* moving toward β is the opposite of a naive "everything tips
to the activator" prediction. Two mechanisms could produce it, and our data cannot yet separate
them: (i) a **feedback brake** — β is itself partly IL-6/STAT-inducible, so a rising β fraction may
be the tissue attempting to restrain runaway α signalling; or (ii) a **spliceosome perturbation** —
a globally dysregulated spliceosome would shift many splice ratios, this one incidentally among
them, without STAT3-specific regulation. The effect size is modest (a ≈ 1.2-percentage-point ratio
shift), so we frame it as *"the splice ratio is measurably perturbed toward β, in a context where
absolute activating-α is strongly up"* — not as evidence that β dominates.

**Result 2 — testing the spliceosome hypothesis via EFTUD2.** EFTUD2 (a core U5 snRNP spliceosome
GTPase, up in our DE table at log2FC +0.51) was our candidate splicing driver. If it drives the β
shift, PSI_β should rise with EFTUD2 expression. The data show a more interesting pattern (Figure
15): **across all 178 samples the correlation is null** (Spearman ρ = −0.02, p = 0.77), but **within
each group higher EFTUD2 tracks *lower* PSI_β** (lesional ρ = −0.27, p = 0.008; normal ρ = −0.21,
p = 0.053). This sign flip between the pooled and within-group analyses is a textbook **Simpson's
paradox**, driven by a confound the boxplot in Figure 15B makes explicit: EFTUD2 is itself strongly
up in lesional skin (p = 7 × 10⁻²⁵), the same direction PSI_β moves, so pooling the groups cancels
the real within-group relationship.

![Figure 15. PSI beta versus EFTUD2 expression. (A) Per-sample scatter with within-group linear fits: higher EFTUD2 is associated with LOWER STAT3-beta inclusion within both groups (Spearman rho -0.27 lesional, -0.21 normal), even though (B) EFTUD2 itself is markedly higher in lesional skin (Wilcoxon p = 7e-25). The opposite signs of the pooled (near zero) and within-group (negative) correlations are a Simpson's paradox.]({{artifact:art_9cebc8d9-7a11-4130-b6ab-d6aed14961a4}})

**What EFTUD2 tells us — and what it rules out.** The within-group direction is the key result:
*more* EFTUD2 goes with *less* β, i.e. EFTUD2 is associated with favouring the **activating α**
isoform, not the truncated one. This argues **against** the simple "EFTUD2 dysregulation makes the
β brake" hypothesis and is more consistent with EFTUD2 supporting efficient production of full-length
activating STAT3 — coherent with a pro-inflammatory keratinocyte state. It is a correlational result
in bulk tissue (confounded by shifting cell-type composition between lesional and normal skin, and
without isoform-resolved protein data), so it constrains rather than proves mechanism; but it is a
concrete, falsifiable constraint that the isoform-level data alone could deliver.

**Where §9 leaves the hypothesis.** We can now answer "which STAT3 is lesional skin making?" with
data rather than assumption: **predominantly the activating α isoform, in strongly increased absolute
amount, with a small but significant relative increase in the β brake, and a within-group EFTUD2
relationship that points toward α rather than β.** The splicing layer neither simply confirms nor
refutes the STAT3-centrality thesis — it sharpens it, and it defines the natural experimental
follow-up (isoform-resolved protein / phospho-STAT3α, and single-cell data to remove the cell-type
confound). The PSI computation is modular in exactly the way the TF-activity step is: one junction
ratio per sample, poolable across studies in the meta-analysis.

---

## 10. Sample clustering and co-expression structure (descriptive)

Every result so far has *tested* a difference we defined in advance (lesional vs normal). This
section does the complementary thing: it asks the data to organise *itself* and checks whether that
unsupervised structure recovers the disease groups we assigned — a quality-control and
sanity-check step, not a new hypothesis test.

**A note on why this is safe (the "double-dipping" caveat).** A well-known statistical trap is to
use the *same* data to both define groups (e.g. by clustering) and then test for differences between
them — the test is then guaranteed to find a difference because the groups were drawn around it
("circular analysis" / double-dipping). **We are not at risk here:** our NN/PP labels come from the
*a-priori* clinical sample annotation, not from clustering. The heatmaps below are therefore purely
**descriptive** — we colour samples by their pre-existing labels and observe whether unsupervised
clustering happens to agree. Agreement is reassurance about data quality; it is not, and is not
presented as, evidence for the group difference (that evidence is the limma test in §5).

**Sample-sample correlation (Figure 16).** We computed all pairwise sample correlations over the
3,477 DE genes, using both **Pearson** (linear) and **Spearman** (rank-based, robust to outliers and
non-linearity), and clustered samples by 1 − correlation (average linkage). Both metrics produce the
same clean **2 × 2 block structure**: lesional samples correlate tightly with each other, normal
samples with each other, and the two groups are visibly less similar across the block boundary. The
unsupervised dendrogram's top split corresponds to the clinical group almost perfectly — independent
confirmation of the PC1 result from §3 (PC1 ≈ 90 % of between-sample variance aligns with disease
status). That Pearson and Spearman agree tells us the separation is not an artefact of a few
high-leverage genes or a non-linear scaling effect; it is a broad, rank-stable signal.

![Figure 16. Sample-sample correlation heatmaps over the 3,477 DE genes, Pearson (left) and Spearman (right). Rows and columns are samples, coloured by clinical group (NN blue, PP red) in the side bars; colour encodes pairwise correlation. Both metrics resolve the same two tight blocks (lesional, normal); the unsupervised dendrogram's primary split matches the clinical labelling. A small number of samples sit nearer the opposite block (the PCA outliers noted in the roadmap).]({{artifact:art_23e6304d-d0c6-4eaf-8878-a6ef3f276648}})

The heatmaps also make the **outliers** from §3 concrete: a few samples sit closer to the opposite
group's block than to their own (visible as off-colour stripes in the annotation bars). These are the
candidates for the meta-analysis QC step — one normal sample sitting in the lesional cloud and one
lesional sample in the normal cloud — worth revisiting for possible mislabelling or atypical biology,
but too few to affect the group-level conclusions.

**DE-gene expression modules (Figure 17).** To see the *genes* rather than the samples, we drew the
top 120 DE genes (60 most up- and 60 most down-regulated) as a gene × sample heatmap of row
z-scored log2-CPM, clustering both axes by correlation distance. The picture is the canonical
psoriasis transcriptome:

- The **up-in-lesional module** is the IL-17/keratinocyte-inflammation program moving as a single
  coherent block: the antimicrobial peptides and S100 alarmins (`DEFB4A/B`, `S100A7/A7A/A8/A9/A12`,
  `PI3/elafin`, `LCN2`), the IL-17-axis cytokines and chemokines (`IL17A`, `IL19`, `IL20`, `IL36A/G`,
  `CXCL1`, `CXCL8`, `CCL20`, `NOS2`), the protease inhibitors (`SERPINB3/4/11`), and the cornified-
  envelope genes (`SPRR2` family, `LCE3A/C/D`). This is the *same program* the pathway analysis (§7)
  found up as a unit and the network analysis (§8) placed under STAT3/NF-κB control — here it is,
  visible directly at the level of individual genes and individual patients.
- The **down-in-lesional module** is dominated by normal-skin lipid metabolism and terminal-
  differentiation / xenobiotic genes (`AWAT1/2`, `DGAT2L6`, `MOGAT1`, `ELOVL3`, `FADS2`, `CYP1A1`,
  `CYP2W1`, `UGT2A1/2`), consistent with the loss of the normal epidermal-barrier lipid program in
  lesional skin.

![Figure 17. Top 120 DE genes (60 up, 60 down in lesional skin) as row z-scored log2-CPM, both axes clustered by correlation distance. Column side bar = clinical group; row side bar = DE direction. The up-in-lesional block is the IL-17/keratinocyte antimicrobial-inflammation program (DEFB4, S100A7-A12, IL17A, IL19/20/36, SERPINB, SPRR/LCE); the down block is normal-skin lipid-metabolism and differentiation genes.]({{artifact:art_5575bd40-f42c-4711-8a92-92e0a2e6b10b}})

**What we deliberately did *not* do yet: `clust`.** The user raised `clust` (Abu-Jamous & Kelly,
2018) — a method that extracts *tight, non-overlapping co-expressed gene modules* and, crucially,
can require that a module be co-expressed **consistently across multiple datasets simultaneously**.
That multi-dataset consistency is its headline strength, and it is precisely what the forthcoming
multi-study meta-analysis needs. In a single two-group study, "co-expression across conditions"
largely collapses into "up together vs down together," which the correlation heatmap above already
shows. We therefore **defer `clust` to the meta-analysis stage**, where several psoriasis cohorts
can be clustered jointly to find modules robust across studies — and where we can then locate STAT3's
module and GO-enrich each module. The descriptive heatmaps here establish the single-study structure
that clust will generalise.

---

## 11. Multi-study meta-analysis: from one cohort to a field-wide answer

Everything up to this point rested on a single study (SRP035988, 178 samples). That study is
large and clean, but one cohort cannot tell us whether a finding is a property of **psoriasis**
or a property of **that particular experiment** - its patients, its sequencing batch, its
biopsy protocol. The remedy is a **meta-analysis**: repeat the analysis in several independent
cohorts and combine the per-study results statistically. A signal that survives is one we can
attribute to the disease. A signal that does not survive is one we should stop trusting, no
matter how good its single-study p-value looked. This section is where the STAT3 hypothesis is
put to that test.

### 11.1 Choosing cohorts: comparable measurements only

recount3 hosts a large number of psoriasis-related RNA-seq projects, but a meta-analysis is only
valid when it pools **like with like**. We began from a curated list of 35 psoriasis studies and
filtered on three axes:

- **Species.** 12 mouse studies were removed - mouse imiquimod models are not human plaque
  psoriasis.
- **Tissue compartment.** Studies profiling sorted cells (keratinocyte monolayers, PBMCs,
  neutrophils, dendritic cells, mesenchymal stem cells), whole blood, or skin explants were
  removed. Our anchor is **whole-skin bulk** RNA-seq; mixing a blood transcriptome into that
  pool would inject a between-tissue difference far larger than the lesional-vs-healthy signal
  we are trying to measure.
- **Disease definition.** One large study (SRP165679) contained both atopic dermatitis and
  psoriasis patients; we retained only its psoriasis and healthy-control samples (recovering the
  diagnosis from the sample titles, since the machine-readable attribute field carried only skin
  type) and dropped its 54 atopic-dermatitis samples.

Five whole-skin human studies survived, together contributing **548 usable samples** - roughly
three times the single study (Figure 22). Their per-condition composition is uneven, which turns
out to matter for which contrasts each can support:

![Figure 22. Included studies and contrast coverage]({{artifact:art_85ca3348-1752-477d-9ac1-45d2b9f406db}})

*Figure 22. (A) The five retained studies, coloured by sample class (lesional PP, peri-lesional
PN, healthy NN). Total n = 548 (341 PP, 82 PN, 125 NN). SRP076982 is flagged as shallow. (B)
Contrast coverage: which study can contribute to each of the three comparisons.*

One quality flag is load-bearing. **SRP076982** (259 samples, the largest cohort and our main
source of peri-lesional tissue) is **shallow**: a median of ~4 million reads per sample, versus
27-41 million in the other four studies. We verified this against the STAR input-read count, so
it is a genuine property of the study, not a parsing error. Shallow sequencing does not
disqualify a study - `filterByExpr` simply drops the genes it cannot measure reliably, and the
inverse-variance pooling below automatically down-weights a noisier cohort through its larger
standard errors - but it means any result leaning heavily on SRP076982 must be read with care.

### 11.2 Three contrasts, by design

Rather than only ask "lesional vs healthy," we follow a three-contrast design that dissects the
disease gradient:

1. **Lesional vs healthy (PP vs NN)** - the main disease signature. Supported by 3 studies.
2. **Peri-lesional vs healthy (PN vs NN)** - does clinically normal-looking skin in a patient
   already carry a molecular signature? Supported by 2 studies.
3. **Peri-lesional vs lesional (PN vs PP)** - the within-patient gradient from uninvolved skin to
   plaque. Supported by 4 studies.

The scientific value of the peri-lesional arm is that it separates "what is different about a
psoriasis patient's skin everywhere" from "what is different about a plaque specifically." Novel
biology often hides in that distinction.

### 11.3 Per-study differential expression on a harmonised gene space

Because every recount3 study is quantified against the **same annotation** (Gencode v26), the
gene space is harmonised by construction - no cross-platform probe mapping is needed, which is
the single biggest source of noise in older microarray meta-analyses. For each study and each
contrast we re-ran the exact pipeline validated in sections 4-6: `filterByExpr` -> TMM ->
`voom` -> `limma` moderated t-tests, collapsing to unique gene symbols. The function reproduced
the single-study anchor result almost exactly (STAT3 logFC +1.156 vs the +1.16 obtained
earlier; 3,478 significant genes vs 3,477 - a one-gene difference from a filtering tie), which
confirms the harmonised pipeline is faithful.

The per-study volcanoes (Figure 23) make the power differences visible: the two well-powered
studies show dense significant signatures with the core genes (STAT3, SOCS3, S100A7, DEFB4A,
IL17A) far up-and-right, while the 4-vs-4 study has essentially no power - its cloud sits at the
origin. This is exactly the situation inverse-variance weighting is built for.

![Figure 23. Per-study lesional-vs-healthy volcanoes]({{artifact:art_5a5f54d2-69c6-4480-9ce2-448c4598fc13}})

*Figure 23. Per-study PP-vs-NN volcano plots. STAT3 and core psoriasis genes (black) are up in
every study that has power to see them. Red = FDR < 0.05 and |logFC| > 1.*

### 11.4 Pooling: random-effects meta-analysis, gene by gene

For every gene measured in at least two studies, we combined the per-study log-fold-changes with
an **inverse-variance random-effects** model (DerSimonian-Laird estimate of the between-study
variance tau-squared). We implemented the closed-form DL estimator directly and validated it
against `metafor::rma` - the two agree to four decimal places on STAT3 and every spot-checked
gene, so the fast vectorised version is exact.

Two terms carry the interpretation:

- **Random-effects vs fixed-effects.** A fixed-effect model assumes every study is estimating
  the *same* true effect and differences are pure sampling noise. A random-effects model allows
  the true effect to *vary* across studies (different patient populations, protocols) and is the
  honest default when cohorts are heterogeneous. It widens the pooled confidence interval to
  reflect that extra uncertainty.
- **I-squared and tau-squared** quantify heterogeneity: tau-squared is the estimated
  between-study variance in the true effect; I-squared is the percentage of total variation
  attributable to that (rather than to within-study sampling). High I-squared means the studies
  disagree about the *magnitude* of an effect - not necessarily its direction.

An honest result emerged immediately, and it is worth stating plainly because it corrects a
naive expectation. **Raw gene count is not where this meta-analysis wins.** Because the anchor
study alone is large and clean, it detects 18,771 genes at FDR < 0.05 - *more* than the pooled
12,036. Random-effects pooling deliberately widens standard errors to absorb heterogeneity,
which regularises borderline single-study calls downward. That is the correct, conservative
behaviour, not a loss of power.

What the meta-analysis genuinely buys is **replication**. A gene significant in the pooled model
is supported by independent cohorts, not one experiment. At a matched strict threshold (FDR <
0.05 and |pooled logFC| > 1), **264 genes are gained** over the single study - genes that were
not individually significant before but that replicate consistently across cohorts, including
immune genes such as CCR1, C1QB, PTPN7, RHOH, and MMP25 (Figure 24A, green).

![Figure 24. Random-effects meta-analysis, PP vs NN]({{artifact:art_12b69b7a-b00b-4359-a6f8-c1feba665163}})

*Figure 24. (A) Volcano of pooled random-effects effects; green marks the 264 cross-study-
replicated genes gained over the single study. (B) The distribution of I-squared among
meta-significant genes is bimodal: 524 genes have I-squared < 25% (near-identical effects across
cohorts - the reproducible core), while 1,279 have I-squared > 75% (direction agrees but
magnitude varies).*

### 11.5 The STAT3 verdict: three levels, all replicated

The forest plot is the canonical way to display a single gene's meta-analysis. For STAT3 across
the three PP-vs-NN studies (Figure 25), the per-study estimates (+1.16, +1.63, +0.42) pool to
**+1.25 [0.81, 1.69]**, with I-squared = 91%, tau-squared = 0.113, and a heterogeneity test
Q(2) = 22.4, p = 1.3e-5. STAT3 is unambiguously up in psoriatic skin across cohorts; the
heterogeneity is in *how much*, not *whether*, and the random-effects diamond honestly rests on
a wide base to reflect that. The pooled *magnitude* deserves one further caveat at k = 3, and the
*direction* one further robustness check; both are supplied in section 11.5a, which re-fits this
forest under the few-study-appropriate method and drops each cohort in turn.

![Figure 25. STAT3 random-effects forest]({{artifact:art_7fbc43c2-dc6d-49fb-a1c7-1c0a2a30d8a9}})

*Figure 25. STAT3 forest plot (PP vs NN). Per-study squares sized by precision; pooled
random-effects estimate as the red diamond. High I-squared reflects the larger effect in
SRP165679, not any disagreement about direction.*

Extending to a panel of key genes (Figure 26) shows the whole disease axis is robust, from the
STAT3-pathway regulators (EFTUD2 +0.5, STAT3 +1.25, STAT1 +1.86, SOCS3 +2.16) up to the
antimicrobial and keratinocyte effectors (IL17A +5.6, S100A7 +8.1, DEFB4A +10.6). The one
instructive non-result is **CXCL8**: a large point estimate (+5.5) but I-squared = 96% and a
confidence interval that crosses zero, so after heterogeneity correction it is not significant -
exactly the kind of over-optimistic single-study call that random-effects meta-analysis
appropriately reins in.

![Figure 26. Pooled effects for key psoriasis genes]({{artifact:art_f1e97683-d501-4339-8ed1-35f9b414b761}})

*Figure 26. Key psoriasis genes ordered by pooled logFC (PP vs NN). All significant except
CXCL8, whose extreme heterogeneity (I-squared 96%) widens its interval across zero.*

### 11.5a Robustness of the STAT3 gene-level estimate

The STAT3 forest pools three studies with considerable heterogeneity (I-squared = 91%). That is
exactly the regime - few studies, large I-squared - where the DerSimonian-Laird (DL) confidence
interval is known to be too narrow: DL builds its interval from a normal quantile while treating
tau-squared as if it were known exactly, but with only three studies tau-squared is itself
estimated with wide uncertainty, so the nominal 95% interval under-covers [44]. The
recommended remedy is the **Hartung-Knapp-Sidik-Jonkman (HKSJ)** modification, which replaces the
normal quantile with a t-distribution on k-1 degrees of freedom and re-weights, propagating the
tau-squared uncertainty into the interval [44, 45]. We therefore report HKSJ as the primary method
for this forest and retain DL as a sensitivity analysis, following the recommendation that analysts
using the modification be ready to justify it and show the conventional result alongside [45].

The point estimate is stable across methods; the interval is not (Figure 31A):

| Method | Estimate | 95% CI | p |
|---|---|---|---|
| DL + z (original) | +1.25 | [+0.81, +1.69] | 3e-08 |
| DL + Knapp-Hartung | +1.25 | [+0.02, +2.48] | 0.049 |
| **REML + HKSJ (primary)** | **+1.21** | **[-0.11, +2.53]** | **0.059** |

Under HKSJ the STAT3 gene-level interval widens roughly threefold and its lower bound crosses
zero (p = 0.059). The original p = 3e-08 was an artifact of the z-based interval at k = 3, not a
robust magnitude claim. The width is driven almost entirely by the smallest cohort (SRP126422,
n = 8, SE approximately 0.48), which carries negligible weight but, through its influence on
tau-squared, lets the t-correction open the interval; the two well-powered cohorts each place
STAT3 clearly up with tight intervals.

Crucially, this tempers the *magnitude* claim without touching the *direction* claim, which is
what the meta-analysis was built to test. Leave-one-out re-pooling (Figure 31B) keeps STAT3
positive in every case - +0.95, +1.39, and, with the large anchor cohort (SRP035988) removed,
**+1.12** - so the up-regulation is genuinely cross-cohort and not an artifact of the single
dominant study. The defensible statement is therefore stronger than a lone p-value: **STAT3's
up-regulation is reproducible across all cohorts (robust to leave-one-out, including
anchor removal); its pooled magnitude is approximately +1.2 but imprecisely estimated at k = 3,
which few-study-appropriate inference makes explicit.**

![Figure 31. STAT3 forest robustness]({{artifact:art_c793f88c-9806-46f9-95d2-0b570a3a8e50}})

*Figure 31. (A) Pooled STAT3 effect under three methods: the estimate is stable near +1.2 but the
interval triples in width and crosses zero once few-study uncertainty is honoured (HKSJ). (B)
Leave-one-out re-pooling (DerSimonian-Laird): every re-pool stays positive (+0.95 to +1.39), and
dropping the large anchor cohort keeps STAT3 up at +1.12 - the direction is not anchor-driven.*

STAT3 being an up-regulated gene is only the first level. The stronger claim - that STAT3 sits
at the **centre of the regulatory network** - is tested by transcription-factor activity. Using
decoupleR's univariate linear model on the CollecTRI regulon (435 STAT3 targets), we scored
STAT3 activity from each study's per-gene t-statistics. STAT3 ranks **19th of 732 TFs in the
anchor and 20th of 727 in the Tsoi cohort** - the top 3% in both well-powered studies - and
top 6% even in the tiny one. Pooled by a sample-size-weighted Stouffer combination, STAT3
activity gives **Z = 10.1, p = 5e-24** (Figure 27).

![Figure 27. STAT3 transcription-factor activity across studies]({{artifact:art_69a9aa93-291f-4497-a8c0-91b3ea747c77}})

*Figure 27. STAT3 regulon activity (decoupleR ULM z-score) per study and pooled. Per-study
ranks among ~730 TFs are annotated; STAT3 is a top-percentile hub in every cohort.*

At the pathway level (Figure 28), per-study CAMERA on the Hallmark collection, Stouffer-combined,
gives a coherent and reproducible landscape: interferon responses (alpha and gamma), proliferation
(MYC, E2F, G2M - the keratinocyte hyperproliferation programme), and the immune-signalling core.
**IL6-JAK-STAT3 ranks 9th of 50 Hallmark sets** (pooled Z = 8.5), up in every study.

![Figure 28. Pooled Hallmark pathway enrichment]({{artifact:art_6f045a3b-09b0-4ea5-8ec2-eba4da141212}})

*Figure 28. Top pooled Hallmark pathways (PP vs NN, CAMERA + Stouffer over 3 studies). IL6-JAK-
STAT3 (red) is firmly in the reproducible enriched set. Z capped at 38 for display.*

Together Figures 25-28 make the STAT3 case at three independent levels - the gene is up, its
regulon is enriched as a pathway, and its inferred TF activity is a top-ranked network hub -
**all replicated across cohorts.** This is the central positive result of the meta-analysis.

### 11.6 The STAT3 isoform switch does NOT replicate

The single-study analysis (section 9) found that lesional skin shifts toward the truncated,
dominant-negative **STAT3-beta** isoform: PSI-beta (the beta fraction of STAT3 splicing) rose
from 4.97% in healthy to 6.14% in lesional skin (Wilcoxon p = 0.017). The meta-analysis lets us
ask whether that isoform switch is a real feature of psoriasis. We extracted the identical
junction signal (shared acceptor chr17:42,317,181; alpha donor 42,316,902; beta donor
42,316,852) from every study and pooled the PP-vs-NN shift, restricting to samples with junction
depth >= 20.

It does not replicate (Figure 29). The anchor reproduces exactly (+1.17 percentage points,
p = 0.017), but the well-powered, deeply sequenced SRP165679 shows **no shift at all**
(-0.07 pp, p = 0.70). Pooled, the effect is **+0.78 pp [-0.05, +1.60], I-squared = 9%,
p = 0.065** - not significant.

![Figure 29. STAT3 beta-inclusion does not replicate across studies]({{artifact:art_3ae96d0f-708e-432b-90d3-37084480a7d2}})

*Figure 29. (A) Forest of the PSI-beta shift (lesional minus healthy). The anchor's positive
effect is not reproduced in the deep replication cohort; the pooled estimate is not significant.
(B) Per-study PSI-beta distributions (depth >= 20); the anchor shows PP > NN, SRP165679 shows no
gap.*

This is the meta-analysis doing its job. The distinction matters for the STAT3 hypothesis and
must be stated carefully:

- **STAT3 the gene and STAT3 the pathway are confirmed and strengthened** by adding cohorts.
- **The specific isoform-switch mechanism** - that psoriatic skin shifts toward dominant-negative
  STAT3-beta - was a single-study observation that **does not generalise.** The anchor's
  p = 0.017 is most likely a study-specific effect (real in that cohort, or a depth/batch
  artifact), not a reproducible property of the disease.

A negative result of this kind is not a disappointment; it is the difference between a robust
target and an over-interpreted detail, and it is precisely why the meta-analysis was worth doing.

### 11.7 A cross-study co-expression module

Finally, `clust` (deferred from section 10) was run to find gene modules co-expressed
**consistently across studies**, not within any single one. Restricting to the 2,380
meta-significant genes and the three well-powered cohorts, clust returned one high-confidence
consensus module of **65 genes, all up in lesional skin** (median logFC +3.43; Figure 30). Its
composition fuses two biological programmes:

- an **IL-17-driven antimicrobial / keratinocyte programme** - S100A7/A7A/A8/A9, DEFB103B, PI3
  (elafin), IL36A/G, CCL20, SERPINB3/4/11/13, the SPRR2 and LCE3 cornified-envelope families,
  KYNU, TGM1; and
- a **mitotic / proliferation cassette** - BIRC5 (survivin), CCNB1, CEP55, PLK1, PTTG1, RRM2.

![Figure 30. Cross-study co-expression module (clust C0)]({{artifact:art_4cde56c8-75bb-4b28-a402-3dc03f0d5275}})

*Figure 30. Consensus module C0: 65 genes co-expressed across all three large studies. The
per-sample module score is lowest in healthy skin, highest in lesional, and - in the one study
with all three arms (SRP165679) - places peri-lesional cleanly in between (NN -0.73 < PN -0.33 <
PP +1.31).*

![Figure 30b. Module logFC per study]({{artifact:art_69731df9-8dd1-49db-9388-04dc28e4fca9}})

*Figure 30b. Per-study PP-vs-NN log-fold-change for each module gene, confirming the module rises
coordinately in each cohort.*

The most striking feature is the **gradient**: the module score orders NN < PN < PP, and in the
cohort that measured all three, peri-lesional skin is a genuine intermediate. This is direct
molecular support for the idea that clinically uninvolved skin in a psoriasis patient is already
a partially-activated state, not simply healthy skin - a hypothesis the three-contrast design was
built to test.

A methodological note: clust is deliberately conservative. Requiring co-expression across all
five studies (including the 6- and 12-sample cohorts) returned zero modules; the one module above
emerged only when we restricted to the three large studies and loosened the tightness parameter.
clust reports what genuinely replicates, and 65 genes moving together across three independent
cohorts is a robust module by that standard.

### 11.8 What the meta-analysis changed

Bringing four additional cohorts to bear did three things. It **confirmed and strengthened** the
central STAT3 result at gene, pathway, and network-hub levels, all now cross-cohort replicated
(Figures 25-28). It **retired** a single-study finding - the STAT3-beta isoform switch - that did
not generalise (Figure 29). And it **added** 264 replicated genes and a coherent 65-gene disease
module with a healthy -> peri-lesional -> lesional gradient (Figures 24, 30). The net effect is a
STAT3 hypothesis that is narrower but far more defensible: STAT3 is a reproducible,
early-activated integrator of the psoriasis inflammatory programme — up in every cohort and robust
to leave-one-out (its pooled magnitude imprecise at k=3, §11.5a; its regulon activity a top-3%
transcription-factor signal) — while the isoform-level mechanism remains, at most, a
cohort-specific observation requiring dedicated deep-junction data to settle.

---

## 12. A molecular staging axis: in what order does the programme assemble?

Part I established *what* is reproducibly deranged in lesional psoriasis and placed STAT3 at the
centre of that programme. It did so with a case-control framing — lesional (PP) versus healthy
(NN) — which answers "what is different in disease" but not "what comes first." Many cohorts,
however, also sample **peri-lesional** skin (PN): clinically uninvolved skin from the same
patients. That third group turns the binary PP-vs-NN contrast into an **ordered axis**,
NN → PN → PP, and lets us ask a different question of the same data: in what order does the
psoriasis programme assemble? This part reframes the analysis around that axis. STAT3 re-enters
here not as the discovery but as **validation** — and the axis adds one mechanistic detail the
pairwise contrast cannot see.

### 12.1 Why a staging axis (and why peri-lesional skin is the key)

A case-control transcriptome answers "what is different in disease." It cannot answer "what
comes first," because it samples only two states. Psoriasis offers a natural way past this:
uninvolved skin from a psoriasis patient is not healthy skin. It carries a subclinical
molecular signature, and it is sampled in several cohorts as the "peri-lesional" or
"non-lesional" group.

If uninvolved skin is genuinely intermediate - if most disease genes sit part-way between
healthy and lesional there - then the three groups form an **ordered axis**, and the *order*
in which genes and pathways move along it is informative about the sequence in which the
disease programme assembles. That is the object of study here.

This is a cross-sectional ordering, not a time-course in a single patient. "Early" and
"late" refer to position along the healthy->lesional axis, not to clock time. But because
the same molecular programme is being caught at three points of completion, the axis behaves
like a pseudo-temporal trajectory, and the ordering it reveals is testable and reproducible.

### 12.2 The staging axis is real: 85% of the psoriasis programme is monotonic

We used the deepest three-group cohort in our recount3 selection, SRP165679 (the Tsoi
atopic-dermatitis/psoriasis cohort [55]), keeping its psoriasis and healthy arms: NN=38
healthy, PN=27 peri-lesional, PP=28 lesional. Expression was the uniformly reprocessed,
normalised log-expression used throughout the meta-analysis of Part I.

For every gene we fit two complementary tests of trend across the ordered stage (coded
NN=0, PN=1, PP=2):

- a **linear-trend model** (limma, moderated t on the numeric stage covariate), and
- a **rank-based monotonicity check** (Spearman correlation of expression with stage).

Restricting to the 4,563 genes that are DE in the lesional-vs-healthy contrast (the
psoriasis programme itself):

- **85%** are **strictly monotonic** by group means (1,506 progressively up, 2,380
  progressively down);
- **99%** show a **clear monotonic rank trend** (|Spearman rho| > 0.3; median rho = 0.62);
- peri-lesional skin sits a **median 16%** of the way from healthy to lesional (IQR 5-23%).

![Figure 32. Molecular staging axis. Top 40 progressively-up and top 40 progressively-down monotonic genes (rows), samples ordered healthy (NN, blue) -> peri-lesional (PN, orange) -> lesional (PP, red). Up-genes transition blue->pale->red across the axis and down-genes do the reverse; the peri-lesional block is visibly intermediate.]({{artifact:art_1d7a9b1d-a5a8-4879-8d07-f4fb97b29efd}})

![Figure 33. Quantifying the axis. (A) 85% of lesional DE genes are strictly monotonic and 99% show a monotonic rank trend across NN<PN<PP. (B) Distribution of the peri-lesional position (fraction of the lesional change already reached at the PN stage): most genes cluster low, median 16%, showing peri-lesional skin is an early molecular stage.]({{artifact:art_5e16b6b5-0a6e-4b3b-adf4-93b6deeb0b98}})

The key methodological point: prior work frames this as three separate pairwise DEG lists
and observes qualitatively that PN "looks intermediate." Testing the ordering *directly*
turns that impression into a number - 85% monotonic, PN at 16% - with a per-gene position on
the axis that the next sections exploit.

### 12.3 A timing taxonomy: early, progressive, and late genes

Because each gene has a "fraction-of-lesional-reached-at-PN" coordinate, we can classify the
programme by *when* along the axis each gene moves:

- **Progressive (2,015 genes)** - smoothly rising through the peri-lesional stage (PN
  fraction 15-50%).
- **Late / PP-specific (1,721 genes)** - almost no change at PN, switching on only at the
  lesional step (PN fraction < 15%).
- **Early / PN-specific (155 genes)** - most of the change already present at peri-lesional
  (PN fraction >= 50%).
- **PN-divergent (672 genes)** - move opposite to the lesional direction at the PN stage.

Enriching each class against the Hallmark collection gives the classes a clear biological
identity.

![Figure 34. Gradient-gene timing taxonomy. (A) Class sizes. (B) Hallmark enrichment of the progressive (early/inflammatory) versus late (proliferation) classes: progressive genes are dominated by interferon-gamma/alpha, IL6-JAK-STAT3 and inflammatory/TNF-NF-kB signalling; late genes are dominated by G2M checkpoint and E2F targets.]({{artifact:art_a2802b70-d309-4d26-a15c-69bbe9c66c63}})

The progressive class - the genes already rising in peri-lesional skin - is an **immune /
interferon** programme (interferon-gamma FDR 7e-24, interferon-alpha 9e-23, IL6-JAK-STAT3
2.5e-5, inflammatory response and TNF-NF-kB). The late class - the genes that appear only in
established lesions - is a **proliferation** programme (G2M checkpoint FDR 1.8e-8, E2F
targets 2.4e-6). This is the central biological claim of the paper, and section 4 confirms
it at the pathway-score level.

### 12.4 Pathway timing: inflammation is early, proliferation is late

To confirm the taxonomy independently of any gene-classification threshold, we computed
per-sample single-sample GSEA (ssGSEA) scores for the key Hallmark programmes in all 93
samples and tracked each programme's mean score across the three stages.

![Figure 35. Pathway activation timing. Mean ssGSEA score (+/- SE) across NN -> PN -> PP for immune/interferon programmes (left) and proliferation programmes (right). Immune programmes lift clearly at the peri-lesional stage; proliferation programmes stay near baseline until the lesional stage.]({{artifact:art_9d639ec8-5300-4226-827f-eb1af9d153f7}})

Testing each programme two ways - an overall linear trend, and specifically whether it is
already elevated at peri-lesional (PN vs NN) - gives a clean split:

- **Already significantly up at peri-lesional:** interferon-alpha (PN-vs-NN FDR 0.009),
  interferon-gamma (0.02), MYC targets (0.009); these reach 18-23% of their lesional level
  by the PN stage.
- **Not yet up at peri-lesional:** G2M checkpoint (PN fraction 3%, PN-vs-NN FDR 0.58) and
  E2F targets (8%, FDR 0.19) - despite both having highly significant *overall* NN->PP
  trends (FDR ~1e-18). Their entire change happens at the lesional step.

The disease therefore has a molecular order that the standard PP-vs-NN contrast collapses:
the interferon/STAT inflammatory programme is the initiating event, detectable in
still-uninvolved skin, and keratinocyte hyperproliferation is a downstream consequence that
appears only once the lesion is established.

### 12.5 STAT3 and the IL-17 axis: validation of the staging framework

Part I's central result (§11) is that STAT3 is a reproducible network hub of the lesional
programme. Here we place the canonical hubs on the staging axis - as a check that
the framework recovers known biology, and to locate STAT3 in the disease sequence.

![Figure 36. Canonical psoriasis hubs on the staging axis. STAT3, STAT1, SOCS3, IL17A, IL23A, S100A7, DEFB4A and IL36G all rise monotonically NN -> PN -> PP, each sitting intermediate at the peri-lesional stage. Every trend is significant (FDR <= 6e-10).]({{artifact:art_4014f5c6-fbf3-4822-8c98-f743145922c7}})

All eight hubs rise monotonically, each intermediate at peri-lesional - the axis reproduces
the established psoriasis biology gene by gene. STAT3 thus enters as validation, not as the
discovery.

The axis also adds one detail a pairwise contrast cannot. Inferring STAT3 transcription-
factor activity (decoupleR ULM on the CollecTRI regulon, relative to the healthy baseline)
across the three stages:

![Figure 37. STAT3 transcription-factor activity on the staging axis. STAT3 activity rises monotonically (NN 0.01 -> PN 2.21 -> PP 7.43; trend p = 4e-16), reaching ~30% of the lesional level already at peri-lesional.]({{artifact:art_25cf5b00-633a-437e-9561-72ade6e1b720}})

STAT3 *activity* reaches ~30% of its lesional level at peri-lesional, while STAT3's own
*transcript* reaches only ~10% (section 2 coordinate). The gap is mechanistically coherent:
STAT3 acts through phosphorylation and target induction, so its regulon can light up in
uninvolved skin before the gene itself is strongly transcribed. On the staging axis, STAT3
activation is an early event.

### 12.6 What this adds to the existing non-lesional-skin literature

Two ideas here already exist, and we claim neither as novel. Peri-lesional skin as a
molecular intermediate is well established - Gudjonsson et al. described decreased lipid
biosynthesis and increased innate immunity in uninvolved skin [46]; a proteomic study is
explicitly titled "Intermediate Stage of Non-Lesional Psoriatic Skin" [47]; reviews describe
uninvolved skin as "somewhere in between" [48]; and single-cell and spatial studies reproduce
the intermediate profile [49, 50]. Cross-study integration of psoriasis transcriptomes also
exists [51, 52, 53].

Against that backdrop, the specific contributions of this analysis are:

1. **An explicit ordering statistic on a three-point axis**, replacing three pairwise DEG
   lists and a qualitative "PN sits between" with a direct trend + monotonicity test (85%
   monotonic, PN at a median 16%), and a per-gene axis coordinate.
2. **A gradient-timing taxonomy** that orders the disease programmes - inflammation early,
   proliferation late - a sharper statement than "both are dysregulated in lesions."
3. **STAT3 positioned on the axis** with the activity-before-transcription signature.
4. **A reproducible, recount3-anchored pipeline** (pipeline STEP 12; see Methods).

**Limitations, stated plainly.** The clean three-point axis is well-powered in one cohort
(SRP165679); other three-group cohorts in our selection are shallow or tiny, so the axis is
*anchored* in one deep study and *supported* by pooled trends rather than independently
replicated at full power. Peri-lesional effect sizes are small (median ~16% of the lesional
change), so the story is strongest at the pathway/trend level. And the axis is a
cross-sectional ordering of three sampled states, not a temporal trajectory in individual
patients.

### 12.7 Staging-axis conclusion

Adding the peri-lesional group turns the psoriasis case-control transcriptome into an ordered
molecular staging axis. Most of the psoriasis programme moves monotonically along it, with
peri-lesional skin an early intermediate; the interferon/STAT inflammatory programme is the
early event and keratinocyte proliferation the late one; and the canonical STAT3/IL-17 hubs
validate the axis while STAT3's activity marks it early. The framework reframes a well-known
lesional signature as a *sequence*, and it does so from public, uniformly reprocessed data
with a reproducible pipeline.

---

## 13. Synthesis: one inflammatory core, replicated and early

The two analyses in this paper view the same five cohorts through different lenses and arrive at
a mutually reinforcing picture.

The **meta-analysis** (Part I) asked whether the lesional psoriasis signature — and STAT3's
central place in it — reproduces across independent cohorts. It does: the interferon /
IL6-JAK-STAT3 / IL-17 inflammatory programme is enriched in every well-powered cohort, STAT3 is
up-regulated in all of them, and STAT3's regulon activity ranks in the top 3% of transcription
factors. The one honest correction is on *magnitude*, not *direction*: at three studies with high
between-cohort heterogeneity, the pooled STAT3 fold-change is imprecise (HKSJ 95% CI crosses
zero) even though its up-regulation is robust to leave-one-out including anchor removal (§11.5a).
The isoform-switch mechanism did not survive replication and was retired.

The **staging axis** (Part II) asked, of the same programme, *what comes first*. Peri-lesional
skin is a genuine early molecular stage (median 16% of the lesional change), most of the
programme moves monotonically along the axis, and the ordering is sharp: the inflammatory /
interferon / STAT programme is the **early** event — already detectable in still-uninvolved skin
— while keratinocyte hyperproliferation is the **late** event, appearing only once the lesion is
established.

Read together, the results converge on a single statement. Psoriasis has a reproducible,
STAT3-associated inflammatory core that (i) replicates across cohorts at the gene, pathway, and
network-activity levels, and (ii) initiates the disease sequence, with proliferation as a
downstream, lesion-stage consequence. STAT3's activity-before-transcription signature on the
axis is consistent with its role as an early-activated integrator of that core. Two cautions
carry through from the analyses themselves: STAT3's status here is as a reproducible, early
*downstream integrator* of upstream IL-23/IL-17 and IL-6 signalling rather than a demonstrated
causal driver, and the quantitative staging axis is anchored in one deep cohort and supported by
pooled trends rather than independently replicated at full power. Both are stated plainly so the
central claim — a replicated, early, STAT3-associated inflammatory core — rests only on what the
data support.

---

## Methods & reproducibility (running note)

All analysis is in **R** using **Bioconductor** packages. Core stack: `recount3` (data
access), `edgeR` (filtering, TMM normalization, CPM), with `limma`/`DESeq2` available for
differential expression, and `ggplot2` for figures.

**Pipeline provenance.** The steps documented here correspond to a numbered R script series
(`Step-1` … `step-12`). Key methodological upgrades made during this walkthrough, relative
to the original scripts, are recorded explicitly:

1. **Gene filtering** changed from a hard-coded "≥ 1 CPM in ≥ 60 samples" to design-aware,
   group-independent **`filterByExpr()`** (§4), for the statistical-validity reasons in [1].
2. **Group assignment** is verified by a group-count table with an explicit `NA` check (§3)
   before any downstream step.

**Exact operations run (SRP035988, gene level):**
- `create_rse(proj_info)` then `assay(rse,"counts") <- transform_counts(rse)` — coverage → counts.
- Group labels from `sra.sample_attributes` via string match (`normal_skin`→NN,
  `Psoriasis_skin`→PP); verified 83 NN / 95 PP / 0 NA.
- Collapse to gene symbols by summing counts of duplicated symbols (`rowsum`).
- `DGEList` → `filterByExpr(group=)` → subset → `calcNormFactors(method="TMM")`.
- `cpm(..., log=TRUE)` for visualization; `prcomp(t(log2cpm), scale.=FALSE)` for PCA.
- PCA–group association: Pearson correlation of PC scores with 0/1 group + Welch t-test +
  linear-model R² for PC1.
- **Pathway analysis:** limma t-statistic ranking of all 24,528 genes → `fgsea::fgsea()`
  (minSize 10, maxSize 500, multilevel/`eps=0`) against MSigDB Hallmark, Reactome, KEGG, and
  GO:BP collections retrieved via `msigdbr` (2026.1.Hs); ORA cross-check by `fisher.test()` on
  2×2 tables over the 3,477-gene DE hit list with the 24,528-gene tested universe as background.
  Ranking vector checkpointed as `gsea_ranks.rds`; all four result tables as
  `gsea_results_all.rds` / `.csv`.
- **Network analysis (TF activity):** CollecTRI signed regulon fetched from the OmniPath HTTP API
  (`omnipathdb.org`, human, 62,411 edges over 1,201 TFs); `decoupleR::run_ulm()` on the limma
  t-statistic matrix (`minsize=5`), BH-adjusted across TFs. Regulon checkpointed as
  `collectri_regulon.rds`, scored TFs as `tf_activity_collectri.rds` / `.csv`.
- **Isoform splicing (STAT3 α/β):** junction-level RSE via `create_rse(proj_info, type="jxn")`
  (1,816,196 junctions × 178). α- vs β-defining junctions identified from Ensembl/GENCODE transcript
  models (α donor chr17:42,316,902; β donor 42,316,852; shared acceptor 42,317,181; minus strand),
  confirmed by translation length (770 vs 722 aa) and the β `FIDAVWK` C-terminus. Per-sample
  PSI_β = β/(α+β); Wilcoxon NN vs PP. EFTUD2 correlation on aligned log2-CPM (Spearman, pooled and
  within-group). Checkpoints `stat3_isoform_psi.rds`, `stat3_isoform_eftud2.rds`.
- **Sample clustering / heatmaps (descriptive):** pairwise sample correlations (Pearson and
  Spearman) over the 3,477 DE genes, hierarchical clustering on 1 − correlation (average linkage),
  `pheatmap`. DE-gene module heatmap on the top 120 DE genes (row z-scored log2-CPM), both axes
  clustered by correlation distance. Groups are a-priori (no clustering-then-testing), so these are
  descriptive/QC figures, not tests. Figures `fig16_sample_corr_heatmap.png`,
  `fig17_degene_heatmap.png`.

**Data/versioning.** recount3 project **SRP035988**; annotation **G026** (Gencode v26 gene
sums, as returned by recount3). The filtered+normalized object is checkpointed as
`dge_filt_norm.rds`, and the full DE table (with BH and Storey columns) as
`de_results_full.rds`, for exact reproduction of every downstream step.

**One-command reproduction.** The entire pipeline — raw recount3 download → filtering + TMM →
limma-voom DE + Storey q-values → volcano and diagnostic figures — is captured in a single
master script, **`psoriasis_pipeline.R`**, with setup instructions and a verification checklist
in **`README.md`**. A reviewer can regenerate every number and figure in this section with
`Rscript psoriasis_pipeline.R`; each step writes a checkpoint so the run can resume after a
disconnect rather than restarting from the download.

---


**Staging axis (§12).** Cohort SRP165679 (recount3, Gencode v26, GRCh38), psoriasis + healthy
arms: NN=38, PN=27, PP=28; normalised log-expression as above. Per-gene trend by limma
moderated-t on numeric stage (NN=0, PN=1, PP=2) plus Spearman of expression versus stage;
monotonicity by strict ordering of group means; PN fraction = (mu_PN − mu_NN)/(mu_PP − mu_NN).
Timing taxonomy by PN fraction: early (≥0.50), progressive (0.15–0.50), late (<0.15),
PN-divergent (<0); Hallmark enrichment by Fisher's exact test against the trend-tested universe.
Pathway timing by per-sample ssGSEA (Barbie 2009) for key Hallmark sets, with a linear trend and
a PN-vs-NN Wilcoxon test per pathway. STAT3 activity by decoupleR ULM on the CollecTRI regulon
per sample relative to the healthy mean. Reproducible as pipeline STEP 12; figures and tables in
staging_axis/.

---

## References

[1] Bourgon R, Gentleman R, Huber W. **Independent filtering increases detection power for
high-throughput experiments.** *Proc Natl Acad Sci USA.* 2010;107(21):9546–9551.
doi:10.1073/pnas.0914005107.

[2] Smyth GK, edgeR maintainer guidance, Bioconductor Support forum — on the statistical
invalidity of group-oriented (group-aware) gene filtering. (Consistent with the
independent-filtering principle established in [1].)

[3] Chen Y, Lun ATL, Smyth GK. **From reads to genes to pathways: differential expression
analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.**
*F1000Research.* 2016;5:1438. doi:10.12688/f1000research.8987.2. *(Method reference for
`filterByExpr` and the edgeR quasi-likelihood workflow.)*

[4] Robinson MD, Oshlack A. **A scaling normalization method for differential expression
analysis of RNA-seq data.** *Genome Biology.* 2010;11:R25. doi:10.1186/gb-2010-11-3-r25.
*(The TMM method.)*

[5] Wagner GP, Kin K, Lynch VJ. **Measurement of mRNA abundance using RNA-seq data: RPKM
measure is inconsistent among samples.** *Theory in Biosciences.* 2012;131(4):281–285.
doi:10.1007/s12064-012-0162-3. *(RPKM inconsistent across samples; advocates TPM.)*

[6] Bullard JH, Purdom E, Dudoit S, Speed TP. **Evaluation of statistical methods for
normalization and differential expression in mRNA-Seq experiments.** *BMC Bioinformatics.*
2010;11:94. doi:10.1186/1471-2105-11-94. *(Total-count/RPKM scaling is sensitive to a few
high-count genes.)*

[7] Anders S, Huber W. **Differential expression analysis for sequence count data.** *Genome
Biology.* 2010;11:R106. doi:10.1186/gb-2010-11-10-r106. *(DESeq median-of-ratios size
factors.)*

[8] Love MI, Huber W, Anders S. **Moderated estimation of fold change and dispersion for
RNA-seq data with DESeq2.** *Genome Biology.* 2014;15:550. doi:10.1186/s13059-014-0550-8.
*(DESeq2.)*

[9] Li B, Dewey CN. **RSEM: accurate transcript quantification from RNA-Seq data with or
without a reference genome.** *BMC Bioinformatics.* 2011;12:323. doi:10.1186/1471-2105-12-323.
*(Origin of the TPM unit.)*

[10] Dillies M-A, Rau A, Aubert J, et al. **A comprehensive evaluation of normalization
methods for Illumina high-throughput RNA sequencing data analysis.** *Briefings in
Bioinformatics.* 2013;14(6):671–683. doi:10.1093/bib/bbs046. *(Benchmark recommending
TMM/DESeq-type methods over total-count/RPKM for DE.)*

[11] Robinson MD, McCarthy DJ, Smyth GK. **edgeR: a Bioconductor package for differential
expression analysis of digital gene expression data.** *Bioinformatics.* 2010;26(1):139–140.
doi:10.1093/bioinformatics/btp616. *(edgeR; NB model.)*

[12] McCarthy DJ, Chen Y, Smyth GK. **Differential expression analysis of multifactor RNA-Seq
experiments with respect to biological variation.** *Nucleic Acids Research.*
2012;40(10):4288–4297. doi:10.1093/nar/gks042. *(NB dispersion estimation, BCV.)*

[13] Law CW, Chen Y, Shi W, Smyth GK. **voom: precision weights unlock linear model analysis
tools for RNA-seq read counts.** *Genome Biology.* 2014;15:R29. doi:10.1186/gb-2014-15-2-r29.
*(voom precision weights.)*

[14] Ritchie ME, Phipson B, Wu D, et al. **limma powers differential expression analyses for
RNA-sequencing and microarray studies.** *Nucleic Acids Research.* 2015;43(7):e47.
doi:10.1093/nar/gkv007. *(limma; moderated t-test.)*

[15] Smyth GK. **Linear models and empirical Bayes methods for assessing differential
expression in microarray experiments.** *Stat Appl Genet Mol Biol.* 2004;3:Article 3.
doi:10.2202/1544-6115.1027. *(Empirical-Bayes moderated t-statistic.)*

[16] Anders S, Huber W. **Differential expression analysis for sequence count data.** *Genome
Biology.* 2010;11:R106. doi:10.1186/gb-2010-11-10-r106. *(DESeq NB model.)*

[17] Soneson C, Delorenzi M. **A comparison of methods for differential expression analysis of
RNA-seq data.** *BMC Bioinformatics.* 2013;14:91. doi:10.1186/1471-2105-14-91. *(DE method
benchmark.)*

[18] Benjamini Y, Hochberg Y. **Controlling the false discovery rate: a practical and powerful
approach to multiple testing.** *J R Stat Soc B.* 1995;57(1):289–300. *(FDR / BH procedure.)*

[19] Storey JD, Tibshirani R. **Statistical significance for genomewide studies.** *Proc Natl
Acad Sci USA.* 2003;100(16):9440–9445. doi:10.1073/pnas.1530509100. *(q-value, pi0.)*

[20] Subramanian A, Tamayo P, Mootha VK, et al. **Gene set enrichment analysis: a
knowledge-based approach for interpreting genome-wide expression profiles.** *Proc Natl Acad Sci
USA.* 2005;102(43):15545–15550. doi:10.1073/pnas.0506580102. *(original GSEA, running ES,
sample-label permutation.)*

[21] Korotkevich G, Sukhov V, Budin N, et al. **Fast gene set enrichment analysis.** *bioRxiv.*
2021. doi:10.1101/060012. *(fgsea multilevel/Monte-Carlo p-value estimation.)*

[22] Liberzon A, Birger C, Thorvaldsdóttir H, et al. **The Molecular Signatures Database (MSigDB)
hallmark gene set collection.** *Cell Syst.* 2015;1(6):417–425. doi:10.1016/j.cels.2015.12.004.

[23] Milacic M, Beavers D, Conley P, et al. **The Reactome Pathway Knowledgebase 2024.** *Nucleic
Acids Res.* 2024;52(D1):D672–D678. doi:10.1093/nar/gkad1025.

[24] Kanehisa M, Goto S. **KEGG: Kyoto Encyclopedia of Genes and Genomes.** *Nucleic Acids Res.*
2000;28(1):27–30. doi:10.1093/nar/28.1.27.

[25] Ashburner M, Ball CA, Blake JA, et al. **Gene Ontology: tool for the unification of
biology.** *Nat Genet.* 2000;25(1):25–29. doi:10.1038/75556. *(GO Biological Process.)*

[26] Müller-Dott S, Tsirvouli E, Vazquez M, et al. **Expanding the coverage of regulons from
high-confidence prior knowledge for accurate estimation of transcription factor activities.**
*Nucleic Acids Res.* 2023;51(20):10934–10949. doi:10.1093/nar/gkad841. *(CollecTRI signed regulon.)*

[27] Badia-i-Mompel P, Vélez Santiago J, Braunger J, et al. **decoupleR: ensemble of computational
methods to infer biological activities from omics data.** *Bioinformatics Advances.*
2022;2(1):vbac016. doi:10.1093/bioadv/vbac016.

[28] Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. **Benchmark and
integration of resources for the estimation of human transcription factor activities.** *Genome
Res.* 2019;29(8):1363–1375. doi:10.1101/gr.240663.118. *(DoRothEA; benchmarking of TF-activity
methods.)*

[29] Alvarez MJ, Shen Y, Giorgi FM, et al. **Functional characterization of somatic mutations in
cancer using network-based inference of protein activity (VIPER).** *Nat Genet.*
2016;48(8):838–847. doi:10.1038/ng.3593.

[30] Szklarczyk D, Kirsch R, Koutrouli M, et al. **The STRING database in 2023.** *Nucleic Acids
Res.* 2023;51(D1):D638–D646. doi:10.1093/nar/gkac1000. *(Protein–protein interaction networks.)*

[31] Schaefer TS, Sanders LK, Nathans D. **Cooperative transcriptional activity of Jun and Stat3β,
a short form of Stat3.** *Proc Natl Acad Sci USA.* 1995;92(20):9097–9101.
doi:10.1073/pnas.92.20.9097. *(First description of the STAT3β splice variant.)*

[32] Caldenhoven E, van Dijk TB, Solari R, et al. **STAT3β, a splice variant of transcription factor
STAT3, is a dominant negative regulator of transcription.** *J Biol Chem.* 1996;271(22):13221–13227.
doi:10.1074/jbc.271.22.13221. *(Definition of α vs β: C-terminal transactivation-domain truncation;
dominant-negative activity.)*

[33] Maritano D, Sugrue ML, Tininini S, et al. **The STAT3 isoforms α and β have unique and specific
functions.** *Nat Immunol.* 2004;5(4):401–409. doi:10.1038/ni1052. *(β is not merely dominant-
negative — it has specific transcriptional functions.)*

[34] Abu-Jamous B, Kelly S. **Clust: automatic extraction of optimal co-expressed gene clusters from
gene expression data.** *Genome Biol.* 2018;19(1):172. doi:10.1186/s13059-018-1536-8. *(Co-expression
module extraction, including consistent clustering across multiple datasets.)*

[35] Goeman JJ, Bühlmann P. **Analyzing gene expression data in terms of gene sets: methodological
issues.** *Bioinformatics.* 2007;23(8):980–987. doi:10.1093/bioinformatics/btm051. *(competitive vs
self-contained null hypotheses; the framework behind §7.1.)*

[36] Wu D, Smyth GK. **Camera: a competitive gene set test accounting for inter-gene correlation.**
*Nucleic Acids Res.* 2012;40(17):e133. doi:10.1093/nar/gks461. *(variance-inflation correction for
co-regulated genes.)*

[37] Wu D, Lim E, Vaillant F, et al. **ROAST: rotation gene set tests for complex microarray
experiments.** *Bioinformatics.* 2010;26(17):2176–2182. doi:10.1093/bioinformatics/btq401.
*(self-contained rotation test; mroast.)*

[38] Hänzelmann S, Castelo R, Guinney J. **GSVA: gene set variation analysis for microarray and
RNA-seq data.** *BMC Bioinformatics.* 2013;14:7. doi:10.1186/1471-2105-14-7. *(per-sample gene-set
enrichment scores.)*

[39] Barbie DA, Tamayo P, Boehm JS, et al. **Systematic RNA interference reveals that oncogenic
KRAS-driven cancers require TBK1.** *Nature.* 2009;462(7269):108–112. doi:10.1038/nature08460.
*(single-sample GSEA, ssGSEA — the per-sample method we implemented in base R.)*

[40] Conway JR, Lex A, Gehlenborg N. **UpSetR: an R package for the visualization of intersecting
sets and their properties.** *Bioinformatics.* 2017;33(18):2938–2940.
doi:10.1093/bioinformatics/btx364. *(UpSet plots for leading-edge overlap.)*

[41] DerSimonian R, Laird N. **Meta-analysis in clinical trials.** *Control Clin Trials.*
1986;7(3):177–188. doi:10.1016/0197-2456(86)90046-2. *(Random-effects model; moment estimator of
between-study variance tau-squared - the pooling used in §11.4-11.6.)*

[42] Viechtbauer W. **Conducting meta-analyses in R with the metafor package.** *J Stat Softw.*
2010;36(3):1–48. doi:10.18637/jss.v036.i03. *(metafor::rma; forest plots; I-squared/tau-squared -
validation reference for the vectorised DL implementation.)*

[43] Whitlock MC. **Combining probability from independent tests: the weighted Z-method is superior
to Fisher's approach.** *J Evol Biol.* 2005;18(5):1368–1373. doi:10.1111/j.1420-9101.2005.00917.x.
*(Stouffer weighted-Z combination used for pooled TF-activity and pathway enrichment, §11.5.)*

[44] IntHout J, Ioannidis JPA, Borm GF. **The Hartung-Knapp-Sidik-Jonkman method for random effects
meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird
method.** *BMC Med Res Methodol.* 2014;14:25. doi:10.1186/1471-2288-14-25. *(HKSJ under-coverage of
DL at few studies / high heterogeneity; primary method for the STAT3 forest, §11.5a.)*

[45] Jackson D, Law M, Rücker G, Schwarzer G. **The Hartung-Knapp modification for random-effects
meta-analysis: A useful refinement but are there any residual concerns?** *Stat Med.*
2017;36(25):3923-3934. doi:10.1002/sim.7411. *(Recommendation to report HKSJ as primary with the
conventional method as a sensitivity analysis, §11.5a.)*

---

*Document status: **sections 1–11 drafted.** Single-study rungs: differential-expression **results**
(3,477 significant genes; STAT3 and its IL-17/SOCS3 neighbourhood up), **pathway analysis** (GSEA
across Hallmark, Reactome, KEGG, GO:BP; the IL6/JAK/STAT3–IL-17–NF-κB–keratinocyte program up as a
unit; GSEA-vs-ORA contrast), **network analysis** (decoupleR/CollecTRI TF-activity: STAT3 an active
regulator at the 97th percentile, in a hub cluster with NF-κB/STAT1), **isoform splicing**
(junction-level α/β: PSI_β up in lesional skin, p = 0.017; EFTUD2 within-group correlation points
toward α), **sample clustering** (Pearson/Spearman correlation heatmaps recover the clinical groups;
DE-gene module heatmap shows the IL-17/keratinocyte block), and **enrichment robustness** (§7.1:
CAMERA/ROAST confirm the STAT3 regulon under correlation-corrected and self-contained tests;
per-sample ssGSEA scores track STAT3 and PSI_β; UpSet shows the immune pathways share one redundant
core). **Multi-study meta-analysis (§11):** 5 curated whole-skin human cohorts (548 samples);
per-study harmonised DE; gene-level random-effects pooling (DerSimonian-Laird, validated vs metafor);
STAT3 confirmed at gene (up in every cohort, pooled approximately +1.2; direction robust to
leave-one-out including anchor removal, magnitude imprecise at k = 3 - HKSJ 95% CI [-0.11, 2.53],
§11.5a), pathway (IL6-JAK-STAT3 rank 9/50), and network-hub
(pooled TF-activity Z = 10.1, top 3% of ~730 TFs) levels, all cross-cohort replicated; the STAT3-beta
isoform switch does NOT replicate (pooled p = 0.065); 264 replicated genes gained; one 65-gene
cross-study co-expression module (clust) with a healthy -> peri-lesional -> lesional gradient.
Figures 1–30 complete. Companion deep-dives: `normalization_deep_dive.md`,
`differential_expression_deep_dive.md`, `pathway_analysis_deep_dive.md`, `network_analysis_deep_dive.md`,
`spliceosome_analysis_deep_dive.md`, `meta_analysis_deep_dive.md`. Reproducible via
`psoriasis_pipeline.R` (STEP 11) + `README.md`.*

[46] Gudjonsson JE, Ding J, Li X, et al. Global gene expression analysis reveals evidence
for decreased lipid biosynthesis and increased innate immunity in uninvolved psoriatic
skin. J Invest Dermatol. 2009;129(12):2795-2804. doi:10.1038/jid.2009.173.

[47] Bata-Csorgo Z et al. Comprehensive Proteomic Analysis Reveals Intermediate Stage of
Non-Lesional Psoriatic Skin and Points out the Importance of Proteins Outside this Trend.
(PMC6684579).

[48] Study of Molecular Mechanisms Involved in the Pathogenesis of Immune-Mediated
Inflammatory Diseases, using Psoriasis As a Model. (PMC3347524).

[49] Profiling Long Noncoding RNA in Psoriatic Skin Using Single-Cell RNA Sequencing.
J Invest Dermatol. 2024. (builds on Tsoi et al. 2019, the SRP165679 cohort).

[50] Spatial transcriptomics stratifies psoriatic disease severity by emergent cellular
ecosystems. (PMC10502701).

[51] Cross-Study Homogeneity of Psoriasis Gene Expression in Skin across a Large Expression
Range. PLoS One. 2013. (PMC3537625).

[52] The integration of large-scale public data and network analysis uncovers molecular
characteristics of psoriasis. Hum Genomics. 2022;16:doi:10.1186/s40246-022-00431-x.

[53] Suarez-Farinas M, Lowes MA, Zaba LC, Krueger JG. Evaluation of the psoriasis
transcriptome across different studies by gene set enrichment analysis (GSEA). PLoS One.
2010;5(4):e10247.

[54] Li B, Tsoi LC, Swindell WR, et al. Transcriptome analysis of psoriasis in a large case-
control sample: RNA-seq provides insights into disease mechanisms. J Invest Dermatol.
2014;134(7):1828-1838. (SRP035988, the anchor cohort).

[55] Tsoi LC, et al. Atopic dermatitis is an IL-13-dominant disease with greater molecular
heterogeneity, contrasted with the IL-17 axis of psoriasis. J Invest Dermatol. 2019.
(SRP165679, the three-group cohort anchoring the staging axis).
