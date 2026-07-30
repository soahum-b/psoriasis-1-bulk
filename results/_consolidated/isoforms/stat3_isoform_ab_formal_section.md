# Formal evaluation of the STAT3 α/β isoform ratio in psoriasis

*Consolidated, current version. Supersedes the 3-study framing in "Case study in honesty"
(meta_analysis_deep_dive.md §8), which reported the original +0.78 pp / p=0.065 pool before the
cohort set was extended.*

## Scope: what bulk RNA-seq can and cannot say about the six isoforms

The six catalogued STAT3 protein isoforms arise by three different mechanisms, and only one leaves
a signature that bulk RNA-seq can quantify:

| Isoform(s) | ~MW | Origin | Measurable from Recount-3? |
|---|---|---|---|
| STAT3α | 92 kDa | full-length canonical transcript | yes — reference transcript |
| STAT3β | 83 kDa | **alternative 3′ splice-site** in the terminal coding exon | **yes — distinct mRNA junction** |
| STAT3γ | 72 kDa | proteolytic cleavage (neutrophils) | no — mRNA identical to α |
| STAT3δ | 64 kDa | proteolytic cleavage (granulocyte maturation) | no — mRNA identical to α |
| STAT3ε, STAT3ζ | — | poorly defined; no agreed transcript model | no — no reliable RNA signature |

Only the **α/β split** is a splicing event with resolvable junctions, so it is the only isoform
question this data type can address directly. γ and δ are post-translational cleavage products with
mRNA identical to α; ε/ζ have no consensus transcript model. This section therefore evaluates the
α/β ratio and states the ceiling explicitly.

## The junction signal

We quantify the β fraction as **PSI-β = β-junction / (α-junction + β-junction)**, using the shared
acceptor chr17:42,317,181 with the α donor at 42,316,902 and the β donor at 42,316,852 (two-base
tolerance for coordinate conventions), extracted identically from each study's junction-level RSE.
Samples are restricted to junction depth ≥ 20, below which the ratio is too noisy to trust. The
tested contrast is the PP−NN (lesional vs normal) shift in PSI-β.

## Result: a real-but-fragile shift that rests on one cohort

**Per-study (Figure, panel a).** Of the studies clearing the depth floor, only the anchor cohort
SRP035988 shows a significant β shift (+1.17 pp, Wilcoxon p = 0.017). The deep, well-powered
SRP165679 shows essentially none (−0.07 pp, p = 0.70); the small cohorts SRP126422 (+1.01 pp,
p = 0.86) and SRP154474 (+2.21 pp, p = 0.37) are underpowered and non-significant individually.

**Estimator dependence (panel b).** Whether the pooled shift is "significant" depends entirely on
the estimator chosen:

| Method | Estimate | p | Verdict |
|---|---|---|---|
| Meta — DL random-effects | +0.87 pp | 0.018 | significant |
| Stouffer weighted-Z (√n) | — | 0.045 | significant |
| Stouffer, unweighted | — | 0.124 | n.s. |
| Meta — HKSJ (few-study honest) | +0.80 pp | 0.132 | n.s. |
| Fisher (unsigned) | — | 0.193 | n.s. |
| Mega — depth-weighted linear | +0.40 pp | 0.211 | n.s. |
| Mega — quasibinomial, 268 samples, study covariate | +0.46 pp | 0.228 | n.s. |

Only the two estimators that up-weight the anchor and use the anti-conservative normal quantile
(DL + z; √n-weighted Stouffer) clear p = 0.05. Every estimator that either propagates few-study
uncertainty (HKSJ), weights studies equally (unweighted Stouffer), or pools at the sample level
(both mega-analyses on all 268 samples) returns non-significant. The HKSJ caveat is the same one
that applies to the STAT3 *gene*-level forest (§11.5a): with k≈3–4 studies and high heterogeneity,
the DL z-interval under-covers, and the t-based HKSJ interval is the honest primary.

**Leave-one-out (panel c).** The dependence localises to a single study. Dropping any one cohort
from the √n-weighted Stouffer combination:

| Dropped | Combined p |
|---|---|
| drop SRP165679 | 0.011 |
| drop SRP126422 | 0.045 |
| drop SRP154474 | 0.066 |
| **drop SRP035988 (anchor)** | **0.913** |

Removing the anchor moves the combined p from ~0.045 to **0.91**. The entire β signal is one
cohort's effect.

**Gradient contrast (panel d).** Along the disease gradient NN→PN→PP in SRP165679, PSI-β is flat
(8.85 → 7.69 → 8.78 %) while STAT3 *gene* expression climbs monotonically (0 → +0.16 → +1.21
log₂FC vs NN). The psoriasis STAT3 signal is transcriptional up-regulation of the whole gene, not a
redistribution toward the β isoform.

## Interpretation (surgical, in both directions)

- It does **not** overturn STAT3 the gene or STAT3 the pathway. Those are confirmed and
  strengthened by the same meta-analysis: STAT3 protein is being made *more* in lesional skin
  (gene-level direction robust to leave-one-out including anchor removal; §11.5a).
- It **does** retire the specific claim that psoriasis shifts the α/β *isoform ratio* toward the
  dominant-negative β. The anchor's p = 0.017 is most parsimoniously a study-specific effect — real
  in that cohort, or a depth/batch artifact — not a reproducible feature of the disease. Under
  few-study-honest inference (HKSJ) and at the sample level (mega-analysis), the shift is not
  significant, and it vanishes entirely when the anchor is removed.

## The measurement ceiling

This is near the limit of what bulk short-read RNA-seq can establish about STAT3 isoforms. γ and δ
(neutrophil/granulocyte cleavage products) and ε/ζ (no consensus transcript) are not measurable
here at any depth. Resolving them would require a different data type: targeted long-read RNA-seq
or high-depth junction capture for a definitive α/β ratio, and Ribo-seq or isoform-resolved
proteomics for the cleavage products γ/δ. Absent that, the honest prior for any future isoform
work is the one above: **the α/β ratio in psoriasis is not reproducibly shifted; the STAT3 signal
is gene-level up-regulation.**
