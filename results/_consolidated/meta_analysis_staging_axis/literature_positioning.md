# Literature positioning: what a cross-study monotonic staging axis adds

## The prior art is real - we position against it honestly

Two ideas central to this reframe already exist in the psoriasis literature, and the
contribution here is **not** to claim either as novel:

**1. Peri-lesional/uninvolved skin is a molecular intermediate.** This is well established.
Gudjonsson et al. (2009) showed uninvolved psoriatic skin carries a distinct signature -
decreased lipid biosynthesis and increased innate immunity - relative to healthy skin,
with 179 dysregulated genes [1]. A proteomic study is titled, explicitly, "Intermediate
Stage of Non-Lesional Psoriatic Skin" and shows most differential proteins place
non-lesional skin between healthy and lesional, while flagging proteins that break the
trend [2]. A review states plainly that uninvolved skin "is somewhere in between the skin
of the healthy individual and the inflamed skin" [3]. Recent single-cell/lncRNA work
reproduces the intermediate profile of PN skin [4], and spatial transcriptomics shows even
distal non-lesional skin stratifies by disease severity [5].

**2. Cross-study integration of psoriasis transcriptomes exists.** Cross-study concordance
of the PP-vs-PN signature has been quantified [6], and a 23-dataset integration built
co-expression network models of lesional and uninvolved skin [7]. Suarez-Farinas et al.
evaluated the psoriasis transcriptome across studies by GSEA [8].

## What is genuinely added here

Against that backdrop, the specific, defensible contributions of this analysis are:

**(a) An explicit ordering statistic on a three-point axis, not a pairwise contrast.**
Prior work almost universally frames the problem as *pairwise* DEG lists (PP-vs-NN,
PN-vs-NN, PP-vs-PN) and then notes, qualitatively, that PN "sits between." Here the
staging structure is tested *directly*: a per-gene linear-trend model across NN=0, PN=1,
PP=2 plus a rank-based monotonicity check (Spearman of expression vs stage). This converts
"PN looks intermediate" into a quantitative statement - **85% of lesional DE genes are
strictly monotonic (NN<PN<PP), 99% show a clear monotonic rank trend, and peri-lesional
skin sits a median 16% of the way to lesional** - with a per-gene position on the axis.

**(b) A gradient-timing taxonomy that separates early from late programs.** Because each
gene gets a "fraction-of-lesional-reached-at-PN" coordinate, the psoriasis program can be
partitioned into progressive (rising through PN), late/PP-specific (switching on only at
the lesional step), and PN-specific classes. This yields a concrete, testable ordering of
disease programs: **the interferon / IL6-JAK-STAT3 / inflammatory axis is the early event
(already significantly elevated in peri-lesional skin), while keratinocyte
cell-cycle/proliferation is the late event (flat at PN, switching on only at the lesional
stage).** This "inflammation precedes proliferation" ordering, recovered from a
cross-sectional three-group transcriptome, is a sharper statement than the standard
observation that both programs are dysregulated in lesions.

**(c) STAT3 positioned on the axis as validation, with an activation-before-transcription
signature.** Placing the canonical hubs on the staging axis shows the framework recovers
known biology (all rise monotonically). It also surfaces a mechanistically coherent detail
that pairwise contrasts miss: STAT3's inferred regulon *activity* reaches ~30% of its
lesional level in peri-lesional skin, while its own *transcript* reaches only ~10% -
consistent with STAT3 being activated (through its targets) before it is strongly
transcribed, and placing STAT3 activation early on the disease trajectory.

**(d) A reproducible, recount3-anchored pipeline.** The staging analysis is anchored in a
uniformly reprocessed recount3 cohort (SRP165679, the deepest three-group study: NN=38,
PN=27, PP=28) and supported by pooled trends across the wider five-study meta-analysis, all
from a single reproducible pipeline.

## Honest limitations

- The clean three-point NN<PN<PP axis is well-powered in **one** cohort (SRP165679).
  SRP076982 has peri-lesional samples but no healthy baseline and is shallow; SRP126422 is
  tiny. The staging axis is therefore *anchored* in one deep cohort and *supported* by
  pooled trends, not independently replicated at full power in three cohorts.
- Peri-lesional effect sizes are small (median ~16% of the lesional change), so the story
  is strongest at the pathway/genome-wide trend level; per-gene PN-stage significance is
  modest for most individual genes.
- The axis is a *cross-sectional* ordering (three sampled states), not a true temporal
  trajectory; "early/late" refers to position along the healthy->lesional axis, not to
  time in an individual patient.

## References

[1] Gudjonsson JE, Ding J, Li X, et al. Global gene expression analysis reveals evidence
for decreased lipid biosynthesis and increased innate immunity in uninvolved psoriatic
skin. J Invest Dermatol. 2009;129(12):2795-2804. doi:10.1038/jid.2009.173.

[2] Bata-Csorgo Z et al. Comprehensive Proteomic Analysis Reveals Intermediate Stage of
Non-Lesional Psoriatic Skin and Points out the Importance of Proteins Outside this Trend.
(PMC6684579).

[3] Study of Molecular Mechanisms Involved in the Pathogenesis of Immune-Mediated
Inflammatory Diseases, using Psoriasis As a Model. (PMC3347524).

[4] Profiling Long Noncoding RNA in Psoriatic Skin Using Single-Cell RNA Sequencing.
J Invest Dermatol. 2024. (builds on Tsoi et al. 2019, the SRP165679 cohort).

[5] Spatial transcriptomics stratifies psoriatic disease severity by emergent cellular
ecosystems. (PMC10502701).

[6] Cross-Study Homogeneity of Psoriasis Gene Expression in Skin across a Large Expression
Range. PLoS One. 2013. (PMC3537625).

[7] The integration of large-scale public data and network analysis uncovers molecular
characteristics of psoriasis. Hum Genomics. 2022;16:doi:10.1186/s40246-022-00431-x.

[8] Suarez-Farinas M, Lowes MA, Zaba LC, Krueger JG. Evaluation of the psoriasis
transcriptome across different studies by gene set enrichment analysis (GSEA). PLoS One.
2010;5(4):e10247.

[9] Li B, Tsoi LC, Swindell WR, et al. Transcriptome analysis of psoriasis in a large case-
control sample: RNA-seq provides insights into disease mechanisms. J Invest Dermatol.
2014;134(7):1828-1838. (SRP035988, the anchor cohort).

[10] Tsoi LC, et al. Atopic dermatitis is an IL-13-dominant disease with greater molecular
heterogeneity, contrasted with the IL-17 axis of psoriasis. J Invest Dermatol. 2019.
(SRP165679, the three-group cohort anchoring the staging axis).
