# A molecular staging axis for psoriasis: healthy -> peri-lesional -> lesional skin as an ordered transcriptional trajectory

*Companion paper to "Meta-analysis of the psoriasis transcriptome centred on STAT3"
(psoriasis_meta_analysis_whitepaper.md). That paper establishes, across five recount3
cohorts, that STAT3 is a reproducible central hub of the lesional psoriasis programme. This
paper reframes the same data around a different question: not "what differs between lesional
and healthy skin," but "in what order does the psoriasis programme assemble," using
peri-lesional skin as an observable intermediate stage.*

---

## 0. Summary

Psoriasis transcriptomics is almost always framed as a two-group contrast: lesional (PP)
versus healthy (NN) skin. But many cohorts also sample **peri-lesional** skin (PN) -
clinically uninvolved skin from the same patients - and that third group turns a binary
comparison into an **ordered axis**: NN -> PN -> PP.

Analysing the deepest three-group recount3 cohort (SRP165679; NN=38, PN=27, PP=28), we show
that the psoriasis transcriptional programme is overwhelmingly **monotonic** along this axis.
Of 4,563 lesional differentially expressed (DE) genes, **85% are strictly monotonic**
(NN<PN<PP or NN>PN>PP) and **99% show a clear monotonic rank trend** (|Spearman rho|>0.3
between expression and stage). Peri-lesional skin sits, on average, only **16% of the way**
from healthy to lesional - a genuine early molecular stage, not a binary switch.

Placing genes on this axis by *how much of the lesional change has already happened at the
PN stage* separates the programme into **early** and **late** modules. The result is a clear
temporal ordering: **the interferon / IL6-JAK-STAT3 / inflammatory programme is the early
event** (already significantly elevated in peri-lesional skin), while **keratinocyte
cell-cycle / proliferation is the late event** (flat at PN, switching on only at the lesional
stage). Inflammation precedes proliferation.

The canonical psoriasis hubs - STAT3, STAT1, SOCS3, IL17A, IL23A and the antimicrobial
keratinocyte genes - all rise monotonically along the axis, confirming that the staging
framework recovers established biology. STAT3 therefore enters this paper as **validation**,
not as the finding; and the axis adds one mechanistic detail a pairwise contrast cannot see:
STAT3's inferred regulon *activity* reaches ~30% of its lesional level in peri-lesional skin
while its own *transcript* reaches only ~10%, consistent with STAT3 being activated through
its targets before it is strongly transcribed - placing STAT3 activation early on the disease
trajectory.

---

## 1. Why a staging axis (and why peri-lesional skin is the key)

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

## 2. The staging axis is real: 85% of the psoriasis programme is monotonic

We used the deepest three-group cohort in our recount3 selection, SRP165679 (the Tsoi
atopic-dermatitis/psoriasis cohort [10]), keeping its psoriasis and healthy arms: NN=38
healthy, PN=27 peri-lesional, PP=28 lesional. Expression was the uniformly reprocessed,
normalised log-expression used throughout the companion meta-analysis.

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

![Figure S5. Molecular staging axis. Top 40 progressively-up and top 40 progressively-down monotonic genes (rows), samples ordered healthy (NN, blue) -> peri-lesional (PN, orange) -> lesional (PP, red). Up-genes transition blue->pale->red across the axis and down-genes do the reverse; the peri-lesional block is visibly intermediate.]({{artifact:art_1d7a9b1d-a5a8-4879-8d07-f4fb97b29efd}})

![Figure S5b. Quantifying the axis. (A) 85% of lesional DE genes are strictly monotonic and 99% show a monotonic rank trend across NN<PN<PP. (B) Distribution of the peri-lesional position (fraction of the lesional change already reached at the PN stage): most genes cluster low, median 16%, showing peri-lesional skin is an early molecular stage.]({{artifact:art_5e16b6b5-0a6e-4b3b-adf4-93b6deeb0b98}})

The key methodological point: prior work frames this as three separate pairwise DEG lists
and observes qualitatively that PN "looks intermediate." Testing the ordering *directly*
turns that impression into a number - 85% monotonic, PN at 16% - with a per-gene position on
the axis that the next sections exploit.

## 3. A timing taxonomy: early, progressive, and late genes

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

![Figure S2. Gradient-gene timing taxonomy. (A) Class sizes. (B) Hallmark enrichment of the progressive (early/inflammatory) versus late (proliferation) classes: progressive genes are dominated by interferon-gamma/alpha, IL6-JAK-STAT3 and inflammatory/TNF-NF-kB signalling; late genes are dominated by G2M checkpoint and E2F targets.]({{artifact:art_a2802b70-d309-4d26-a15c-69bbe9c66c63}})

The progressive class - the genes already rising in peri-lesional skin - is an **immune /
interferon** programme (interferon-gamma FDR 7e-24, interferon-alpha 9e-23, IL6-JAK-STAT3
2.5e-5, inflammatory response and TNF-NF-kB). The late class - the genes that appear only in
established lesions - is a **proliferation** programme (G2M checkpoint FDR 1.8e-8, E2F
targets 2.4e-6). This is the central biological claim of the paper, and section 4 confirms
it at the pathway-score level.

## 4. Pathway timing: inflammation is early, proliferation is late

To confirm the taxonomy independently of any gene-classification threshold, we computed
per-sample single-sample GSEA (ssGSEA) scores for the key Hallmark programmes in all 93
samples and tracked each programme's mean score across the three stages.

![Figure S3. Pathway activation timing. Mean ssGSEA score (+/- SE) across NN -> PN -> PP for immune/interferon programmes (left) and proliferation programmes (right). Immune programmes lift clearly at the peri-lesional stage; proliferation programmes stay near baseline until the lesional stage.]({{artifact:art_9d639ec8-5300-4226-827f-eb1af9d153f7}})

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

## 5. STAT3 and the IL-17 axis: validation of the staging framework

The companion paper's central result is that STAT3 is a reproducible network hub of the
lesional programme. Here we place the canonical hubs on the staging axis - as a check that
the framework recovers known biology, and to locate STAT3 in the disease sequence.

![Figure S4. Canonical psoriasis hubs on the staging axis. STAT3, STAT1, SOCS3, IL17A, IL23A, S100A7, DEFB4A and IL36G all rise monotonically NN -> PN -> PP, each sitting intermediate at the peri-lesional stage. Every trend is significant (FDR <= 6e-10).]({{artifact:art_4014f5c6-fbf3-4822-8c98-f743145922c7}})

All eight hubs rise monotonically, each intermediate at peri-lesional - the axis reproduces
the established psoriasis biology gene by gene. STAT3 thus enters as validation, not as the
discovery.

The axis also adds one detail a pairwise contrast cannot. Inferring STAT3 transcription-
factor activity (decoupleR ULM on the CollecTRI regulon, relative to the healthy baseline)
across the three stages:

![Figure S4b. STAT3 transcription-factor activity on the staging axis. STAT3 activity rises monotonically (NN 0.01 -> PN 2.21 -> PP 7.43; trend p = 4e-16), reaching ~30% of the lesional level already at peri-lesional.]({{artifact:art_25cf5b00-633a-437e-9561-72ade6e1b720}})

STAT3 *activity* reaches ~30% of its lesional level at peri-lesional, while STAT3's own
*transcript* reaches only ~10% (section 2 coordinate). The gap is mechanistically coherent:
STAT3 acts through phosphorylation and target induction, so its regulon can light up in
uninvolved skin before the gene itself is strongly transcribed. On the staging axis, STAT3
activation is an early event.

## 6. What this adds to the existing non-lesional-skin literature

Two ideas here already exist, and we claim neither as novel. Peri-lesional skin as a
molecular intermediate is well established - Gudjonsson et al. described decreased lipid
biosynthesis and increased innate immunity in uninvolved skin [1]; a proteomic study is
explicitly titled "Intermediate Stage of Non-Lesional Psoriatic Skin" [2]; reviews describe
uninvolved skin as "somewhere in between" [3]; and single-cell and spatial studies reproduce
the intermediate profile [4,5]. Cross-study integration of psoriasis transcriptomes also
exists [6,7,8].

Against that backdrop, the specific contributions of this analysis are:

1. **An explicit ordering statistic on a three-point axis**, replacing three pairwise DEG
   lists and a qualitative "PN sits between" with a direct trend + monotonicity test (85%
   monotonic, PN at a median 16%), and a per-gene axis coordinate.
2. **A gradient-timing taxonomy** that orders the disease programmes - inflammation early,
   proliferation late - a sharper statement than "both are dysregulated in lesions."
3. **STAT3 positioned on the axis** with the activity-before-transcription signature.
4. **A reproducible, recount3-anchored pipeline** (companion pipeline STEP 12).

**Limitations, stated plainly.** The clean three-point axis is well-powered in one cohort
(SRP165679); other three-group cohorts in our selection are shallow or tiny, so the axis is
*anchored* in one deep study and *supported* by pooled trends rather than independently
replicated at full power. Peri-lesional effect sizes are small (median ~16% of the lesional
change), so the story is strongest at the pathway/trend level. And the axis is a
cross-sectional ordering of three sampled states, not a temporal trajectory in individual
patients.

## 7. Conclusion

Adding the peri-lesional group turns the psoriasis case-control transcriptome into an ordered
molecular staging axis. Most of the psoriasis programme moves monotonically along it, with
peri-lesional skin an early intermediate; the interferon/STAT inflammatory programme is the
early event and keratinocyte proliferation the late one; and the canonical STAT3/IL-17 hubs
validate the axis while STAT3's activity marks it early. The framework reframes a well-known
lesional signature as a *sequence*, and it does so from public, uniformly reprocessed data
with a reproducible pipeline.

---

## Methods (brief)

- **Cohort.** SRP165679 (recount3, Gencode v26, GRCh38), psoriasis + healthy arms: NN=38,
  PN=27, PP=28. Normalised log-expression as in the companion meta-analysis.
- **Trend test.** Per-gene limma moderated-t on numeric stage (NN=0,PN=1,PP=2); Spearman of
  expression vs stage. Monotonicity by strict ordering of group means. PN fraction =
  (mu_PN - mu_NN)/(mu_PP - mu_NN).
- **Taxonomy.** Classes by PN fraction: early (>=0.50), progressive (0.15-0.50), late
  (<0.15), PN-divergent (<0). Hallmark enrichment by Fisher's exact test against the
  trend-tested gene universe.
- **Pathway timing.** Base-R ssGSEA (Barbie 2009) per sample for key Hallmark sets; linear
  trend and PN-vs-NN Wilcoxon per pathway.
- **STAT3 activity.** decoupleR ULM on the CollecTRI regulon, per sample relative to the
  healthy mean.
- **Reproducibility.** All steps in psoriasis_pipeline.R STEP 12; figures and tables in
  staging_axis/.

## References

See staging_axis/literature_positioning.md for the full annotated reference list [1]-[10].
Key: [1] Gudjonsson 2009 JID; [2] Intermediate-stage proteomics PMC6684579; [3] PMC3347524;
[9] Li 2014 JID (SRP035988 anchor); [10] Tsoi 2019 JID (SRP165679, staging cohort).
