# References (durable library)

Per paper: full citation + DOI, why it is in the project, and the SPECIFIC
claim/parameter/method taken from it. Organised by the role each plays in the
argument, not by date.

---

## 1. Methods & data sources

### Sun, Guan, Moran, … Xia 2022 — Scissor
- **DOI:** 10.1038/s41587-021-01091-3 · *Nature Biotechnology* 40, 527–538
- **Why:** the core method — identifying phenotype-associated single-cell
  subpopulations by integrating bulk + single-cell data.
- **Used:** network-regularized elastic-net objective (reimplemented on glmnet);
  Gaussian/continuous phenotype mode for the ordinal NN→PN→PP gradient.

### Wilks et al. 2021 — recount3
- **DOI:** 10.1186/s13059-021-02533-6 · *Genome Biology* 22, 323
- **Why:** uniform reprocessing (Monorail) of all bulk SRA studies used.
- **Used:** `create_rse_manual` gene-level counts, gencode v26, for all
  meta-analysis studies. Note: recount3 is **bulk-only** — no single-cell.

### DerSimonian & Laird 1986 — random-effects meta-analysis
- **DOI:** 10.1016/0197-2456(86)90046-2 · *Control Clin Trials* 7(3):177–188
- **Why:** the pooling estimator for all bulk meta-analysis.
- **Used:** closed-form method-of-moments τ² (`meta_dl`), vectorised across
  ~27,000 genes; verified identical to `metafor::rma(method="DL")` to machine
  precision (max Δ logFC 2.7×10⁻¹⁵; see `meta_dl_vs_metafor_concordance.csv`).
  Small-k caveat → HKSJ sensitivity analysis (deck S23).

### Ma et al. 2023 — GSE173706 single-cell psoriasis atlas
- **Accession:** GEO GSE173706
- **Why:** the single-cell reference for the Scissor arm and the deconvolution
  signature.
- **Used:** 10x droplet scRNA-seq of NN/PN/PP skin (96,088 cells → 89,058 after
  QC). Platform limit: 10x short-read **cannot resolve STAT3 isoforms per cell**.

### Tsoi et al. 2019 — SRP165679 bulk RNA-seq
- **Accession:** recount3 SRP165679
- **Why:** tier-balanced bulk phenotype anchor (Scissor arm) + one of the k=4
  meta-analysis studies.
- **Used:** NN=38 / PN=27 / PP=28 biopsies; the **only** recount3 study carrying
  both healthy and uninvolved controls (anchors the control-type moderator).

---

## 2. The vascular thesis — external support for the project's central finding

*Our result: endothelium is the top phenotype-tracking lineage (OR 5.24) and
rises monotonically across tiers (q=1.3×10⁻¹⁹). These papers independently
support a vascular-first reading of psoriasis.*

### Heidenreich, Röcken & Ghoreschi 2009 — Angiogenesis drives psoriasis pathogenesis
- *Int. J. Experimental Pathology* 90(3):232–248 · DOI 10.1111/j.1365-2613.2009.00669.x
- **Why:** the review-level statement of the vascular thesis; directly frames our
  endothelial finding as mechanism rather than curiosity.
- **Used / key claims:**
  - "Formation of new blood vessels **starts with early psoriatic changes and
    disappears with disease clearance**" — vascular change is early *and*
    reversible, i.e. it tracks disease activity (exactly what our gradient
    analysis measures).
  - "**Normalization of the superficial microvascular dermal plexus proceeds
    normalization of the epidermal structure**" (Braverman & Sibley 1982) — the
    vasculature resolves *before* the epidermis.
  - "Vascular proliferation **precedes** changes in epidermal keratin" (the
    'active edge' of plaque psoriasis).
  - Normal arterial capillary loops in dermal papillae convert to a **venous
    phenotype** in psoriatic skin — consistent with our top gradient genes being
    post-capillary venule markers (**ACKR1** +2.64, **PLVAP** +2.44).
  - Pro-angiogenic mediators enriched in psoriatic skin: VEGF, HIF,
    angiopoietins, TNF, TGF-α, **IL-8**, IL-17 — our bulk k=4 has CXCL8/IL-8
    **+4.78** (FDR 0.04).

### Garshick et al. 2020 — Activated platelets induce endothelial inflammation via COX-1
- *Arterioscler. Thromb. Vasc. Biol.* 40(5):1340–1351 · DOI 10.1161/ATVBAHA.119.314008 · PMID 32131611
- **Why:** supplies the **causal, interventional** arm our correlative data lacks,
  and an external driver of endothelial inflammation.
- **Used / key claims:**
  - Platelets from psoriasis patients (n=45 vs 18 controls) show increased
    activation **correlating with skin severity**; 2–3× increased adhesion to
    human aortic endothelial cells.
  - Platelets induce endothelial pro-inflammatory transcripts — **IL8, IL1β,
    COX-2** — i.e. a circulating driver of the endothelial IL-1 biology we see.
  - Platelet RNA-seq: interferon signature; **COX-1 correlates with severity
    (r=0.83, p=0.01)**.
  - **Randomised trial:** 2 weeks low-dose aspirin (81 mg) reduced serum TxB2 and
    cut brachial-vein endothelial pro-inflammatory transcripts **>70%** vs no
    treatment (NCT03228017) — the vascular compartment is *pharmacologically
    addressable*.
  - Our overlap: CXCL8/IL-8 **+4.78** (FDR 0.04) and TBXAS1 (thromboxane
    synthase) **+0.45** (FDR 0.017) in bulk k=4. Caveat: COX-1/COX-2 and platelet
    markers are flat in whole skin — expected, since their measurements are in
    blood and brachial-vein endothelium.

### Ponikowska et al. 2026 — Cardiovascular Disease and Psoriasis (review)
- *Dermatol Ther (Heidelb)* 2026;16:155–169 · DOI 10.1007/s13555-025-01566-0
- **Why:** the current, citable synthesis of the psoriasis→CV link; supplies the
  Background claim in the conference abstract and independently names the same
  mechanistic genes our data flags.
- **Used / key claims:**
  - Psoriasis carries an **~50% increased relative risk of major cardiovascular
    events**; a 2022 meta-analysis (31 cohorts, >665,000 patients vs ~18M
    controls) gives pooled RR ≈1.17 myocardial infarction, 1.19 stroke, **1.46
    cardiovascular mortality**.
  - Endothelial cells isolated from psoriasis patients show **2–8-fold
    upregulation of pro-inflammatory transcripts**, resembling endothelium
    stimulated by TNF-α + IL-17A + IFN-γ — a direct vascular impact of psoriatic
    inflammation. (Independent support for our endothelium-first result.)
  - **Inflammasome signalling is the most differentially expressed systemic
    pathway in psoriasis** and plays a causal role in atherosclerosis — directly
    relevant to our IL1R1⁺ endothelial subset being GSDMD/PYCARD/CASP1-high.
  - S100A7A and S100A8/A9 alarmins correlate with psoriasis severity and with
    atherosclerotic plaque instability.
  - Biologics (TNF-α, IL-17, IL-23 inhibitors) associate with **~50% reductions
    in myocardial infarction incidence** vs topical regimens.
- **Concordance with our bulk k=4 meta-analysis** (all lesional vs healthy):
  S100A7A **+7.78**, S100A7 **+6.37**, S100A9 **+6.08**, S100A8 **+6.01**,
  CXCL10 **+3.01**, PYCARD **+1.09**, GSDMD **+0.75**, CASP1 **+0.59** — every
  gene the review names is significantly up in our data (VCAM1 and NLRP3 are the
  two exceptions, both n.s.).

### → Concordance: the review's "emerging mediators" vs our k=4 meta-analysis

Ponikowska 2026 §Emerging Molecular Pathways nominates three CV-linked mediators
beyond TNF-α/IL-17/IFN-γ. **All three are detectable in our bulk meta-analysis,
and two rank among the top effects genome-wide** — and notably, the strongest of
them is IL-1 superfamily. (Full table: `results/il1_superfamily_bulk_k4.csv`.)

| Review's mediator | Our k=4 result (lesional vs healthy) |
|---|---|
| **IL-36 family** (IL-1 superfamily; keratinocyte activation, vascular injury) | **IL36A +9.03** (FDR 1.8×10⁻¹⁹) — **rank 2 of 28,339 genes** by effect size (behind DEFB4A +9.84); IL36G +4.66, IL36RN +2.79, IL36B +1.05 (all ***) |
| **S100 alarmins** (S100A7A, S100A8/A9; severity + plaque instability) | **S100A7A +7.78**, S100A12 +7.44, S100A7 +6.37, S100A9 +6.08, S100A8 +6.01 (all ***) |
| **GPR15** (epithelial immunity, lymphocyte trafficking) | **+1.84** (FDR 0.028 *) — significant but weakest: k=2, I²=86%. Treat as suggestive, not established. `GPR15L` not in annotation. |

**Caveat:** this concordance is **tissue-level (bulk)**. IL36A/IL36B/S100A9/
S100A12/GPR15 are absent from the single-cell tested family (below expression
thresholds), and IL36G/S100A8 are *lower* in disease-associated cells — a
composition effect, since those cells are endothelium-enriched while IL-36/S100
are keratinocyte genes. The concordance does not extend to the cell-level arm.

### The IL-1 superfamily in psoriasis — our full profile

Because IL-36 is IL-1 superfamily, the project's IL-1 story is best read across
the whole family. Our k=4 meta-analysis separates it cleanly into an **active
arm and a silent arm**:

**Strongly up — the IL-36 axis and inflammasome effectors**
- IL36A **+9.03** (rank 2 genome-wide), IL36G **+4.66**, IL36RN **+2.79**, IL36B **+1.05** (all ***)
- CASP5 **+3.54**, AIM2 **+2.84**, PYCARD **+1.09**, CASP4 **+0.88**,
  GSDMD **+0.75**, CASP1 **+0.59** (all *** / **)
- IL1RN (antagonist) **+0.83** ***; IL33 +0.50 ***

**Not significant / down — the classical IL-1β axis**
- **IL1B +1.09 n.s.** (I²=94 — among the most heterogeneous in the panel, with IL18 and IL37 at I²=95)
- IL18 −0.05 n.s.; IL1A **−0.47** **; IL37 **−1.91** *
- Receptors: IL1R1 −0.02 **n.s.**, IL1R2 −0.34 n.s., IL1RAP +0.20 n.s.
  (IL1RL2, the IL-36 receptor, is up +0.29 **)
- NLRP3 +0.52 n.s.

**Reading:** in whole skin it is the **IL-36 sub-family, not IL-1β, that carries
the IL-1 superfamily signal**, together with the downstream inflammasome
effectors (caspases, PYCARD, GSDMD, AIM2). The classical IL-1β ligand/receptor
axis is flat at tissue level. This is consistent with (i) the disappointing
performance of IL-1β blockade in plaque psoriasis, (ii) DITRA — IL-36 receptor
antagonist deficiency causing severe pustular disease, and (iii) our cell-level
finding that IL-1 *responsiveness* (IL1R1⁺ endothelium) is an early, declining
feature rather than a maintained one.

### STAT3 and the JAK–STAT axis — and why it is NOT IL-1 superfamily

**Classification (important, easy to conflate):** STAT3 is **not** a member of
the IL-1 superfamily. The IL-1 superfamily comprises secreted ligands
(IL-1α/β, IL-18, IL-33, IL-36α/β/γ, IL-37, IL-38) and their receptors —
extracellular signalling proteins. **STAT3 is an intracellular transcription
factor**, phosphorylated by JAK kinases downstream of **gp130-family** receptors
(IL-6, IL-11, IL-22, IL-23, LIF, OSM).

**The real link is indirect and downstream:**
`IL-1 → IL1R1 → MyD88 → NF-κB → IL-6 → IL6R/gp130 → JAK → STAT3`
IL-1 can therefore *drive* STAT3, but through an IL-6 intermediate — they are
different molecule classes on the same pathway, not the same family.

**Our k=4 profile** (full table: `results/jakstat_axis_bulk_k4.csv`):

| Layer | Result |
|---|---|
| STATs | **STAT1 +1.70***, **STAT3 +1.06***, STAT2 +0.80***; STAT5B −0.25***, STAT4/5A n.s. |
| JAKs | **JAK3 +1.12*** (I²=95), TYK2 +0.23***; JAK1 **−0.23*** |
| gp130 cytokines (true STAT3 upstream) | **IL23A +1.50***; IL6 −1.37 **n.s.** (I²=99, k=2), IL6R/IL6ST/IL22RA1 n.s., **LIF −1.50*** |
| Negative regulators | SOCS3 +0.86 **n.s.** (I²=98), SOCS1 n.s., PTPN2 +0.46** |

**Reading:** STAT3 mRNA is reproducibly elevated (+1.06, all four cohorts), but
its canonical upstream ligands are **not** — IL6 is non-significant with extreme
heterogeneity (I²=99), IL6R/IL6ST flat, LIF *down*. The one significant upstream
signal is **IL23A (+1.50)**, consistent with the IL-23/Th17 axis rather than
IL-6 driving STAT3 here. Combined with our single-cell result (STAT3 not
significant per-cell, q=0.13; 66% vs 40% enriched specifically in IL1R1⁺
endothelium), the tissue-level STAT3 elevation looks **compositional and
IL-23-associated**, not a cell-intrinsic IL-6→STAT3 program.

*Caution for writing:* do not describe STAT3 as IL-1 superfamily. The
defensible phrasings are "STAT3, downstream of IL-1-induced IL-6" or "the
JAK–STAT axis," and our IL1R1⁺-endothelium result is precisely where the two
pathways intersect at cell level.

### Psoriatic march / endothelial dysfunction → atherosclerosis
- *Endothelial Dysfunction in Psoriasis: An Updated Review* — Front. Med. 2022;9:864185
- *Psoriasis and Cardiovascular Disease: Novel Mechanisms* — PMC9744099
- **Used:** the ~50% increased CV risk; psoriatic endothelium adopts a
  pro-inflammatory, pro-atherogenic phenotype (IL-1, IL-6, TNF-α, VCAM1/ICAM1/
  E-selectin). Frames the abstract's closing significance statement.

---

## 3. The IL-1 axis — and the direction question

*Our result: within endothelium, IL1R1 is **highest in healthy skin** and falls
with progression (51%→29%→27%, NN→PN→PP; log2FC −1.12), while IL1RN (antagonist)
rises (+4.10 endothelium; +0.83 bulk). Interpretation: IL-1 initiates rather than
maintains vascular involvement.*

### Cooper et al. 1990 — IL-1 activity in psoriatic skin
- *J. Invest. Dermatol.* (attached scan)
- **Why:** the historical precedent for the counter-intuitive direction we found.
- **Used:** IL-1 **bioactivity is reduced** in lesional psoriatic skin despite
  active inflammation, with a shift toward IL-1 inhibitory/antagonist activity —
  a decades-old observation that anticipates our IL1R1↓ / IL1RN↑ result. This is
  the citation that turns our finding from "contradicts the literature" into
  "recapitulates a known paradox at cell-type resolution."

### Cai et al. 2019 — A critical role of the IL-1β–IL-1R pathway in psoriasis
- *J. Invest. Dermatol.* · DOI 10.1016/j.jid.2018.07.025 · PMC6392027 · PMID 30120937
- **Why:** the standard citation for IL-1R being *required* in psoriasis; must be
  read carefully so it is not mistaken for a contradiction of our result.
- **Used / precise scope:** IL-1β is elevated in lesional skin; IL-1R **pathway
  activity** (GSEA) tracks glucocorticoid response and recurrence in patients;
  IL-1R on **γδ T cells and keratinocytes** is required for IMQ-induced
  inflammation (mouse KO).
  ⚠ **Not** a measurement of IL1R1 transcript in **endothelium** — different
  compartment, pathway-activity vs transcript abundance, mouse-necessity vs
  human-abundance. No conflict with our finding.

### IL-36 / early IL-23
- *IL-36 signaling in keratinocytes controls early IL-23 production* —
  Life Sci. Alliance · DOI 10.26508/lsa.202000688 · PMC7190273
- **Used:** supports "IL-1-family as an early signal"; and the IL-36 axis is the
  dominant significant IL-1-family signal in our bulk data (IL36A +9.03,
  IL36RN +2.79) — not IL1B (n.s., I²=94%).

### S1PR3–STAT3 feedback
- *Cell Death & Disease* · s41419-025-07358-w
- **Used:** early + prolonged STAT3 activation in psoriatic keratinocytes;
  context for STAT3 being reproducible in bulk (+1.06) yet not a per-cell
  gradient marker (q=0.13).

---

## 4. Theory 2 — mosaic loss of Y (deferred; needs blood)

- *Mosaic loss of chromosome Y and cardiovascular disease* — Nat. Rev. Cardiol.
  2023 · DOI 10.1038/s41569-023-00976-x
- *Hematopoietic mLOY → cardiac fibrosis and heart failure* — Science 2022 ·
  DOI 10.1126/science.abn3100
- **Used:** grounds theory 2 (mLOY worse in males). mLOY is hematopoietic →
  cannot be called from skin biopsy; our skin Y-expression check was a clean
  negative at the dropout floor (`notes/lossofY_skin_check.md`).

---

## Cross-references
- Arm-A detail: `IL1B_STAT3_expanded_analysis.md`
- Arm-C detail: `theory1_endothelial_IL1_and_blood_arm.md`
- Method validation: `meta_dl_vs_metafor_concordance.csv`
- Study selection: `study_exclusion_table.csv`, `all_human_studies_sample_metadata.csv`
