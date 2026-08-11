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
