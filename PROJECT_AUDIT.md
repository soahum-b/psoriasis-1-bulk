# Project audit — psoriasis meta-analysis → STAT3 / pathways → targets

*Independent audit, 2026-07-17. Scope: all project artifacts + three granted local
folders + the cluster git repo `soahum-b/psoriasis-1-bulk` (CSB node n013). This document
records the state of each analytical arm, what is strong, what is missing, and a prioritised
list of concrete next actions. It is a companion to CHECKPOINT.md (resume-from-here state)
and DECISIONS_LOG.md (choice history); it does not replace them.*

---

## 1. Executive summary

The project is a mature, multi-arm transcriptomic dissection of psoriasis built on recount3
whole-skin bulk RNA-seq plus one single-cell reference (GSE173706, Ma 2023). Its central,
well-supported result is that **STAT3 is a reproducible hub of the lesional psoriasis
programme at the gene and pathway level**, embedded in an IL-17 / interferon / IL6-JAK-STAT3
inflammatory core with downstream S100 / defensin / chemokine output.

The work's defining strength is **statistical honesty**: random-effects meta-analysis with
few-study-corrected intervals, retained negative results (the STAT3 α/β isoform switch and
the single-cell STAT3-as-gradient-marker claim were both correctly retired), and thorough
decision logging.

The thinnest part relative to the biology is the **target-selection layer** — the step from
"these pathways/proteins are up" to "these are the ranked, tractable, druggable nodes." That
is the largest scientific opportunity. The largest operational item is **figure/state
housekeeping**, most of which is already closer to done than the notes suggest.

---

## 2. Arm-by-arm status

### Arm A — Bulk meta-analysis (backbone) · SOLID
- 3→4-study random-effects (DerSimonian-Laird) meta-analysis of lesional (PP) vs normal (NN),
  harmonised on Gencode v26 via recount3. Contrasts PP-NN, PN-NN, PN-PP.
- STAT3 pools to **+1.06 (k=4) / +1.25 (k=3)** log2FC, robust to leave-one-out including
  anchor removal. Direction is not anchor-driven.
- Heterogeneity handled correctly: HKSJ (t-based, few-study-honest) reported as primary at
  low k, DL as sensitivity. I² interpreted as *magnitude* not *direction* disagreement.
- Pathway (GSEA/ORA/CAMERA/ROAST), TF-activity (decoupleR/CollecTRI), and a cross-study
  co-expression module (clust) all converge on the IL-17/IFN/IL6-JAK-STAT3 + proliferation
  programme.
- Assembled into `psoriasis_integrated_whitepaper.md` (13 sections + staging axis).
- **Verdict:** publication-grade. No methodological gap identified.

### Arm B — Staging axis (NN→PN→PP as ordered trajectory) · SOLID
- 85% of the lesional programme is monotonic across the three stages; peri-lesional skin sits
  a median 16% of the way to lesional — an early molecular stage, not a binary switch.
- Timing taxonomy: **inflammation/interferon is early** (already elevated in peri-lesional),
  **keratinocyte proliferation is late** (switches on only at the lesional step).
- STAT3 *activity* (regulon) reaches ~30% of lesional level at PN while its *transcript*
  reaches only ~10% — activity-before-transcription, placing STAT3 activation early.
- **Verdict:** a genuine analytical contribution over the standard pairwise contrast.
  Limitation correctly stated: well-powered in one cohort (SRP165679), supported by pooled
  trends elsewhere.

### Arm C — Scissor single-cell arm · SOLID, with honest self-correction
- Ordinal biopsy-site phenotype projected onto 89,058 cells; non-circular by construction
  (clinical label, not molecular score).
- **Endothelial cells dominate** the gradient-tracking population (OR 5.24 at full census,
  BH q≈0), validated by two permutation nulls (reliability p=0.000, selection null p=0.000)
  and orthogonal NNLS deconvolution (3/3 strongest lineages concordant).
- STAT3 as a *gradient marker* was **correctly downgraded to n.s.** at full census
  (log2FC 0.30, BH q=0.13) after being significant on the 20k backbone — a model example of
  not over-claiming a small-subset effect.
- **Verdict:** structurally sound. See §4 for the one open housekeeping item (figures).

### Arm D — Theory-1: IL1R1⁺ endothelium / psoriatic march · SOLID as an exploratory result
- DC→endothelial paracrine IL-1 model: IL-1β is DC-restricted (48% of DCs express it, 0.2%
  of endothelium), but endothelium carries the receptor (IL1R1 ~29%).
- IL1R1⁺ endothelial cells are 2.3× enriched in the gradient-tracking set and co-express more
  STAT3/IL6/NFKB1/CASP1/GSDMD/PYCARD — the receptor→program wiring made visible at cell level.
- Receptor is early/constitutive (IL1R1 down in lesional endothelium), consistent with IL-1 as
  an initiating signal.
- **Verdict:** a specific, literature-grounded angle within the established psoriatic-march
  theme. Appropriately framed as exploratory.

### Dropped / deferred (correctly)
- **Blood bulk arm — dropped.** ERP110814-vs-GTEx is perfectly study-batch confounded
  (disease status 100% aliased with lab/globin-prep). Correct call.
- **mLOY (loss of Y) — deferred.** Genomic mLOY needs blood/DNA; an RNA Y-expression proxy in
  skin was a clean negative (dropout-floor artifact, no disease trend).
- **STAT3 isoform (α/β) switch — retired.** Real only in the anchor cohort; vanishes under
  HKSJ, mega-analysis, and anchor-drop. Near the ceiling of what short-read bulk can resolve.

---

## 3. Prioritised gaps and opportunities

Ordered by expected scientific impact on the stated goal (genes → pathways → **proteins/targets**).

### P1 — Systematic target/druggability layer *(largest scientific gap)*
The project converges beautifully on pathways but stops short of a ranked target table. Missing:
- **Open Targets** query of the meta-analysis hit list: genetic association (psoriasis GWAS),
  tractability bucket (small-molecule / antibody), known/approved drugs, safety.
- **ChEMBL / DGIdb** drug-gene interaction mapping for the top pathway proteins.
- An explicit **abundance ≠ validity** treatment: TYK2 is flat at the transcript level
  (+0.12, n.s.) yet deucravacitinib (TYK2i) is an approved psoriasis drug. JAK3 is strongly up
  (+2.04) and drugged (tofacitinib); JAK1 is *down* (−0.32). The target argument must be made
  on pathway position + tractability, not transcript fold-change alone.
- **Deliverable:** a ranked target table (protein · pathway role · meta log2FC · gradient
  timing · genetic evidence · tractability · existing drugs).

### P2 — Orthogonal, drug-first target evidence (CMap / LINCS)
Test whether the psoriasis meta-signature is **reversed** by JAK/TYK2 inhibitors and other
perturbagens. This nominates/prioritises targets from perturbation data independently of
expression, and cross-checks P1.

### P3 — Protein / activation-state evidence for STAT3
Every conclusion rests on mRNA. The staging axis already hints STAT3 acts via activity before
transcription. Strengthen with a literature-anchored pSTAT3 statement, and (if any proteomic
psoriasis dataset is reachable) a protein-level cross-check.

### P4 — Sex-stratified skin re-analysis *(queued, feasible now)*
Sex is already called per donor (21M/12F). A sex-stratified re-run of the skin arm connects to
the male-severity thread and is the natural next increment with existing data. No new data
needed.

### P5 — Reproducibility housekeeping
- **Git push:** RESOLVED on inspection — cluster `main` (ee40082) == GitHub `origin/main`
  (public, pushed 2026-07-17). The CHECKPOINT "blocked on push credentials" note is stale and
  should be corrected. *(This audit corrects it.)*
- **§4 figure regen:** the whitepaper §4 figures still depict the 20k backbone; regenerate from
  the 89k full-census object (`results_full/reference_scissor_full.rds`, cluster).
- **Compiled-Scissor (APML1) cross-check:** still staged; validates the glmnet solver port.
- **Benchmarked deconvolution (deconvBenchmarking):** staged Stage-1 upgrade to replace the
  transparent-but-unbenchmarked NNLS.

### P6 — Close the circulation / mLOY thread (data acquisition)
The blood arm and mLOY can only close with a **single study containing both plaque-psoriasis
and healthy blood processed together** (not available in recount3 — needs GEO/ArrayExpress).
Until then both threads stay open. A recount3 title-keyword scan for plaque-psoriasis PBMC is
a noted open task before concluding none exists.

---

## 4. Things confirmed already handled well (no action needed)
- Random-effects model choice and few-study interval correction.
- Study eligibility screen (23-study recount3 audit; defensible inclusions/exclusions).
- Scissor non-circularity argument.
- Retention of negative results rather than burial.
- Living decision log with scope tags.

---

## 5. One-line recommendation
The science is strong and honest; the missing half-step is turning the converged pathway
biology into a **ranked, tractable target table** (P1) with an orthogonal drug-first
cross-check (P2). Everything else is either already done (git) or incremental
(figures, sex split).

---

## 6. Detailed walkthrough — how the project unfolds, step by step

*This section traces the analysis in the order the story actually developed, so the audit is
self-contained. It condenses `NARRATIVE_walkthrough.md` (the full guided tour) and
`METHODS.md` (the reproducibility record); consult those for every intermediate number and
pin. Each step names the input, the operation, the output artifact, and the finding.*

### The question and the arc
The stated goal is **genes → pathways → proteins/targets**, with STAT3 as the lead
hypothesis. What the data delivered is more nuanced than a single-gene confirmation, and the
intellectual arc is the value:

```
STAT3 as lead target (hypothesis)
   │
   Arm A (bulk):        STAT3/JAK-STAT robustly up  ──► looks like a target ✓
   │                    IL1B unreliable; IL-36 dominant; inflammasome up
   ▼
   Arm B (single-cell): STAT3 does NOT track progression at full census ✗
   │                    the progression-tracking cells are ENDOTHELIAL (vascular-led)
   ▼
   Arm C (theory 1):    the vascular program is an IL-1-RESPONSIVE endothelium
   │                    IL1R1⁺ vessels run STAT3 + inflammasome (paracrine, DC-sourced IL-1β)
   ▼
   Reframed target:     not STAT3 alone, not IL-1β alone —
                        the IL-1-primed vascular circuit (STAT3/inflammasome downstream),
                        with a direct line to the psoriatic-march cardiovascular comorbidity
```

### Step 1 — Data source and harmonisation
- **Input:** recount3 (uniformly reprocessed SRA, Gencode v26 / GRCh38) for bulk; GSE173706
  (Ma 2023, 33 count matrices, 22 donors) for single-cell.
- **Why recount3:** uniform annotation harmonises the gene space by construction — no
  cross-platform probe mapping, the historic bane of expression meta-analysis.
- **Output:** harmonised per-study count matrices; a 96,088→89,058-cell (post-QC) Seurat
  reference annotated to 9 lineages.

### Step 2 — Study eligibility screen (Arm A foundation)
- **Operation:** 23 human recount3 psoriasis studies screened on sample-level metadata;
  eligibility = whole-skin biopsy, contains lesional (PP) + normal (NN), treatment-free,
  adequate depth.
- **Output:** `study_eligibility_PPvsNN.csv`. Final **k=4**: SRP035988, SRP165679, SRP126422,
  SRP065812 (141 NN + 145 PP). ERP110816 excluded (all on etanercept, which suppresses the
  exact TNF→IL-1→STAT3 axis under study); in-vitro/cell-line/blood studies excluded.

### Step 3 — Per-study differential expression
- **Operation:** identical pipeline per study — `filterByExpr` → TMM (edgeR) → voom → lmFit →
  eBayes (limma). Effect = log2FC, SE = log2FC/t.
- **Rationale:** model-then-merge (meta-), not merge-then-model (mega-), so each study keeps
  its own normalization and dispersion and heterogeneity stays inspectable.

### Step 4 — Random-effects meta-analysis (the STAT3 gene verdict)
- **Operation:** closed-form DerSimonian-Laird random-effects pooling, vectorised across
  ~24k genes, validated to 4 dp against `metafor::rma`. Heterogeneity via Q / I² / τ²;
  few-study interval via Hartung-Knapp-Sidik-Jonkman (primary), DL as sensitivity.
- **Output:** `meta_de_PPvsNN_4study.csv`; `fig25_stat3_forest.png`; `stat3_forest_robustness.csv`.
- **Findings:**
  - **STAT3/JAK-STAT is the tightest signal** — STAT3 +1.06, STAT1 +1.70, JAK3 +1.12,
    TYK2 +0.23 (all sig, k=4); JAK1 −0.23 (constitutive). STAT3 gene-level up-regulation
    survives leave-one-out including anchor removal (+1.12 without the big anchor).
  - **IL1B is the least reliable panel member** (+1.09, n.s., I²=94%; the 4th cohort pulls
    it down) — a well-replicated negative. The inflammasome machinery around it (CASP1,
    PYCARD, IL1RN, AIM2, CASP5) is significant.
  - **IL-36 dwarfs IL-1β** (IL36A +9.0, IL36G +4.7); the shared downstream output
    (S100A7/8/9, DEFB4A, LCN2, CCL20) is the strongest signal overall.
  - Adding the 4th study shrank inflated fold-changes toward honest pooled values and exposed
    (did not create) between-study disagreement in SOCS3/IL18/CXCL2.

### Step 5 — Pathway, TF-activity, and co-expression (genes → pathways)
- **Operations:** GSEA + ORA (`gsea_ora_merged.csv`), CAMERA/ROAST competitive gene-set tests
  (`meta_pathway_camera_PPvsNN.csv`), decoupleR TF-activity on CollecTRI
  (`tf_activity_collectri.csv`), and cross-study co-expression via **clust** (consensus
  modules replicating across cohorts).
- **Findings:** convergent IL-17 / interferon-α,γ / IL6-JAK-STAT3 / TNF-NF-κB inflammatory
  core plus a keratinocyte-proliferation cassette (E2F, G2M). clust returns one high-confidence
  65-gene module fusing the antimicrobial/IL-17 keratinocyte programme with a mitotic cassette
  — the two plaque hallmarks recovered as one co-regulated unit with no disease labels used.

### Step 6 — STAT3 isoform question (retired, correctly)
- **Operation:** PSI-β = β-junction/(α+β) from junction-level RSEs (depth ≥ 20), PP−NN shift
  pooled across cohorts, seven estimators compared.
- **Output:** `stat3_isoform_ab_formal_section.md`; `fig40_stat3_isoform_formal_evaluation.png`.
- **Finding:** the anchor's apparent β-switch (p=0.017) **does not replicate** — it vanishes
  under HKSJ, unweighted Stouffer, and both mega-analyses, and disappears entirely on
  anchor-drop (combined p 0.045 → 0.91). The psoriasis STAT3 signal is **gene-level
  up-regulation, not an isoform redistribution**. Near the ceiling of what short-read bulk
  can resolve.

### Step 7 — Staging axis (Arm B foundation: order of assembly)
- **Operation:** in the deepest three-group cohort (SRP165679; NN=38/PN=27/PP=28), fit a
  linear-trend model + Spearman monotonicity per gene across NN=0/PN=1/PP=2; classify by
  "fraction of lesional change already reached at PN."
- **Findings:** 85% of the lesional programme is strictly monotonic; PN sits a median 16% of
  the way to lesional. **Inflammation/interferon is early** (elevated already at PN),
  **proliferation is late** (flat at PN, switches on only at PP). STAT3 *activity* leads its
  *transcript* (~30% vs ~10% of lesional level at PN).

### Step 8 — Scissor single-cell projection (pathways → cells)
- **Operation:** project the ordinal biopsy-site phenotype (NN<PN<PP) onto 89,058 cells via a
  network-regularized elastic net (symmetric-normalized graph Laplacian, sparse edge-difference
  augmentation, solved on glmnet — a validated pure-R re-implementation of Scissor's APML1).
  Non-circular by construction (clinical label, not a molecular score).
- **Significance controls:** reliability test (100 label permutations, p=0.000) + selection
  permutation null (30 reruns, p=0.000).
- **Findings:** selected cells are monotonic on the gradient and peak at the peri-lesional
  tier; **endothelial cells dominate** (OR 5.24 at full census); the program is vascular-led
  (CCL14, ACKR1, RAMP3, PLVAP). **STAT3 does not survive as a gradient marker** at full census
  (BH q=0.13) — the 20k-backbone significance was a small-subset effect. The bulk arm's
  "STAT3 up" and the single-cell arm's "STAT3 doesn't track progression" are reconciled:
  STAT3 is expressed nearly everywhere and is a passenger in a vascular-led program.

### Step 9 — Orthogonal deconvolution (composition vs state)
- **Operation:** NNLS deconvolution of the real SRP165679 bulk against a Ma-reference
  signature; test monotonic proportion trends NN→PN→PP.
- **Finding:** the 3 strongest lineages are concordant with Scissor direction — endothelial
  rises (0→0.2→3.8%), fibroblast and melanocyte fall — confirming the endothelial signal is a
  **real compositional increase**, not solely a state change. Discordant cases (NK
  infiltration, keratinocyte dominance) are mechanistically informative.

### Step 10 — Theory-1: the IL-1-primed vascular circuit (cells → mechanism → target)
- **Operation:** on the full-census object, localise IL-1β source vs receptor; within
  endothelium, Fisher-test IL1R1⁺ vs IL1R1⁻ against Scissor+; compare downstream co-expression.
- **Output:** `theory1_endothelial_IL1_and_blood_arm.md`; `fig_theory1_endothelial_IL1.png`.
- **Findings:** IL-1β is **DC-restricted** (48% of DCs, 0.2% of endothelium) → a **paracrine
  model**. IL1R1⁺ endothelial cells are **2.3× enriched** in the gradient-tracking set
  (OR 2.32, p=1.4×10⁻³⁸) and co-express markedly more STAT3/IL6/NFKB1/CASP1/GSDMD/PYCARD —
  the receptor→program wiring made visible. Receptor is early/constitutive (down in lesional
  endothelium), consistent with IL-1 as an initiating signal. Ties directly to the
  psoriatic-march cardiovascular comorbidity.

### Step 11 — Side investigations (negatives worth keeping)
- **Blood arm — dropped:** ERP110814-vs-GTEx is perfectly study-batch confounded (disease
  100% aliased with lab/globin-prep). No covariate fixes a confounder collinear with the
  exposure.
- **Loss-of-Y in skin — clean negative:** Y-silent-cell fraction is depth-driven (31.5%→1.4%)
  and does not track disease tier. Genomic mLOY needs blood/DNA.

### Where each arm's target-relevant output lands
| Arm | Target-relevant readout | Status |
|---|---|---|
| A (bulk meta) | STAT3/JAK-STAT up, druggable via JAK/TYK2i; IL-36 > IL-1β | solid |
| B (Scissor) | endothelial/vascular program is the progression driver | solid |
| C (theory-1) | IL-1-responsive IL1R1⁺ endothelium as the integrated target | exploratory, grounded |
| Staging | inflammation early, proliferation late; STAT3 activation early | solid |

**Reframed target statement:** not STAT3 alone and not IL-1β alone, but the **IL-1-primed
vascular circuit** with STAT3/inflammasome activity downstream — which is also where the
missing systematic druggability pass (§3, P1) would land its ranking.
