# Methods / Reproducibility (living)

## Arm A — Bulk multi-study PP-vs-NN meta-analysis (STAT3 / IL-1β)
1. **Study screen** — 23 human recount3 psoriasis studies screened against sample-level metadata; eligibility = contains lesional (PP) and normal (NN) skin, whole-skin biopsy, treatment-free, adequate depth. Result: `study_eligibility_PPvsNN.csv`. Final k=4: SRP035988, SRP165679, SRP126422, SRP065812 (141 NN + 145 PP).
2. **Per-study DE** — `filterByExpr` → TMM (edgeR) → voom → lmFit → eBayes (limma). Effect = log2FC, SE = log2FC/t. Applied identically per study on the harmonised gene space (`de_one`).
3. **Random-effects meta** — closed-form DerSimonian-Laird (`meta_dl`): inverse-variance weights, τ² by method-of-moments, I² heterogeneity, BH-FDR. Validated to 4 dp vs `metafor::rma`.
4. **Panel figure** — IL-1β/inflammasome + STAT3/JAK-STAT + IL-36 + Th17 + output modules; forest plot with 95% CI, per-gene k flags, I²≥75% markers.

## Arm B — Scissor-on-gradient (single-cell)
1. **Reference build** — GSE173706 33 count matrices → Seurat 5; Ensembl→symbol (org.Hs.eg.db); QC → 89,058 cells; cluster → 9 lineages.
2. **Bulk anchor** — SRP165679 ordinal phenotype (NN<PN<PP), logCPM (edgeR).
3. **Scissor solve** — network-regularized elastic net; symmetric-normalized graph Laplacian `L_sym = I − D^(−1/2) A D^(−1/2)` encoded as sparse edge-difference augmentation; solved on glmnet (`scissor_glmnet_solver.R`), λ-path walked to target selected fraction; alpha=0.40.
4. **Significance** — reliability test (100 label permutations) + selection permutation null (30 shuffled-label reruns).
5. **Interpretation** — cell-type enrichment (Fisher) of Scissor± vs background; gradient-program DE (Scissor+ vs Background, `FindMarkers`, logfc.threshold=0.1, min.pct=0.1); orthogonal NNLS bulk deconvolution concordance.

## Arm C — Theory 1: IL-1-responsive endothelium (single-cell, full census)
1. **Receptor localization** — % expressing IL1B (ligand) vs IL1R1/IL1R2/IL1RAP (receptor complex) per cell type (FetchData over the 89k object); establishes source (DC) vs responder (endothelium/fibroblast).
2. **Gradient enrichment** — within endothelium, Fisher exact of IL1R1+ vs IL1R1− against Scissor+ vs Background (OR + p).
3. **Conditional co-expression** — % expressing STAT3/IL6/NFKB1/CASP1/GSDMD/PYCARD in IL1R1+ vs IL1R1− endothelial cells.
4. **Direction** — endothelium PP vs NN (`FindMarkers`, no thresholds) over a 16-gene IL-1/adhesion/atherosclerosis set; BH-FDR over the set.
   Outputs: `figures/fig_theory1_endothelial_IL1.png`, `notes/theory1_endothelial_IL1_and_blood_arm.md`.

## Arm D — Blood / circulation (planned)
1. **recount3 blood screen** — psoriasis blood/PBMC studies; known set (SRP173379/78, SRP132160) all GPP → excluded (subtype mismatch).
2. **ERP110814 baseline** — join ArrayExpress SDRF (E-MTAB-6555) `Factor Value[time]` to recount3 `external_id` (ENA run); keep week-0 (10 tx-naive plaque-psoriasis blood). No healthy ctrl → pair with an external healthy-blood dataset (to be sourced).
3. **Sex inference** (theory 2) — per-donor Y-gene (RPS4Y1/DDX3Y/UTY/EIF1AY/KDM5D/USP9Y) vs XIST expression; sex-stratified skin analysis feasible now, mLOY needs blood.

## Multiple-testing correction (protocol — applies to ALL arms)
- **BH-FDR** (`p.adjust(method="BH")`), q < 0.05, is the significance rule throughout: bulk meta-analysis (Arm A), single-cell gradient-program DE and cell-type enrichment (Arm B), and deconvolution trend tests.
- **Seurat caveat:** `FindMarkers` returns `p_val_adj` as **Bonferroni** (p × n_genes), NOT BH. For the single-cell arm we therefore recompute `fdr_BH = p.adjust(p_val, "BH")` over the tested gene family and report that, not Seurat's `p_val_adj`. (Full-census gradient program: 3,903/4,541 genes at BH q<0.05; STAT3 raw p=0.117 → BH q=0.13, n.s.)
- Permutation-based tests (reliability, selection null) are single global statistics, not a testing family — no FDR correction applies.

## Versions / pins
- Cluster: env `scissor-r`, R 4.5, Seurat 5 (see COMPUTE.md for full list).
- Laptop bulk: env `psoriasis-r`, R 4.5.3.

## Results index
| File | What it is |
|---|---|
| `study_eligibility_PPvsNN.csv` | 23-study screen + per-study decisions |
| `meta_de_PPvsNN_4study.csv` | genome-wide k=4 meta-analysis table |
| `meta_PPvsNN_k3_vs_k4_panel.csv` | k=3 vs k=4 panel comparison (Δlog2FC, I², flips) |
| `fig_il1b_stat3_NNvsPP_4study.png` | k=4 both-modules forest plot |
| `IL1B_STAT3_expanded_analysis.md` | Arm-A writeup section |
| `results_full/reference_scissor_full.rds` | full-census Seurat object w/ Scissor selections (1.2 GB, cluster) |
| `results_full/gradient_program_DE.csv` | full-census gradient program |
| `WHITEPAPER.md` | Arm-B Scissor progression-gradient paper |
| `HANDOFF.md` | Arm-B status / reproduce / next steps |
| `results_full/gradient_program_DE_BH.csv` | full-census gradient program with BH-FDR column (4,541 genes) |
| `results_full/celltype_enrichment_BH.csv` | full-census cell-type enrichment, BH-FDR |
| `figures/fig_theory1_endothelial_IL1.png` | Arm-C: IL1R1+ endothelium enrichment + downstream program + PP-vs-NN panel |
| `notes/theory1_endothelial_IL1_and_blood_arm.md` | Arm-C writeup + blood-arm feasibility |
| `results/blood_arm_feasibility_recount3.csv` | recount3 blood screen (GPP-only) |
| `code/regen_fig4_fullcensus.R` | regenerates WHITEPAPER §4 figures from the 89k full-census object → `figures_full/` |
| `figures_full/*.png` | full-census §4 figures (alpha tuning, Scissor UMAP, composition, significance controls, gradient-program volcano, STAT3) — replace the 20k-backbone figures |
| `PROJECT_AUDIT.md` | independent project audit: arm status, prioritised gaps, step-by-step walkthrough |
