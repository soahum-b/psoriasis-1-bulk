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
5. **Interpretation** — cell-type enrichment (Fisher) of Scissor± vs background; gradient-program DE; orthogonal NNLS bulk deconvolution concordance.

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
