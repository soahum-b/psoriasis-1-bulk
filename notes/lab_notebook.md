# Lab Notebook (dated, append-only)

Chronological record of what was tried, what happened, and WHICH artifact/figure it produced.

## 2026-07-13/14 — Scissor full-census cluster run
- **Did:** submitted `code/run_full_census.sbatch` (SLURM job 56882314) on CSB, full 89,058-cell reference.
- **Saw:** job ran to completion (log reaches `perm 100/100`). Final saved selection (from `results_full` object, `table(scissor)`): **6,837 Scissor+ / 6,228 Scissor− = 14.67% selected** (backbone was 14.84%). Verified `scissor_tuning.rds$chosen$alpha = 0.20` — the full-census tuning re-ran and settled at alpha=0.20 (the 20k backbone used 0.40); the 0.20 grid row is exactly 6,837/6,228/14.67%. (My first note mis-quoted the 6,422/5,981 alpha=0.40 grid line as the selection.) Reliability real CV-MSE 0.135 vs null ~0.835. NOTE: STAT3 does NOT survive at full census (log2FC 0.30, p=0.12, padj=1; 42.7% vs 44.9% expressing) — the backbone's "STAT3 sig, padj 0.018" was a subset effect; the program stays vascular-led (endothelial OR 5.24, top lineage).
- **Produced:** `results_full/{reference_scissor_full.rds (1.2 GB), gradient_program_DE.csv, reliability_test.rds, permutation_null.rds}`.
- **Refs:** Sun 2022 (Scissor).

## 2026-07-16 — Bulk IL-1β/STAT3 + k=4 expansion (this session)
- **Did:** screened 23 recount3 studies; added SRP065812 (16 NN + 18 pre-adalimumab PP) → k=4; excluded ERP110816 (anti-TNF, no timepoint field). Re-fit per-study DE (identical pipeline) + DL meta.
- **Saw:** STAT3 stays +1.06 FDR<0.001 (k=4); IL1B n.s. (+1.09, I²=94%), pulled down by new cohort; IL-36 dominant (IL36A +9.0). 3 significance flips (SOCS3, IL18, CXCL2) — all I²=95–98% = genuine between-study disagreement exposed.
- **Produced:** artifacts `meta_de_PPvsNN_4study.csv` (1ec80908), `meta_PPvsNN_k3_vs_k4_panel.csv` (34f3b909), `fig_il1b_stat3_NNvsPP_4study.png` (44e10919), `study_eligibility_PPvsNN.csv` (8099b3fa), `IL1B_STAT3_expanded_analysis.md` (1ad8ca91).
- **Refs:** recount3, Tsoi 2019.

## 2026-07-16 — Isoform feasibility + project scaffolding
- **Did:** confirmed GSE173706 is 10x droplet → isoform-by-cell-type infeasible per-cell. Applied project-setup skill (additive) to the cluster repo.
- **Saw:** full-census results already present on n003; no re-run needed.
- **Produced:** `RESUME_il1b_stat3_session.md` (7dd05578); docs/ + notes/ + cluster/ skeleton (this scaffolding).

## 2026-07-16/17 — Full-census validation, BH-FDR audit, theory 1, blood-arm scoping
- **Did:** validated full-census object; applied protocol BH-FDR to single-cell arm (was Seurat Bonferroni); corrected WHITEPAPER §4 STAT3 to n.s.; ran theory-1 (IL-1 receptor localization, IL1R1×Scissor enrichment, endothelium PP-vs-NN); screened recount3 blood (GPP-only, excluded); examined ERP110814/ERP110816 baseline via ArrayExpress SDRF; sex-inference feasibility.
- **Saw:** STAT3 BH q=0.13 (n.s.); gradient program 3,903/4,541 at BH<0.05. Theory-1: IL1R1+ endothelium is the gradient-tracking, STAT3/inflammasome-active population (Fisher OR 2.32, p 1.4e-38); IL1R1 early-high (down in lesional, log2FC −1.12). Blood: known recount3 psoriasis-blood = GPP; ERP110814 = plaque-psoriasis blood, wk-0 baseline recoverable via SDRF `Factor Value[time]` (corrects earlier "no timepoint field") but no healthy ctrl. Sex cleanly inferable (Y-genes vs XIST, ~20M/~10F).
- **Produced:** `results_full/gradient_program_DE_BH.csv` (77c061df), `celltype_enrichment_BH.csv` (29ba858c), WHITEPAPER.md v2 (63b78805), METHODS.md v2 (76c95368), `fig_theory1_endothelial_IL1.png` (93ad7fe6), `theory1_endothelial_IL1_and_blood_arm.md` (1a3e3584), `blood_arm_feasibility_recount3.csv` (45450d34). Cluster commits e470eca, 0aa3093, 31b6a3c.
- **Refs:** psoriatic march (Front Med 2022;864185); IL-1 early / IL-1R progression (PMC6392027); S1PR3-STAT3 (s41419-025-07358-w); mLOY-CV (Nat Rev Cardiol 2023; Science abn3100).

## 2026-07-17 (pm) — Blood arm dropped; loss-of-Y skin check (negative)
- **Blood arm:** ingested ERP110814 baseline (10 tx-naive plaque-psoriasis blood, wk-0 via ArrayExpress SDRF E-MTAB-6555). Only in-recount3 healthy control is GTEx whole blood → perfectly study-batch confounded (disease aliased with lab/prep; globin 0.05% vs 37%, ~10yr age gap). **DROPPED** per user (uninterpretable cross-study case/control). Commit 5e56b4a.
- **Loss-of-Y (RNA proxy) in skin:** 21 male donors / 65,524 male cells. Y-silent fraction is depth-driven (31.5%→1.4% low→high depth = dropout) and does not track disease tier. Clean negative; genomic mLOY still needs blood/DNA. See `notes/lossofY_skin_check.md`.
