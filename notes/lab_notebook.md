# Lab Notebook (dated, append-only)

Chronological record of what was tried, what happened, and WHICH artifact/figure it produced.

## 2026-07-13/14 — Scissor full-census cluster run
- **Did:** submitted `code/run_full_census.sbatch` (SLURM job 56882314) on CSB, full 89,058-cell reference.
- **Saw:** job ran to completion (log reaches `perm 100/100`). Final saved selection (from `results_full` object, `table(scissor)`): **6,837 Scissor+ / 6,228 Scissor− = 14.67% selected** (backbone was 14.84%). NOTE: this differs from the 6,422/5,981 alpha=0.40 tuning-grid log line I first mis-quoted; the object is authoritative. The exact alpha the final run selected is NOT yet verified (tuning.rds read was blocked by a cluster-proxy outage) — confirm from `scissor_tuning.rds` before citing an alpha. Reliability real CV-MSE 0.135 vs null ~0.835. NOTE: STAT3 does NOT survive at full census (log2FC 0.30, p=0.12, padj=1; 42.7% vs 44.9% expressing) — the backbone's "STAT3 sig, padj 0.018" was a subset effect; the program stays vascular-led (endothelial OR 5.24, top lineage).
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
