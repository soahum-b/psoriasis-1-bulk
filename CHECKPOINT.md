# CHECKPOINT (resume-from-here)

Lightweight state snapshot — overwrite each checkpoint. Large serialized state lives in the
artifact store / `results_full/`, referenced by ID below.

## As of 2026-07-16
- **Running now:** nothing for this project. (Unrelated `md_prod` GPU array 56887578 is a different project.)
- **Last completed:** Arm-A bulk k=4 IL-1β/STAT3 meta-analysis + writeup; project-setup scaffolding (additive) into this repo.
- **Next:** retrieve + validate Scissor full-census numbers vs the 20k backbone (job 56882314 already finished on n003); then update WHITEPAPER.md §4 and re-render.
- **Key IDs:**
  - Scissor full-census SLURM job: **56882314** (completed) → `results_full/`
  - Cluster repo: `~/Soahum/Project/psoriasis/psoriasis-1-bulk` (n003, git)
  - Artifacts: IL1B_STAT3_expanded_analysis.md `1ad8ca91`, k=4 forest `44e10919`, k=4 meta CSV `1ec80908`, k3-vs-k4 panel `34f3b909`, eligibility `8099b3fa`, resume note `7dd05578`
  - Full-census headline (verified from `results_full` object): **alpha=0.20, 6,837 Scissor+ / 6,228 Scissor− = 14.67% selected**; reliability CV-MSE 0.135 vs null ~0.835. See `fullcensus_vs_backbone.csv` (artifact 223ced02). Note: full-census tuning settled at alpha=0.20 (backbone used 0.40) — verified `scissor_tuning.rds$chosen$alpha`.
  - **STAT3 correction:** at full census STAT3 is n.s. in the Scissor+ program (log2FC 0.30, p=0.12, padj=1) — the backbone's "sig, padj 0.018" was a subset effect. WHITEPAPER §4 STAT3 sentence must be corrected on the next update; endothelial-led thesis is unchanged (OR 5.24, top lineage).
- **Blocked on:** git push credentials unresolved (local commits only, no push). Cluster node name changes per session — re-register.
