# CHECKPOINT (resume-from-here)

Lightweight state snapshot — overwrite each checkpoint. Large serialized state lives in the
artifact store / `results_full/`, referenced by ID below.

## As of 2026-07-17 (audit + reproducibility sync)
- **Independent project audit produced** → `PROJECT_AUDIT.md` (+ PDF): arm-by-arm status,
  prioritised gaps (P1 = systematic target/druggability table is the thinnest layer), and a
  step-by-step walkthrough. Copied to local `Recount-3/PROJECT_AUDIT.md`.
- **Git push: NOT blocked (stale note corrected).** Cluster `main` (ee40082) == GitHub
  `origin/main` (public, pushed 2026-07-17). The three "local-only" commits (e470eca, 0aa3093,
  31b6a3c) are ancestors of the pushed HEAD — already on GitHub. The earlier "blocked on push
  credentials" note was stale.
- **§4 figures regenerated from the 89k full-census object** (DONE). `regen_fig4_fullcensus.R`
  (SLURM job 56889747, n011) → `figures_full/{alpha_tuning, scissor_umap, scissor_composition,
  significance, gradient_program, stat3}.png`. All six saved as artifact v2 over the backbone
  figures, so WHITEPAPER §4 inline figures now resolve to full-census. Selection-null gap now
  0.382 vs null 0.019; STAT3 panel carries the n.s. (BH q=0.13) verdict.
- **Stray lowercase `checkpoint.md` removed** (superseded by this CHECKPOINT.md).
- **Still-open (unchanged):** systematic target/druggability table (P1); CMap/LINCS reversal;
  sex-stratified skin pass; APML1 cross-check; benchmarked deconvolution; Sei-LLRA arm; blood arm
  (needs non-recount3 internally-controlled dataset).

## As of 2026-07-17 (pm — pausing, resume tomorrow)
- **Blood bulk arm: DROPPED** (ERP110814-vs-GTEx study-batch confounded). **Loss-of-Y in skin: clean NEGATIVE** (RNA proxy at dropout floor, no disease trend; genomic mLOY needs blood). Both documented + committed (5e56b4a + this doc sync).
- **Resume tomorrow with:** Theory-2 Part A — sex-stratified analysis of the SKIN data (feasible now, no blood; sex already called: 21M/12F donors). Prompt drafted: `PROMPT_theory2_sex_mLOY.md` (artifact f202208b) if run in a fresh chat.
- **Defensible arms standing:** A (bulk k=4 PP-vs-NN meta), B (Scissor gradient), C (Theory-1 IL1R1+ endothelium). All skin-based.
- **Still-open (from HANDOFF):** WHITEPAPER HTML/PDF re-render; §4 figures regen from 89k object; APML1 cross-check; benchmarked deconvolution (Stage 1); Sei-LLRA sequence arm.

## (earlier 2026-07-17)
- **Running now:** nothing for this project. (Unrelated `md_prod` GPU array 56887578 is a different project.)
- **Last completed:** full-census validation + protocol BH-FDR audit (STAT3 q=0.13 n.s.; WHITEPAPER §4 corrected); Theory-1 endothelial IL-1 analysis (IL1R1+ endothelium = gradient-tracking/STAT3-inflammasome-active, OR 2.32); blood-arm scoping (GPP excluded; ERP110814 baseline usable; needs healthy-blood ctrl); all .md docs synced.
- **Next:** (1) include ERP110814 wk-0 baseline (SDRF→recount3 join) + source a healthy-blood control; (2) screen GEO/CELLxGENE for plaque-psoriasis PBMC single-cell; (3) sex-stratified skin pass (theory 2 down payment); (4) re-render WHITEPAPER HTML/PDF; (5) regenerate §4 figures from 89k object.
- **Key IDs:**
  - Scissor full-census SLURM job: **56882314** (completed) → `results_full/`
  - Cluster repo: `~/Soahum/Project/psoriasis/psoriasis-1-bulk` (n003, git)
  - Artifacts: IL1B_STAT3_expanded_analysis.md `1ad8ca91`, k=4 forest `44e10919`, k=4 meta CSV `1ec80908`, k3-vs-k4 panel `34f3b909`, eligibility `8099b3fa`, resume note `7dd05578`
  - Full-census headline (verified from `results_full` object): **alpha=0.20, 6,837 Scissor+ / 6,228 Scissor− = 14.67% selected**; reliability CV-MSE 0.135 vs null ~0.835. See `fullcensus_vs_backbone.csv` (artifact 223ced02). Note: full-census tuning settled at alpha=0.20 (backbone used 0.40) — verified `scissor_tuning.rds$chosen$alpha`.
  - **STAT3 status (DONE):** at full census STAT3 is n.s. under protocol BH-FDR (log2FC 0.30, raw p=0.12, **BH q=0.13**; Seurat's padj=1 was Bonferroni). WHITEPAPER §4 corrected (v2, 63b78805); METHODS codifies BH rule. Endothelial-led thesis unchanged (OR 5.24, top lineage). Gradient program 3,903/4,541 at BH<0.05.
  - **Theory-1 artifacts:** fig `93ad7fe6`, writeup `1a3e3584`, blood feasibility `45450d34`, BH tables `77c061df`/`29ba858c`.
  - **Blood arm:** GPP studies excluded (user); ERP110814 = plaque-psoriasis blood, wk-0 recoverable via ArrayExpress SDRF E-MTAB-6555 `Factor Value[time]` joined to recount3 ENA-run ids (10 tx-naive PP, no healthy ctrl). ERP110816 skin baseline has no NN → PP-vs-PN sensitivity only.
  - **Cluster commits (local):** e470eca (BH-FDR), 0aa3093 (theory-1), 31b6a3c (blood-claim softening).
- **Blocked on:** ~~git push credentials unresolved~~ **RESOLVED** — cluster `main` == GitHub
  `origin/main` (ee40082, public). Nothing unpushed. Cluster node name still changes per session
  — re-register each session.
