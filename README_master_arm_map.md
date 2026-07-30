# psoriasis-1-bulk — master project map

*Reproducible-project consolidation map. Generated from a full lineage audit of all 300 project artifacts across 13 sessions. Companion to `PROJECT_AUDIT.md` (scientific state), `REPRODUCIBILITY_GAP.md` (per-artifact code recoverability), `NARRATIVE_walkthrough.md` (guided tour), and `DECISIONS_LOG.md` (choice history).*

## What this project is

A multi-arm transcriptomic dissection of psoriasis (genes → pathways → druggable proteins), built on recount3 whole-skin bulk RNA-seq plus one single-cell reference (GSE173706). Central result: a **vascular/endothelial IL-1-responsive program** dominates lesional progression; STAT3 is a reproducible bulk hub but a *passenger* at single-cell resolution. Hypothesis-generating; not clinical guidance.

## Reproducibility status

- **264/300 artifacts (88%) have recoverable generating code** in lineage; **13/14 arms fully reproducible**.
- Toolchains: `psoriasis-r` (recount3/DESeq2/edgeR/limma), `python`, `scissor-r` (Seurat5/glmnet), `dock-md` (rdkit/biopython).
- Environment manifests to capture per arm: see table below.

## Arm map

| Arm | Description | Artifacts | Repro | Toolchain | Already in repo | Stage target (additive) |
|-----|-------------|----------:|-------|-----------|-----------------|-------------------------|
| A1_meta_analysis | Bulk meta-analysis (k=4 PP-vs-NN, DerSimonian-Laird) | 101 | FULL | psoriasis-r, python | code/run_pipeline.R, 00_download_data.R; results/, figures/ | results/meta_analysis/, figures/meta_analysis/ |
| A2_peer_review | Peer review of methodology | 4 | FULL | psoriasis-r | — (artifact only) | docs/peer_review/ |
| A3_pathway_TF | Pathway / TF-activity (GSEA/ORA/CollecTRI) | 14 | FULL | psoriasis-r, python | — (artifact only) | code/pathway_TF/, results/pathway_TF/, figures/pathway_TF/ |
| A4_study_exclusion | Study-exclusion documentation (S1–S9) | 23 | FULL | psoriasis-r | (partly in docs) | results/study_exclusion/ |
| B_scissor | Scissor single-cell progression-gradient | 44 | FULL | scissor-r, python | code/scissor_*.R, run_full_census*; results_full/, figures_full/ | (present) |
| C_stat3_isoforms | STAT3 six-isoform (α/β/γ/δ/ε/ζ) characterization | 3 | FULL | python | — (artifact only) | code/isoforms/, results/isoforms/ |
| D_docking_md | Molecular docking + MD simulation | 50 | FULL | dock-md, python | — (artifact only) | code/docking_md/, results/docking_md/, figures/docking_md/ |
| E_genome_browser | Genome-browser comparison tool | 12 | FULL | psoriasis-r, python | — (artifact only) | code/genome_browser/, results/genome_browser/ |
| F_presentation | Figure presentation compile | 1 | FULL | python | — (artifact only) | docs/presentation/ |
| G_audit | Independent project audit | 2 | PARTIAL | — | PROJECT_AUDIT.md | (present) |
| H_session_current | Theory-1 IL-1/endothelium + session docs | 29 | FULL | python | docs/, notes/ (theory1, LOY, decisions, methods) | results/theory1_il1/, figures/theory1_il1/ |
| I_lncRNA | lncRNA differential expression | 11 | FULL | psoriasis-r, python | — (artifact only) | code/lncRNA/, results/lncRNA/, figures/lncRNA/ |
| J_theory2_plan | Theory-2 sex-stratified plan | 1 | FULL | — | docs/plans/ (pending) | docs/plans/ |
| K_refs_deconv | Reference PDFs + deconvolution reading | 5 | FULL | — | — (artifact only) | docs/references/ |

## How to reproduce each arm

Every code-bearing artifact carries its generating code in lineage (recover via `host.lineage[version_id]['code']`). For each arm, the per-arm REPRODUCE note (to be staged under each results/ subdir) lists: input artifacts, the script(s), the conda env, and the output artifacts. The audit CSV `recoverability_audit_full.csv` is the authoritative index (columns: arm, filename, vid, content_type, mb, class, code_len, env).

## Consolidation status (as of this session)

- **Phase 0 (audit): DONE** — `REPRODUCIBILITY_GAP.md` + `recoverability_audit_full.csv`.
- **Phase 1 (stage each arm into repo): PENDING** — blocked on cluster access (n013 SSH channel unresponsive this session).
- **Phase 2 (commit + push): PENDING** — same block.
- Resume: restore n013 (or any CSB node), then stage each arm per the table above (additive only — never overwrite existing files), write per-arm REPRODUCE notes, commit, and push to `github.com/soahum-b/psoriasis-1-bulk`.
