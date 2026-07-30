# Reproducibility gap report — psoriasis project

*Generated 2026-07-30. Audit of all 300 latest-version artifacts across 14 arms (13 sessions). For each artifact, lineage was extracted and classified by whether generating code is recoverable.*

## Headline

- **264/300 artifacts (88%) have recoverable generating code** in lineage.
- 30 authored documents (markdown writeups/READMEs — no generating code by nature).
- 5 raw inputs (user-uploaded reference PDFs).
- 1 terminal output (a rendered audit PDF) with no recoverable source.
- **13 of 14 arms are FULLY reproducible**; the audit arm has one terminal PDF (its .md source is present).

## Per-arm verdict

| Arm | Description | Artifacts | Code-recoverable | Toolchain | Verdict |
|-----|-------------|----------:|-----------------:|-----------|---------|
| A1_meta_analysis | Bulk meta-analysis (k=4 PP-vs-NN) | 101 | 92 | psoriasis-r(73), python(18) | FULLY reproducible |
| A2_peer_review | Peer review of methodology | 4 | 2 | psoriasis-r(2) | FULLY reproducible |
| A3_pathway_TF | Pathway / TF-activity | 14 | 13 | psoriasis-r(11), python(1) | FULLY reproducible |
| A4_study_exclusion | Study-exclusion documentation (S1–S9) | 23 | 23 | psoriasis-r(23) | FULLY reproducible |
| B_scissor | Scissor single-cell gradient | 44 | 43 | scissor-r(31), python(9), psoriasis-r(2) | FULLY reproducible |
| C_stat3_isoforms | STAT3 six-isoform characterization | 3 | 2 | python(2) | FULLY reproducible |
| D_docking_md | Molecular docking + MD | 50 | 50 | dock-md(38), python(10) | FULLY reproducible |
| E_genome_browser | Genome-browser tool | 12 | 10 | psoriasis-r(8), python(2) | FULLY reproducible |
| F_presentation | Figure presentation compile | 1 | 1 | python(1) | FULLY reproducible |
| G_audit | Project audit | 2 | 0 | — | PARTIAL (1 terminal) |
| H_session_current | Theory-1 / IL-1 / session docs | 29 | 16 | python(15) | FULLY reproducible |
| I_lncRNA | lncRNA differential expression | 11 | 11 | psoriasis-r(8), python(2) | FULLY reproducible |
| J_theory2_plan | Theory-2 sex-stratified plan | 1 | 1 | — | FULLY reproducible |
| K_refs_deconv | Reference PDFs + deconv reading | 5 | 0 | — | FULLY reproducible |

## Toolchains to capture (environment.yml / renv per arm)

- **psoriasis-r** — 127 artifacts
- **python** — 60 artifacts
- **dock-md** — 38 artifacts
- **scissor-r** — 31 artifacts

## Interpretation

This is a high-recoverability project: nearly every computed result can be regenerated from code captured in artifact lineage. The consolidation task is therefore primarily *organization* (staging each arm's code+outputs into the repo with a manifest), not *reconstruction*. The only genuinely non-regenerable items are the authored prose docs (which are themselves the deliverable) and user-uploaded reference PDFs (external, by design).