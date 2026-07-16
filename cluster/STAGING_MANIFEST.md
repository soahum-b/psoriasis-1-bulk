# Staging Manifest

The DEFINITIVE list of every file a remote job needs. Check each off before submitting.

## Job: Scissor full-census (SLURM 56882314 — COMPLETED 2026-07-14)
| File | Staged? | Notes |
|---|---|---|
| `code/run_full_census_cluster.R` | [x] | driver: build ref, tune alpha, solve, significance |
| `code/scissor_glmnet_solver.R` | [x] | pure-R network elastic-net solver |
| `code/scissor_run.R` | [x] | Scissor orchestration |
| `code/scissor_reliability.R` | [x] | reliability + permutation-null controls |
| `code/run_full_census.sbatch` | [x] | SLURM script; `--mem`/`--cpus`; mail directives |
| `data/GSE173706_RAW.tar` | [x] | single-cell reference (252 MB) |
| `data/rse_SRP165679.rds` | [x] | bulk phenotype anchor |
| `data/sample_classification.rds` | [x] | tier labels |
| env `scissor-r` | [x] | `conda env create -f environment.yml` |
| **Outputs** → `results_full/` | [x] | reference_scissor_full.rds, gradient_program_DE.csv, reliability_test.rds, permutation_null.rds |

## Job: Compiled-Scissor (APML1) cross-check — NOT YET STAGED
| File | Staged? | Notes |
|---|---|---|
| Scissor source tarball (APML1) | [ ] | needs build on node; was blocker locally |
| `data/rse_SRP165679.rds` | [x] | reuse |
| reference Seurat object | [x] | reuse `results_full/reference_scissor_full.rds` |
| cross-check script | [ ] | compare β / selection vs glmnet port |
