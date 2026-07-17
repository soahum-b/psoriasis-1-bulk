# Access & Data Boundaries

## In-scope host folders
| Path | Mode (ro/rw) | Purpose |
|---|---|---|
| `~/Soahum/Project/psoriasis/psoriasis-1-bulk` (n003 / CSB) | rw | **canonical project root** (git repo) |
| `~/Soahum/Recount3_project/Multi_Study_Psoriasis` (n003) | rw | per-study bulk recount3 pipeline outputs |
| Claude Science artifact store (`proj_e636854c6fa3`) | rw | all session-produced figures/tables/reports |

> Files must only be written under the paths listed here. The other granted local
> folders (`ethan-soahum-project`, `Sex_Melanoma`, `Recount-3`/Kostka) are UNRELATED
> projects and are OFF-LIMITS for this work.

## Inputs
| Source | ID / path | Notes |
|---|---|---|
| Single-cell reference | **GSE173706** (Ma et al. 2023) | 10x droplet; 96,088 cells → 89,058 QC; NN/PN/PP tiers; `data/GSE173706_RAW.tar` |
| Bulk phenotype anchor (Scissor arm) | **SRP165679** (Tsoi et al. 2019, recount3) | only tier-balanced study; `data/rse_SRP165679.rds` |
| Bulk meta-analysis studies (PP-vs-NN) | **SRP035988, SRP165679, SRP126422, SRP065812** | recount3; all treatment-free; k=4 set |
| Bulk candidate pool | 23 human recount3 psoriasis studies | screened in `study_eligibility_PPvsNN.csv` |
| Blood arm candidate | **ERP110814** (E-MTAB-6555, Tsoi/QMUL) | blood, PLAQUE psoriasis; wk-0 baseline (10 tx-naive) recoverable via ArrayExpress SDRF; no healthy ctrl |
| Blood arm — EXCLUDED | SRP173379 / SRP173378 / SRP132160 | all Generalized Pustular Psoriasis (GPP) — different subtype; dropped per user |
| Blood — needed | plaque-psoriasis PBMC single-cell + healthy-blood control | to be screened in GEO/CELLxGENE (recount3 is bulk-only) |

## Outputs
- Cluster repo: `results/` (small CSVs tracked), `results_full/` (large RDS, gitignored), `figures/`.
- Claude Science artifact store: canonical home for all deliverable figures/tables/reports.

## Sensitivity / sharing constraints
- All inputs are public (GEO / recount3-SRA); no controlled-access or patient-identifiable data.
- Large binaries (`*.tar`, `*.rds`, `results_full/`) are gitignored — regenerable; never committed.
- Clinical interpretation is preliminary; outputs are not for patient care.
