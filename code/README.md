# code/ — scripts by arm (what actually reproduces each result)

All scripts are human-authored and derive from open public data (recount3, GSE173706).
Each arm's real scripts are listed below with their source location. Where an arm has
**no committed pipeline script**, that is stated plainly — its outputs live under
`results/_consolidated/<arm>/` but the generating code was not committed as a numbered script.

## Arms WITH real numbered/authored scripts

| Arm | Scripts | Location | How to run |
|-----|---------|----------|------------|
| Bulk meta-analysis (PP-vs-NN, STAT3/IL-1) | `step-1..step-11`, `psoriasis_pipeline.R` (master) | `code/pipeline_scripts/` | `make pipeline` (see Makefile), env `psoriasis-r` |
| STAT3 isoforms / splicing | `step-7-exon-junc-load.R`, `step-8_spliceosome.R`, `step-9_stat3_splicing.R` | `code/pipeline_scripts/` | part of the pipeline sequence |
| Pathway / TF-activity | within `psoriasis_pipeline.R` + `step-5.R` | `code/pipeline_scripts/` | pipeline; GSEA/ORA/CAMERA/decoupleR |
| NF-κB / IL-17 | `step-11_nfkb_il17.R` | `code/pipeline_scripts/` | pipeline |
| Cross-study co-expression (clust) | `run_clust.py` + `make clust` | `code/pipeline_scripts/` | env `clust-env` |
| Staging axis (NN→PN→PP) | `staging_figures.R` | `code/staging_axis/` | env `psoriasis-r` |
| Scissor single-cell gradient | `scissor_run.R`, `scissor_glmnet_solver.R`, `scissor_reliability.R`, `run_full_census_cluster.R`, `deconv_validation.R` | `code/` (top level) | cluster; see `run_full_census.sbatch` |

Setup for the above: `code/pipeline_scripts/setup.sh` + `environment.yml` / `environment-clust.yml`,
and `Makefile` targets (`make setup / pipeline / clust / all / verify`). See
`code/pipeline_scripts/README.md`, `SETUP.md`, `RESTART_HERE.md`.

## Arms WITHOUT a committed pipeline script (outputs only)

These arms were run in separate sessions; their figures/tables are archived under
`results/_consolidated/<arm>/`, but a clean numbered generating script is **not** committed here.
Reproducing them would require re-deriving the code (the analysis is described in each arm's
`results/_consolidated/<arm>/REPRODUCE.md` and in `PROJECT_AUDIT.md`).

| Arm | Outputs | Note |
|-----|---------|------|
| Molecular docking + MD | `results/_consolidated/docking_md/` | pdb/sdf/pdf/csv present; no committed pipeline script |
| lncRNA differential expression | `results/_consolidated/lncRNA/` | tables/figure present; no committed pipeline script |
| Genome-browser tool | `results/_consolidated/genome_browser/` | figures present; no committed pipeline script |

## Correction note
An earlier attempt auto-extracted per-cell code snapshots from the platform's lineage and
committed them as scripts; those were **removed** (commit `db577e7`) because (a) the real
numbered pipeline already existed and (b) language/env auto-detection was unreliable. The real
scripts above supersede them.
