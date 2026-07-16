# psoriasis-multiomics

**Question:** Which deranged genes, pathways, and druggable proteins drive psoriatic lesion
formation — integrating bulk multi-study meta-analysis with single-cell phenotype mapping?

**This project does NOT:** make clinical, diagnostic, or causal-inference claims. It is a
**target-discovery / hypothesis-generating** program on public data.

## Canonical location
- Project root: `~/Soahum/Project/psoriasis/psoriasis-1-bulk` (CSB cluster, git repo).
- NEVER write project files outside this folder (or the paths in `docs/ACCESS.md`).

## Arms
- **Arm A — Bulk meta-analysis** (STAT3 / IL-1β): multi-study recount3 PP-vs-NN DE + random-effects meta.
- **Arm B — Scissor single-cell**: ordinal biopsy-site phenotype (NN<PN<PP) projected onto GSE173706.
- **Downstream**: deconvolution (two-stage design), druggability/docking, sequence-level (Sei-LLRA), isoforms.

## Canonical layout (project-setup skill)
| Dir / file | Contents |
|---|---|
| `docs/ACCESS.md` | host folders, inputs, outputs, data boundaries |
| `docs/COMPUTE.md` | CSB cluster + `scissor-r` env pins (reproducibility-critical) |
| `docs/DECISIONS_LOG.md` | every non-trivial choice, scoped |
| `docs/METHODS.md` | full pipeline (Arm A + Arm B) + results index |
| `notes/references.md` | paper library (DOI + specific claim used) |
| `notes/lab_notebook.md` | dated, append-only activity log |
| `notes/PLAN.md` | static mirror of the live plan |
| `cluster/STAGING_MANIFEST.md` | every file each SLURM job needs |
| `CHECKPOINT.md` | resume-from-here snapshot |
| `code/` `data/` `results/` `results_full/` `figures/` `logs/` | as before |

## Old → new file map (additive retrofit — nothing was moved or deleted)
| Existing file | Role | Canonical counterpart |
|---|---|---|
| `README.md` | Arm-B (Scissor) readme | this `PROJECT_OVERVIEW.md` is the program-level readme |
| `HANDOFF.md` | Arm-B status / reproduce / next | complemented by `CHECKPOINT.md` (+ `docs/`) |
| `checkpoint.md` (lowercase) | earlier Arm-B resume note | superseded by `CHECKPOINT.md` |
| `WHITEPAPER.md` / `.pdf` | Arm-B paper | referenced from `docs/METHODS.md` results index |
| `environment.yml` | env spec | pinned in `docs/COMPUTE.md` |

> Retrofit was **additive**: the docs/notes/cluster skeleton was added; no existing file was
> renamed, moved, or overwritten. Consolidation (e.g. merging the two checkpoint files) is left
> as an explicit future choice.
