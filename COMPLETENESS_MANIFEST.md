# Project-folder completeness manifest

*Audit date: 2026-08-03. Verifies that every project file — scripts, figures, data
(downloaded or referenced), documents — exists in the project folder on disk or in
git, not only in the ephemeral artifact store.*

## Result: no content gap

326 project artifacts were audited against the repo (git-tracked) + cluster disk
(`/net/dali/home/mscbio/sba50/Soahum/Project/psoriasis/psoriasis-1-bulk` on the CSB
node). Every human-readable file and every data object now lives in the project
folder. The only store-only items are session scaffolding (transfer tarballs, a
plan JSON, a screenshot) that are not project content.

## Where each category lives

| Category | Location | In git? | Notes |
|---|---|---|---|
| Scripts (`.R`/`.py`/`.sh`) | `code/`, `code/pipeline_scripts/`, `code/staging_axis/` | ✅ tracked | Real numbered pipeline (step-1..11, psoriasis_pipeline.R, Makefile) + Scissor scripts |
| Figures (`.png`) | `figures/`, `results/_consolidated/*/` | ✅ tracked | Incl. significance-annotated fig_s4, fig30a, fig14 |
| Whitepapers | `WHITEPAPER.md`, `docs/whitepaper_pdfs/*.pdf`, `results/_consolidated/*/psoriasis_*_whitepaper.md` | ✅ tracked | MD sources + figure-embedded PDFs |
| Presentation | `results/_consolidated/presentation/*_UPDATED.pptx` | ✅ tracked | With per-figure script refs |
| Small result tables (`.csv`) | `results/`, `results/_consolidated/` | ✅ tracked | `!results/*.csv` in .gitignore |
| Docs (methods, decisions, plan, notebook) | `docs/`, `notes/` | ✅ tracked | |
| **Downloaded data** — GSE173706 (10x) | `data/GSE173706_RAW.tar`, `data/raw/` | ⬛ on disk, gitignored | Re-fetch: `code/00_download_data.R` |
| **Downloaded data** — recount3 SRP objects | `data/rse_SRP165679.rds` + recount3 cache | ⬛ on disk, gitignored | Re-fetch via recount3 in pipeline |
| **Downloaded dependency** — GenomeInfoDbData | `data/GenomeInfoDbData_1.2.13.tar.gz` | ⬛ on disk, gitignored | Bioconductor annotation package |
| Reference objects (Scissor) | `results/reference_*.rds` (2 GB) | ⬛ on disk, gitignored | reference_processed 1.2GB, reference_raw 512MB, subset20k 276MB |
| Regenerable checkpoints (`.rds`) | `results/*.rds` (40 files) | ⬛ on disk, gitignored | voom_fit, per_study_de, dge_filt_norm, gsea_*, deconv_*, clust_*, etc. |
| Full-census Scissor outputs | `results_full/` (1.3 GB) | ⬛ on disk, gitignored | reference_scissor_full.rds, DE tables, tuning/reliability |

## Policy (why bulk data is on disk but not in git)

The `.gitignore` deliberately keeps large regenerable objects (`.rds`, `.tar`,
`results_full/`, downloaded raw data) **on disk in the project folder** but **out of
git history**, to avoid a multi-GB repo. This satisfies "nothing ephemeral, all in
the project folder" — the files exist on the cluster filesystem — while keeping the
git repo clonable. Every such file is either re-fetchable (downloaded data, via
`code/00_download_data.R` / recount3) or regenerable (pipeline checkpoints, via the
numbered scripts). The artifact store holds backup copies of all of them.

## Verification performed

- 269/326 artifacts confirmed already in repo/disk by basename.
- 37 regenerable checkpoints (123 MB) that had been swept from disk during an idle
  cleanup were restored from the artifact store to `results/` and `data/`.
- 5 whitepaper/README markdown "gaps" were confirmed byte-identical (md5) to
  canonical repo copies under different names — no duplication committed.
- 2 rendered PDF deliverables (PROJECT_AUDIT.pdf, docking_whitepapers_combined.pdf)
  added to `docs/`.
- Big reference `.rds` (2 GB) confirmed present at `results/reference_*.rds` by
  `find` + byte size.
