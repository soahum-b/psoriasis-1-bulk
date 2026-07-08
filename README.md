# Scissor-on-gradient: peri-lesional psoriasis

Projecting a bulk RNA-seq **ordinal biopsy-site phenotype** (normal `NN` <
peri-lesional `PN` < lesional `PP`) onto a single-cell reference with
**Scissor** (Gaussian mode) to find the cell states that track the
psoriasis-progression gradient.

> **Start with [`HANDOFF.md`](HANDOFF.md)** — it is the single source of truth
> for project status, results, how to reproduce/scale, and next steps.

## Quick start

```bash
conda env create -f environment.yml && conda activate scissor-r
Rscript code/00_download_data.R          # fetch GSE173706 + SRP165679
# full-census run on a cluster:
sbatch --mem=256G --cpus-per-task=16 --wrap "Rscript code/run_full_census_cluster.R"
```

## Headline result (20k-cell backbone)

- Gradient-tracking cells (**Scissor+**) are **monotonic on NN<PN<PP** and
  **peak at the peri-lesional (PN) tier** — the intermediate-state hypothesis.
- Both significance controls pass at **p = 0.000**: reliability test (real
  CV-MSE 0.147 vs null 0.779) and selection permutation null (directionality
  gap 0.638 → −0.305 under shuffled labels).
- **Endothelial cells** are the dominant gradient-trackers (5.2× enriched,
  OR 11.3); the 1,861-gene program is vascular-led. **STAT3** is significantly
  up in Scissor+ (log2FC 0.43, padj 0.018) as part of that program.

## Method

Faithful Scissor pipeline with the compiled `APML1` solver replaced by a
**pure-R glmnet network-regularized elastic net** (sparse edge augmentation of
the symmetric-normalized graph Laplacian). See `code/scissor_glmnet_solver.R`
and HANDOFF.md §3 for the derivation and the caveat about the canonical solver.

## Layout

`code/` scripts · `results/` tables (CSV tracked, RDS regenerable) ·
`figures/` 8 PNGs · `environment.yml` · `HANDOFF.md` (status & next steps).
