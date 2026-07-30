# Environment setup

The pipeline needs two environments: an **R stack** for the main analysis and a
small **Python env** for the `clust` co-expression tool (STEP 11.8), which
requires `numpy<2` and is therefore kept separate.

## Quick start (conda, recommended)

```bash
bash setup.sh
```

This creates `psoriasis-r` and `clust-env`, then verifies both. It uses `mamba`
if available (much faster for the Bioconductor stack), otherwise `conda`.

To run the analysis:

```bash
conda activate psoriasis-r
Rscript psoriasis_pipeline.R
```

## Manual setup

Create the environments individually:

```bash
conda env create -f environment.yml         # -> psoriasis-r
conda env create -f environment-clust.yml   # -> clust-env
```

## Plain R (no conda)

If you already have R (>= 4.4) and prefer not to use conda:

```bash
Rscript setup.R
```

`setup.R` installs any missing CRAN and Bioconductor packages via
`install.packages()` and `BiocManager::install()`, then verifies they all load.
You still need `clust` separately: `pip install clust` in a `numpy<2` env.

## What gets installed

**R (Bioconductor):** recount3, SummarizedExperiment, edgeR, limma, fgsea,
clusterProfiler, decoupleR, qvalue, org.Hs.eg.db
**R (CRAN):** msigdbr, metafor, UpSetR, data.table, dplyr, ggplot2, ggrepel,
patchwork, pheatmap, RColorBrewer, gridExtra, BiocManager
**Python (clust-env):** clust, numpy<2, scipy, scikit-learn, pandas, portalocker, joblib

## Notes

- GSVA is **not** required: per-sample pathway scores are computed with a
  transparent base-R ssGSEA implementation inside the pipeline (STEP 10), so
  there is no dependency on the (build-heavy) Bioconductor `GSVA` package.
- `metafor` and `UpSetR` are pure-R CRAN packages loaded lazily where used.
- The pipeline was developed on R 4.5.3 with Bioconductor 3.22; any R >= 4.4
  with the matching Bioconductor release should reproduce it.
