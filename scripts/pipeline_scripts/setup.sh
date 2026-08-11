#!/usr/bin/env bash
# One-command environment setup for the psoriasis meta-analysis.
#
# Usage:   bash setup.sh
#
# Creates two conda environments:
#   - psoriasis-r : the R analysis stack (recount3, edgeR, limma, fgsea, ...)
#   - clust-env   : Python env for the `clust` co-expression tool (numpy<2)
#
# Requires conda (or mamba) on PATH. If you have mamba, it is used automatically
# because it resolves the Bioconductor stack much faster.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Prefer mamba if available (faster solver)
if command -v mamba >/dev/null 2>&1; then
  CONDA=mamba
elif command -v conda >/dev/null 2>&1; then
  CONDA=conda
else
  echo "ERROR: neither mamba nor conda found on PATH." >&2
  echo "Install Miniforge (https://github.com/conda-forge/miniforge) and re-run." >&2
  exit 1
fi
echo "Using solver: $CONDA"

echo "==> Creating R environment (psoriasis-r) ..."
$CONDA env create -f "$HERE/environment.yml" || \
  $CONDA env update -f "$HERE/environment.yml"

echo "==> Creating clust environment (clust-env) ..."
$CONDA env create -f "$HERE/environment-clust.yml" || \
  $CONDA env update -f "$HERE/environment-clust.yml"

echo "==> Verifying R packages ..."
conda run -n psoriasis-r Rscript -e '
  pkgs <- c("recount3","edgeR","limma","fgsea","clusterProfiler","decoupleR",
            "qvalue","org.Hs.eg.db","SummarizedExperiment","msigdbr","metafor",
            "UpSetR","data.table","dplyr","ggplot2","ggrepel","patchwork",
            "pheatmap","RColorBrewer")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) { cat("MISSING:", paste(miss, collapse=", "), "\n")
    cat("Run: conda run -n psoriasis-r Rscript setup.R  to install the rest.\n") }
  else cat("All", length(pkgs), "R packages OK.\n")
'

echo "==> Verifying clust ..."
conda run -n clust-env python -c "import clust; print('clust import OK')" || \
  echo "clust check failed - see environment-clust.yml"

cat <<'MSG'

Setup complete.

To run the pipeline:
  conda activate psoriasis-r
  Rscript psoriasis_pipeline.R

The clust step (STEP 11.8) runs from the clust-env:
  conda activate clust-env
  python -m clust clust_input -o clust_out -t 0
MSG
