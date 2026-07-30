#!/usr/bin/env Rscript
# Fallback installer: installs any missing packages the pipeline needs.
# Use this if you prefer plain R over conda, or to top up a partial install.
#   Rscript setup.R
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c("msigdbr","metafor","UpSetR","data.table","dplyr",
               "ggplot2","ggrepel","patchwork","pheatmap","RColorBrewer",
               "gridExtra","BiocManager")
bioc_pkgs <- c("recount3","SummarizedExperiment","edgeR","limma","fgsea",
               "clusterProfiler","decoupleR","qvalue","org.Hs.eg.db")

install_if_missing <- function(pkgs, installer) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing: ", p)
      installer(p)
    } else message("Already present: ", p)
  }
}

install_if_missing(cran_pkgs, function(p) install.packages(p))
# BiocManager pins to the Bioc release matching your R version
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install_if_missing(bioc_pkgs, function(p) BiocManager::install(p, update = FALSE, ask = FALSE))

# Final verification
all_pkgs <- c(cran_pkgs, bioc_pkgs)
missing  <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Still missing after install: ", paste(missing, collapse = ", "))
} else {
  message("\nAll ", length(all_pkgs), " R packages are installed and loadable.")
}
