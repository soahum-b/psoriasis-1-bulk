#!/usr/bin/env Rscript
# 00_download_data.R
# Fetch the raw inputs for the Scissor-on-gradient pipeline.
#  (1) GSE173706 scRNA-seq (Ma et al. 2023) — 33 per-sample count matrices.
#  (2) SRP165679 bulk RNA-seq (Tsoi 2019) via recount3 — the phenotype anchor.
# Run once from the repo root:  Rscript code/00_download_data.R
# Outputs land in data/. Large; re-runnable and idempotent.

suppressMessages({library(GEOquery); library(recount3); library(SummarizedExperiment)})
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

## (1) GSE173706 scRNA-seq raw counts -----------------------------------------
tarball <- "data/GSE173706_RAW.tar"
if (!file.exists(tarball)) {
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/",
                "?acc=GSE173706&format=file")
  message("Downloading GSE173706_RAW.tar (~264 MB) ...")
  download.file(url, tarball, mode = "wb")
}
if (length(list.files("data/raw", pattern = "csv.gz$")) == 0) {
  message("Extracting per-sample count matrices ...")
  untar(tarball, exdir = "data/raw")
}
message("scRNA files: ", length(list.files("data/raw", pattern = "csv.gz$")))

## (2) SRP165679 bulk RNA-seq via recount3 ------------------------------------
if (!file.exists("data/rse_SRP165679.rds")) {
  message("Fetching SRP165679 from recount3 ...")
  av   <- available_projects(organism = "human")
  proj <- subset(av, project == "SRP165679" & project_type == "data_sources")
  rse  <- create_rse(proj)
  assay(rse, "counts") <- transform_counts(rse)
  saveRDS(rse, "data/rse_SRP165679.rds")
}
message("Bulk RSE ready: data/rse_SRP165679.rds")
message("Done. Next: Rscript code/run_pipeline.R")
