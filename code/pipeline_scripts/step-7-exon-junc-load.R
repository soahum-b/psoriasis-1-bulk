# 7. Load Exon- and Junction-Level Data ----
# Purpose: extend the gene-level RSE from Step 1 with exon and junction
# resolution objects from the same recount3 project. This is the foundation
# for all downstream alternative splicing analysis (Steps 8-10).

# 1. Setup & Load Packages -----
library(recount3)
library(SummarizedExperiment)

# `human_projects` and `proj_info` are still in your environment from Step 1.
# If you've restarted R, re-run Step 1 first.

# 2. Build Exon-Level RSE -----
rse_exon <- create_rse(proj_info, type = "exon")
assay(rse_exon, "counts") <- transform_counts(rse_exon)

# 3. Build Junction-Level RSE -----
# Note: junction RSEs from recount3 are already in raw count format,
# so transform_counts() is NOT needed here.
rse_jxn <- create_rse(proj_info, type = "jxn")

# Keep a friendly handle on the count matrix
exon_counts <- assay(rse_exon, "counts")
jxn_counts  <- assay(rse_jxn,  "counts")

# 4. Sanity Checks -----
cat("=== Dimension check (rows = features, cols = samples) ===\n")
cat(sprintf("Genes:     %s\n", paste(dim(rse_gene), collapse = " x ")))
cat(sprintf("Exons:     %s\n", paste(dim(rse_exon), collapse = " x ")))
cat(sprintf("Junctions: %s\n", paste(dim(rse_jxn),  collapse = " x ")))

cat("\n=== Available assays ===\n")
cat("Gene:     ", paste(assayNames(rse_gene), collapse = ", "), "\n")
cat("Exon:     ", paste(assayNames(rse_exon), collapse = ", "), "\n")
cat("Junction: ", paste(assayNames(rse_jxn),  collapse = ", "), "\n")

# 5. Verify Sample Alignment -----
# All three RSEs MUST have the same sample order, otherwise downstream
# contrasts will silently produce garbage.
stopifnot(identical(colnames(rse_gene), colnames(rse_exon)))
stopifnot(identical(colnames(rse_gene), colnames(rse_jxn)))

cat("\nSample alignment across gene / exon / junction RSEs: OK\n")

# 6. Inspect Metadata Structure -----
# Useful when you later need to map exons/junctions back to genes
cat("\n=== rowData columns ===\n")
cat("Exon: ", paste(colnames(rowData(rse_exon)), collapse = ", "), "\n")
cat("Jxn:  ", paste(colnames(rowData(rse_jxn)),  collapse = ", "), "\n")

# Save objects so downstream steps don't re-download
save(rse_exon, rse_jxn, file = "rse_exon_jxn.RData")
cat("\nSaved: rse_exon_jxn.RData\n")
