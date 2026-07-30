# Auto-extracted generating script
# Produces: reference_raw.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 77f289a2-f05c-4ecc-a0b8-bf157880b0b0
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(Seurat)
library(Matrix)
library(data.table)
library(R.utils)

setDTthreads(4)

man <- read.csv("scissor_repo/data/sample_manifest.csv", stringsAsFactors=FALSE)
rawdir <- "scissor_repo/data/raw"

read_one <- function(fn){
  dt <- fread(file.path(rawdir, fn), showProgress=FALSE)
  genes <- dt[[1]]; m <- as.matrix(dt[,-1]); rownames(m) <- genes
  as(m, "CsparseMatrix")
}

mats <- vector("list", nrow(man))
for (i in seq_len(nrow(man))) {
  m <- read_one(man$file[i]); colnames(m) <- paste0(man$sample[i], "_", colnames(m)); mats[[i]] <- m
}
stopifnot(length(unique(lapply(mats, rownames)))==1)
counts <- do.call(cbind, mats)
cat("merged:", nrow(counts), "genes x", ncol(counts), "cells\n")

cell_samp <- sub("_[ACGT]+-\\d+$", "", colnames(counts))
meta <- data.frame(row.names=colnames(counts), sample=cell_samp,
                   condition=man$condition[match(cell_samp, man$sample)],
                   tier=man$tier[match(cell_samp, man$sample)],
                   donor=sub("^(NS|PN|PP)-", "", cell_samp))
meta$condition <- factor(meta$condition, levels=c("NN","PN","PP"))
print(table(meta$condition, useNA="ifany")); cat("n donors:", length(unique(meta$donor)), "\n")

so <- CreateSeuratObject(counts=counts, meta.data=meta, project="Ma2023_psoriasis")
saveRDS(so, "scissor_repo/results/reference_raw.rds")
cat("saved. dim:", nrow(so), "x", ncol(so), "\n")