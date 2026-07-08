#!/usr/bin/env Rscript
# deconv_validation.R
# ---------------------------------------------------------------------------
# Orthogonal validation of the Scissor gradient result by BULK DECONVOLUTION.
# Question: do the cell types Scissor flags as gradient-tracking also show
# monotonic PROPORTION shifts across NN -> PN -> PP in the real SRP165679 bulk?
#
# Method: NNLS deconvolution against a cell-type signature built from the Ma
# 2023 single-cell reference. This is the transparent core shared by
# CIBERSORT/MuSiC-style regression deconvolution. For a benchmarked choice of
# method+reference, see the deconvBenchmarking note in HANDOFF.md (this script
# is the direct real-bulk estimate; deconvBenchmarking is the method-selection
# layer that would justify which deconvolution method to trust here).
#
# Inputs : results/reference_processed.rds (sc), data/rse_SRP165679.rds (bulk)
# Outputs: results/deconv_*.{rds,csv}, figures/fig_deconv_validation.png
# ---------------------------------------------------------------------------
suppressMessages({library(Seurat); library(Matrix); library(nnls);
                  library(SummarizedExperiment); library(ggplot2); library(reshape2)})

## 1. Build cell-type signature from the single-cell reference ----------------
so  <- readRDS("results/reference_processed.rds")
cts <- GetAssayData(so, layer = "counts")
ctype <- factor(so$celltype)
keep_ct <- names(table(ctype))[table(ctype) >= 200]      # drop tiny clusters
cellkeep <- ctype %in% keep_ct
cts <- cts[, cellkeep]; ctype <- droplevels(ctype[cellkeep])
cpm <- t(t(cts) / Matrix::colSums(cts)) * 1e6
G <- sparse.model.matrix(~0 + ctype); colnames(G) <- levels(ctype)
sig <- as.matrix((cpm %*% G) %*% Diagonal(x = 1 / Matrix::colSums(G)))
colnames(sig) <- levels(ctype)
fold <- log1p(sig) - rowMeans(log1p(sig))
mk <- unique(unlist(lapply(colnames(sig), \(ct) rownames(sig)[order(-fold[,ct])][1:150])))
S <- sig[intersect(mk, rownames(sig)), ]

## 2. Bulk CPM (93 classified samples), Ensembl -> symbol --------------------
rse  <- readRDS("data/rse_SRP165679.rds")
anch <- readRDS("results/bulk_anchor.rds")
rse2 <- rse[, match(anch$samples, colData(rse)$external_id)]
B <- t(t(assay(rse2,"counts")) / colSums(assay(rse2,"counts"))) * 1e6
gn <- rowData(rse2)$gene_name; o <- order(rowMeans(B), decreasing = TRUE)
B <- B[o,]; gn <- gn[o]; keep <- !duplicated(gn) & !is.na(gn) & gn != ""
B <- B[keep,]; rownames(B) <- gn[keep]

## 3. NNLS deconvolution per sample ------------------------------------------
common <- intersect(rownames(S), rownames(B)); S <- S[common,]; B <- B[common,]
props <- as.data.frame(t(apply(B, 2, \(b){ w <- nnls(S,b)$x; w/sum(w) })))
colnames(props) <- colnames(S)
props$tier <- anch$cls; props$tord <- c(NN=0,PN=1,PP=2)[anch$cls]

## 4. Monotonic trend test per cell type -------------------------------------
trend <- do.call(rbind, lapply(colnames(S), \(ct){
  co <- summary(lm(props[[ct]] ~ props$tord))$coefficients[2,]
  data.frame(celltype=ct, slope_per_tier=co[1], p=co[4]) }))
trend$padj <- p.adjust(trend$p, "BH"); trend <- trend[order(trend$p),]
write.csv(props, "results/deconv_proportions.csv", row.names = FALSE)
write.csv(trend, "results/deconv_trend.csv", row.names = FALSE)
saveRDS(list(props=props, trend=trend, signature=S), "results/deconv_result.rds")
print(trend)
