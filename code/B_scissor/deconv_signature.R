# Auto-extracted generating script
# Produces: deconv_signature.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): reference_processed.rds
# Source artifact version: fd8dc30a-5768-4849-ae7b-9777dbfac5a9
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(Seurat); library(Matrix)})
so <- readRDS("reference_processed.rds")
# signature = mean CPM-normalized expression per cell type, over marker-informative genes.
# Use counts -> CPM per cell -> mean per celltype (standard signature construction).
cts <- GetAssayData(so, layer="counts")
ctype <- factor(so$celltype)
# drop tiny Tcell cluster (62 cells) - too few for a stable signature
keep_ct <- names(table(ctype))[table(ctype) >= 200]
cat("cell types kept (>=200 cells):", paste(keep_ct, collapse=", "), "\n")
cellkeep <- ctype %in% keep_ct
cts <- cts[, cellkeep]; ctype <- droplevels(ctype[cellkeep])

# CPM normalize each cell
libsize <- Matrix::colSums(cts)
cpm <- t(t(cts)/libsize)*1e6
# mean CPM per cell type (sparse-friendly via aggregation matrix)
G <- sparse.model.matrix(~0+ctype); colnames(G) <- levels(ctype)
n_per <- Matrix::colSums(G)
sig <- as.matrix((cpm %*% G) %*% Diagonal(x=1/n_per))
colnames(sig) <- levels(ctype)
cat("signature dim:", paste(dim(sig),collapse=" x "), "\n")

# select informative genes: top marker genes per cell type by fold vs rowmean
logsig <- log1p(sig)
fold <- logsig - rowMeans(logsig)
marker_genes <- unique(unlist(lapply(colnames(sig), function(ct){
  rownames(sig)[order(-fold[,ct])][1:150]
})))
marker_genes <- intersect(marker_genes, rownames(sig))
cat("marker genes for signature:", length(marker_genes), "\n")
sig_m <- sig[marker_genes,]
saveRDS(list(signature=sig_m, celltypes=colnames(sig_m), all_sig=sig),
        "deconv_signature.rds")
cat("saved deconv_signature.rds\n")