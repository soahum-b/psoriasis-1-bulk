# Auto-extracted generating script
# Produces: reference_processed.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): reference_pca.rds
# Source artifact version: 701f5d00-6716-4e99-a423-4e21a3c7b061
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(Seurat)})

so2 <- readRDS("reference_pca.rds")
# canonical broad-lineage markers for skin
markers <- list(
  Keratinocyte = c("KRT14","KRT5","KRT10","KRT1","KRT15","KRTDAP","SBSN"),
  Fibroblast   = c("COL1A1","COL1A2","COL3A1","PDGFRA","LUM","DCN"),
  Myeloid      = c("LYZ","CD68","CD14","AIF1","ITGAX","FCGR3A"),
  DC           = c("CD1C","CLEC9A","LAMP3","CD207"),
  Tcell        = c("CD3D","CD3E","CD2","TRAC","IL7R"),
  NK           = c("NKG7","GNLY","KLRD1"),
  Bcell        = c("MS4A1","CD79A","IGHG1","MZB1"),
  Endothelial  = c("PECAM1","VWF","CLDN5","CCL21"),
  Melanocyte   = c("MLANA","PMEL","TYRP1","DCT"),
  Mast         = c("TPSAB1","TPSB2","CPA3"),
  Muscle       = c("ACTA2","TAGLN","MYH11","DES")
)
markers <- lapply(markers, function(g) intersect(g, rownames(so2)))
# mean scaled expression per cluster per lineage
avg <- AverageExpression(so2, features=unique(unlist(markers)), group.by="seurat_clusters", layer="data")$RNA
# z-score genes across clusters
zz <- t(scale(t(as.matrix(avg))))
score <- sapply(markers, function(g) colMeans(zz[intersect(g,rownames(zz)),,drop=FALSE], na.rm=TRUE))
assign <- colnames(score)[apply(score,1,which.max)]
names(assign) <- rownames(score)  # cluster ids

# strip the 'g' prefix AverageExpression added
clust_ids <- sub("^g","",names(assign))
names(assign) <- clust_ids
so2$celltype <- factor(unname(assign[as.character(so2$seurat_clusters)]))
cat("any NA in celltype:", sum(is.na(so2$celltype)), "\n")
cat("cell-type counts:\n"); print(sort(table(so2$celltype), decreasing=TRUE))
cat("\ncelltype x tier:\n"); print(table(so2$celltype, so2$condition))
saveRDS(so2, "reference_processed.rds")
cat("\nsaved reference_processed.rds\n")