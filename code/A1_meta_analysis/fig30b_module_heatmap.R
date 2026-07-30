# Auto-extracted generating script
# Produces: fig30b_module_heatmap.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds
# Source artifact version: 98f69636-4e48-4e09-bd05-559ada4d0517
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(pheatmap)

de <- readRDS("per_study_de.rds")

mod <- c("DEFB4A","DEFB4B","S100A7","S100A7A","S100A8","S100A9","S100A12","PI3","LCN2",
         "IL17A","IL17C","IL17F","IL19","IL20","IL36A","IL36G","CXCL1","CXCL8","CCL20",
         "NOS2","SERPINB3","SERPINB4","SERPINB11","SPRR2A","SPRR2B","SPRR2C","SPRR2D",
         "SPRR2E","SPRR2F","SPRR2G","LCE3A","LCE3C","LCE3D","KRT16","KRT17","KRT6A",
         "KRT6B","KRT6C","MKI67","TOP2A","CCNB1","CCNB2","CDC20","AURKA","AURKB",
         "PLK1","BUB1","BUB1B","CDKN1A","CDKN2A","MCM2","MCM4","MCM5","MCM6",
         "PCNA","TYMS","RRM2","TK1","DHFR","E2F1","E2F2","FOXM1","MYBL2","BIRC5",
         "CENPA","NDC80","KIF20A")

studies_use <- names(de$PPvsNN)

lfc_mat <- matrix(NA_real_, length(mod), length(studies_use), dimnames=list(mod, studies_use))
for (s in studies_use) {
  x <- de$PPvsNN[[s]]; x <- x[!duplicated(x$gene),]
  lfc_mat[, s] <- x$logFC[match(mod, x$gene)]
}
lfc_mat <- lfc_mat[order(-rowMeans(lfc_mat, na.rm=TRUE)),,drop=FALSE]

png("fig30b_module_heatmap.png", width=560, height=1450, res=170)
pheatmap(lfc_mat, cluster_cols=FALSE, cluster_rows=FALSE, fontsize_row=5, fontsize_col=8,
         color=colorRampPalette(c("white","#FDBB84","#C44E52","#7F0000"))(50),
         main="C0 module logFC (PP vs NN)", na_col="grey90", breaks=seq(0,8,length.out=51))
dev.off()