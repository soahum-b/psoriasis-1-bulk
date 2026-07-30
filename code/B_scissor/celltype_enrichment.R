# Auto-extracted generating script
# Produces: celltype_enrichment.csv
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 3b9e6a06-7fc6-429c-90b4-052f37332290
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(Seurat); library(ggplot2); library(patchwork)})
sub <- readRDS("scissor_repo/results/reference_subset20k.rds")
res <- readRDS("scissor_repo/results/scissor_result.rds")

if (!"umap" %in% names(sub@reductions)) {
  sub <- RunUMAP(sub, dims=1:20, verbose=FALSE)
}
lab <- res$scissor_label[colnames(sub)]
sub$scissor <- factor(lab, levels=c("Scissor-","Background","Scissor+"))

saveRDS(sub, "scissor_repo/results/reference_subset20k_umap.rds")

sub <- readRDS("scissor_repo/results/reference_subset20k_umap.rds")
md <- sub@meta.data
overall <- prop.table(table(md$celltype))
enr <- data.frame(celltype=names(overall))
for (grp in c("Scissor+","Scissor-")){
  p <- prop.table(table(md$celltype[md$scissor==grp]))
  enr[[grp]] <- as.numeric(p[enr$celltype]/overall[enr$celltype])
}
enr$pos_OR <- NA; enr$pos_p <- NA
for (i in seq_len(nrow(enr))){
  ct <- enr$celltype[i]
  tt <- table(md$celltype==ct, md$scissor=="Scissor+")
  if (all(dim(tt)==2)){ ft <- fisher.test(tt); enr$pos_OR[i] <- ft$estimate; enr$pos_p[i] <- ft$p.value }
}
enr <- enr[order(-enr$`Scissor+`),]
cat("Cell-type enrichment (fold vs overall abundance):\n")
print(format(enr, digits=2))
write.csv(enr, "celltype_enrichment.csv", row.names=FALSE)
cat("\nsaved celltype_enrichment.csv\n")