# Auto-extracted generating script
# Produces: gradient_program_DE.csv
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 3ae9be22-9daa-47b7-85ea-59ae05a6a0bb
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(Seurat)})
sub <- readRDS("scissor_repo/results/reference_subset20k_umap.rds")
Idents(sub) <- sub$scissor

de <- FindMarkers(sub, ident.1="Scissor+", ident.2="Background",
                  logfc.threshold=0.1, min.pct=0.1, test.use="wilcox")
de$gene <- rownames(de)
de <- de[order(-de$avg_log2FC),]
write.csv(de, "scissor_repo/results/gradient_program_DE.csv", row.names=FALSE)