# Auto-extracted generating script
# Produces: gsea_meta_primary_jointBH.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_meta_both.rds
# Source artifact version: 041bf1d6-9620-476b-8d7d-7c0a0d4f32d6
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(data.table); library(fgsea)})
gm <- readRDS("gsea_meta_both.rds")

g <- copy(gm$k4_primary)
setDT(g)

g[, padj_percollection := padj]
g[, padj_joint := p.adjust(pval, method="BH")]

fwrite(g[order(padj_joint)], "gsea_meta_primary_jointBH.csv")