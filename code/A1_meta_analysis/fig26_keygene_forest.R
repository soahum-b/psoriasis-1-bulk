# Auto-extracted generating script
# Produces: fig26_keygene_forest.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): meta_de_results.rds
# Source artifact version: 4b5b258a-aaea-4b6e-894e-6315bc5ee9ba
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(metafor)
library(data.table)

mp <- as.data.table(readRDS("meta_de_results.rds")$PPvsNN)

kg <- c("STAT3","STAT1","SOCS3","IL17A","DEFB4A","S100A7","S100A9","IL36G","CXCL8","CCR1","C1QB","EFTUD2")
kgd <- mp[gene %in% kg][match(kg, gene)]
kgd <- kgd[!is.na(gene)]
kgd[, ci.lb := logFC-1.96*SE][, ci.ub := logFC+1.96*SE]

png("fig26_keygene_forest.png", width=1450, height=1050, res=190)
par(mar=c(5,4,3,2))
ord <- order(kgd$logFC)
kk <- kgd[ord]
forest(x=kk$logFC, sei=kk$SE, slab=kk$gene,
       xlab="pooled log2 fold-change (random-effects, PP vs NN)",
       header=c("Gene","logFC [95% CI]"), col="#C44E52",
       annotate=TRUE, cex=0.9, psize=1.2, pch=18,
       ilab=sprintf("k=%d, I2=%.0f%%", kk$k, kk$I2), ilab.xpos=-8)
title("Figure 26. Pooled effects for key psoriasis genes", cex.main=1.0)
dev.off()