# Auto-extracted generating script
# Produces: targets_gradient_position.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gradient_gene_classes.csv
# Source artifact version: 910debc9-5d2f-4d98-8059-981a176a46f0
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)

gc <- read.csv("gradient_gene_classes.csv")
setDT(gc)

targets <- c("STAT3","JAK1","JAK2","JAK3","TYK2","IL23A","IL12B","IL17A","IL17F","IL36G","IL36A",
             "TNF","IL6","IFNG","S100A7","S100A8","S100A9","DEFB4A","PI3","CXCL8","CCL20","SOCS3")
tg <- gc[gene %in% targets]
tg[, npn := (mu_PN-mu_NN)/(mu_PP-mu_NN)]
tg[, lfc := mu_PP - mu_NN]

act <- tg[, .(gene, cls, lfc, npn, FDR=trend_FDR)]
act <- act[lfc > 0.5 & FDR < 0.05]
act[, grp := ifelse(npn < 0.15, "late_PP (plaque-specific)", "progressive (rises through PN)")]
setorder(act, npn)

write.csv(act, "targets_gradient_position.csv", row.names=FALSE)