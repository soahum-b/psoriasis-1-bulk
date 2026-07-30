# Auto-extracted generating script
# Produces: targets_jak_tyk2_abundance.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gradient_gene_classes.csv
# Source artifact version: 55b1d748-d371-4b77-bf65-761ad77b5284
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

sig_nodes <- data.frame(
  gene=c("TYK2","JAK1","JAK2","JAK3"),
  lfc=c(0.123, -0.320, -0.169, tg[gene=="JAK3",lfc]),
  FDR=c(0.44, 6.1e-8, 0.33, tg[gene=="JAK3",trend_FDR]),
  drug=c("deucravacitinib","JAK1i (upadacitinib etc.)","JAK2i","JAK3/pan-JAKi (tofacitinib)"))
sig_nodes$sig <- sig_nodes$FDR < 0.05

write.csv(sig_nodes, "targets_jak_tyk2_abundance.csv", row.names=FALSE)