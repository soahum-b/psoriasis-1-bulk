# Auto-extracted generating script
# Produces: fig25_stat3_forest.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds
# Source artifact version: 45be30ba-2771-4578-a18f-d666cbb8f5a7
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

wslib <- "/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/.r-libs/psoriasis-r"
.libPaths(c(wslib, .libPaths()))
suppressMessages({library(metafor); library(data.table)})
de <- readRDS("per_study_de.rds")

get_gene_studies <- function(contrast, g) {
  rbindlist(lapply(names(de[[contrast]]), function(s){
    x <- as.data.table(de[[contrast]][[s]]); r <- x[gene==g]
    if(nrow(r)) data.table(study=s, yi=r$logFC, sei=r$SE, n1=r$n1, n2=r$n2) else NULL
  }))
}
s3 <- get_gene_studies("PPvsNN","STAT3")
res <- rma(yi=yi, sei=sei, data=s3, method="DL")

png("fig25_stat3_forest.png", width=1500, height=950, res=190)
par(mar=c(5,4,3,2))
slab <- sprintf("%s  (%d PP vs %d NN)", s3$study, s3$n1, s3$n2)
forest(res, slab=slab, xlab="log2 fold-change (lesional vs healthy)",
       header=c("Study","logFC [95% CI]"), col="#C44E52", border="#C44E52",
       cex=0.95, mlab="RE model (DerSimonian-Laird)", psize=1.3, xlim=c(-4.5,4.2))
text(-4.5, -1.6, pos=4, cex=0.8, col="#333333",
     bquote(paste("Heterogeneity:  ", I^2, " = ", .(sprintf("%.0f",res$I2)), "%,   ",
       tau^2, " = ", .(sprintf("%.3f",res$tau2)), ",   Q(", .(res$k-1), ") = ",
       .(sprintf("%.1f",res$QE)), ",  p = ", .(sprintf("%.1e",res$QEp)))))
title("Figure 25. STAT3 random-effects forest (PP vs NN)", cex.main=1.0)
dev.off()