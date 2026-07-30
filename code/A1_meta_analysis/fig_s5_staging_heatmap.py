# Auto-extracted generating script
# Produces: fig_s5_staging_heatmap.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: aa4cf1a4-3349-43f6-9c23-2be752b9b5ce
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(pheatmap)
library(grid)

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")

mat <- fread("clust_input/SRP165679.tsv")
X <- as.matrix(mat[,-1])
rownames(X) <- mat$Genes

cls <- as.data.table(readRDS("sample_classification.rds"))
lab <- cls[srp=="SRP165679"][class %in% c("NN","PN","PP")]
X <- X[, lab$external_id]
grp <- factor(lab$class, levels=c("NN","PN","PP"))

g <- readRDS("staging_axis/results/gradient_de_genes.rds")
g[, monotonic_flag := ifelse(mono_up,"up", ifelse(mono_down,"down","non"))]
up  <- g[monotonic_flag=="up"][order(-trend_slope)][1:40, gene]
dn  <- g[monotonic_flag=="down"][order(trend_slope)][1:40, gene]
sel <- c(up, dn)
M <- X[sel, ]
Z <- t(scale(t(M)))
ord <- order(grp, colMeans(Z[1:40,,drop=FALSE]))
Z <- Z[, ord]
grp_ord <- grp[ord]
Z[Z>2.5] <- 2.5
Z[Z< -2.5] <- -2.5

ann_col <- data.frame(Stage=grp_ord)
rownames(ann_col) <- colnames(Z)
ann_row <- data.frame(Direction=rep(c("up NN<PN<PP","down NN>PN>PP"),c(40,40)))
rownames(ann_row) <- rownames(Z)
anncols <- list(Stage=c(NN="#4C72B0",PN="#DD8452",PP="#C44E52"),
                Direction=c("up NN<PN<PP"="#C44E52","down NN>PN>PP"="#4C72B0"))

png("fig_s5_staging_heatmap.png", width=1500, height=1150, res=150)
pheatmap(Z, cluster_cols=FALSE, cluster_rows=FALSE,
  annotation_col=ann_col, annotation_row=ann_row, annotation_colors=anncols,
  show_colnames=FALSE, fontsize_row=5.5, gaps_row=40,
  color=colorRampPalette(c("#3B4CC0","white","#B40426"))(101),
  main="Molecular staging axis: top monotonic genes ordered NN -> PN -> PP (SRP165679)")
dev.off()