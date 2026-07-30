# Auto-extracted generating script
# Produces: fig17_degene_heatmap.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds, dge_filt_norm.rds
# Source artifact version: 58a5459a-9bca-4099-b72b-7be108adebff
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({
  library(edgeR)
  library(pheatmap)
  library(RColorBrewer)
})

NNblue <- "#4C72B0"; PPred <- "#C44E52"

dge <- readRDS("dge_filt_norm.rds")
de_tab <- readRDS("de_results_full.rds")

logcpm <- cpm(dge, log=TRUE)

sig <- subset(de_tab, adj.P.Val < 0.05 & abs(logFC) > 1)

sig_up   <- head(sig$gene[sig$logFC > 0][order(-sig$logFC[sig$logFC>0])], 60)
sig_dn   <- head(sig$gene[sig$logFC < 0][order(sig$logFC[sig$logFC<0])], 60)
topg <- c(sig_up, sig_dn)
Z <- t(scale(t(logcpm[topg, , drop=FALSE])))
Z[Z >  3] <-  3; Z[Z < -3] <- -3

div_pal <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
brks_z  <- seq(-3, 3, length.out = 101)

ann <- data.frame(group = dge$samples$group)
rownames(ann) <- colnames(logcpm)

rowann <- data.frame(direction = ifelse(topg %in% sig_up, "up in PP", "down in PP"))
rownames(rowann) <- topg
ann_colors2 <- list(group = c(NN = NNblue, PP = PPred),
                    direction = c("up in PP" = PPred, "down in PP" = NNblue))

png("fig17_degene_heatmap.png", width=11, height=10, units="in", res=150)
pheatmap(Z, color = div_pal, breaks = brks_z,
         annotation_col = ann, annotation_row = rowann,
         annotation_colors = ann_colors2,
         show_rownames = TRUE, show_colnames = FALSE,
         fontsize_row = 5.5,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         clustering_method = "average",
         treeheight_row = 25, treeheight_col = 30,
         main = "Top 120 DE genes (60 up + 60 down in lesional skin), row z-scored log2-CPM",
         annotation_names_row = FALSE, annotation_names_col = FALSE)
dev.off()