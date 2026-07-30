# Auto-extracted generating script
# Produces: fig16_sample_corr_heatmap.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds, dge_filt_norm.rds
# Source artifact version: 56523034-6cf8-4afb-8e04-f750b4a0ee7c
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({
  library(edgeR)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

NNblue <- "#4C72B0"
PPred <- "#C44E52"

dge <- readRDS("dge_filt_norm.rds")
de_tab <- readRDS("de_results_full.rds")

logcpm <- cpm(dge, log = TRUE)

sig <- subset(de_tab, adj.P.Val < 0.05 & abs(logFC) > 1)
sig_genes <- sig$gene
X <- logcpm[sig_genes, , drop = FALSE]

ann <- data.frame(group = dge$samples$group)
rownames(ann) <- colnames(logcpm)
ann_colors <- list(group = c(NN = NNblue, PP = PPred))

cor_p <- cor(X, method = "pearson")
cor_s <- cor(X, method = "spearman")

rng <- range(c(cor_p, cor_s))
brks <- seq(floor(rng[1] * 100) / 100, 1, length.out = 101)
pal <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

png("fig16_sample_corr_heatmap.png", width = 13, height = 6.6, units = "in", res = 150)
pushViewport(viewport(layout = grid.layout(1, 2)))
draw_hm <- function(mat, ttl, col_i) {
  ph <- pheatmap(mat, color = pal, breaks = brks,
                 annotation_col = ann, annotation_row = ann,
                 annotation_colors = ann_colors,
                 show_rownames = FALSE, show_colnames = FALSE,
                 clustering_distance_rows = as.dist(1 - mat),
                 clustering_distance_cols = as.dist(1 - mat),
                 clustering_method = "average",
                 main = ttl, silent = TRUE, legend = TRUE,
                 annotation_names_row = FALSE, annotation_names_col = FALSE)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = col_i))
  grid.draw(ph$gtable)
  popViewport()
}
draw_hm(cor_p, "Pearson correlation (linear)", 1)
draw_hm(cor_s, "Spearman correlation (rank)", 2)
popViewport()
dev.off()