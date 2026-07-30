# Auto-extracted generating script
# Produces: fig_scissor_umap.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 67454b2a-d9a6-41f4-b7a3-83b8fb45e255
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)

OUT <- "figures_full"
dir.create(OUT, showWarnings = FALSE)

theme_set(theme_bw(base_size = 12))
pal <- c("Scissor-"="#2c7fb8","Background"="#bdbdbd","Scissor+"="#d7301f")

so <- readRDS("results_full/reference_scissor_full.rds")

if (!"umap" %in% Reductions(so)) {
  so <- RunUMAP(so, dims=1:20, verbose=FALSE)
}

df <- data.frame(so@reductions$umap@cell.embeddings, scissor=so$scissor)
names(df)[1:2] <- c("UMAP1","UMAP2")
df <- df[order(df$scissor),]
pB <- ggplot(df, aes(UMAP1, UMAP2, color=scissor)) +
  geom_point(size=0.15, alpha=0.6) + scale_color_manual(values=pal) +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1))) +
  labs(title="Scissor selection on the full-census reference (89,058 cells)", color=NULL) +
  theme(legend.position="bottom")
ggsave(file.path(OUT,"fig_scissor_umap.png"), pB, width=7, height=7, dpi=150)