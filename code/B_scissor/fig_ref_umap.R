# Auto-extracted generating script
# Produces: fig_ref_umap.png
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: bb001486-de6c-4552-8d17-d5eca80418c9
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(Seurat); library(ggplot2); library(patchwork)})
so2 <- readRDS("scissor_repo/results/reference_processed.rds")

ct_levels <- names(sort(table(so2$celltype), decreasing=TRUE))
so2$celltype <- factor(so2$celltype, levels=ct_levels)
ct_pal <- setNames(c("#4C72B0","#DD8452","#55A868","#C44E52","#8172B3",
                     "#937860","#DA8BC3","#8C8C8C","#CCB974")[seq_along(ct_levels)], ct_levels)
pal_tier <- c(NN="#4C9BD4", PN="#E8A33D", PP="#C43C4E")

emb <- as.data.frame(Embeddings(so2,"umap")); colnames(emb) <- c("UMAP1","UMAP2")
emb$celltype <- so2$celltype; emb$tier <- so2$condition
set.seed(1); emb <- emb[sample(nrow(emb)),]  # avoid overplot bias

# label positions = cluster centroids
cent <- aggregate(cbind(UMAP1,UMAP2)~celltype, emb, median)

base <- theme_void(base_size=11) + theme(plot.title=element_text(size=11,hjust=0))
pA <- ggplot(emb, aes(UMAP1,UMAP2,color=celltype)) +
  geom_point(size=.15, alpha=.5) +
  geom_text(data=cent, aes(label=celltype), color="black", size=3, fontface="bold") +
  scale_color_manual(values=ct_pal, guide="none") +
  annotate("segment", x=min(emb$UMAP1), xend=min(emb$UMAP1)+3, y=min(emb$UMAP2), yend=min(emb$UMAP2),
           arrow=arrow(length=unit(.2,"cm")), linewidth=.4) +
  labs(title="Cell types (broad lineage)") + base

pB <- ggplot(emb, aes(UMAP1,UMAP2,color=tier)) +
  geom_point(size=.15, alpha=.5) +
  scale_color_manual(values=pal_tier, name="Tier",
                     guide=guide_legend(override.aes=list(size=2,alpha=1))) +
  labs(title="Ordinal tier (NN < PN < PP)") + base +
  theme(legend.position=c(.08,.2))

fig <- pA + pB + plot_annotation(tag_levels="a")
ggsave("scissor_repo/figures/fig_ref_umap.png", fig, width=11, height=5, dpi=200)
cat("saved fig_ref_umap.png\n")