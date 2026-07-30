# Auto-extracted generating script
# Produces: reference_subset20k_umap.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): reference_subset20k.rds
# Source artifact version: ea99c227-abdd-41bd-a571-3fa4a28675ce
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(Seurat); library(ggplot2); library(patchwork)})
sub <- readRDS("reference_subset20k.rds")
res <- readRDS("scissor_repo/results/scissor_result.rds")

# need UMAP on the subset (was recomputed for SNN but not UMAP). Compute it.
if (!"umap" %in% names(sub@reductions)) {
  sub <- RunUMAP(sub, dims=1:20, verbose=FALSE)
}
lab <- res$scissor_label[colnames(sub)]
sub$scissor <- factor(lab, levels=c("Scissor-","Background","Scissor+"))

emb <- as.data.frame(Embeddings(sub,"umap")); colnames(emb) <- c("UMAP1","UMAP2")
emb$scissor <- sub$scissor; emb$celltype <- sub$celltype; emb$tier <- sub$condition
pal_s <- c("Scissor-"="#4C9BD4","Background"="grey85","Scissor+"="#C43C4E")

# plot background first, then selected on top
ord <- order(match(emb$scissor, c("Background","Scissor-","Scissor+")))
emb <- emb[ord,]
base <- theme_void(base_size=11)+theme(plot.title=element_text(size=11))
pA <- ggplot(emb, aes(UMAP1,UMAP2,color=scissor))+
  geom_point(size=.25, alpha=ifelse(emb$scissor=="Background",.35,.8))+
  scale_color_manual(values=pal_s, name=NULL,
                     guide=guide_legend(override.aes=list(size=2,alpha=1)))+
  labs(title="Scissor-selected cells on the ordinal gradient")+
  base+theme(legend.position=c(.1,.2))
ggsave("scissor_repo/figures/fig_scissor_umap.png", pA, width=6.5, height=5.5, dpi=200)

# composition: cell-type makeup of Scissor+ vs Scissor- vs background
comp <- as.data.frame(prop.table(table(sub$scissor, sub$celltype),1)*100)
names(comp) <- c("scissor","celltype","pct")
ct_ord <- names(sort(table(sub$celltype),decreasing=TRUE))
comp$celltype <- factor(comp$celltype, levels=rev(ct_ord))
comp$scissor <- factor(comp$scissor, levels=c("Scissor-","Background","Scissor+"))
pB <- ggplot(comp, aes(pct, celltype, fill=scissor))+
  geom_col(position="dodge")+
  scale_fill_manual(values=pal_s, name=NULL)+
  labs(title="Cell-type composition of the gradient-tracking fraction",
       x="% of selection class", y=NULL)+
  theme_classic(base_size=10)+theme(legend.position="top")
ggsave("scissor_repo/figures/fig_scissor_composition.png", pB, width=7, height=5, dpi=200)
cat("saved fig_scissor_umap.png, fig_scissor_composition.png\n")

# enrichment table: which cell types are over-represented among Scissor+/-
tab <- table(sub$scissor, sub$celltype)
cat("\ncell-type x selection counts:\n"); print(tab)
saveRDS(sub, "reference_subset20k_umap.rds")