# Auto-extracted generating script
# Produces: fig4_pc1_loadings.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): pca_diag.rds
# Source artifact version: 5b0cb529-4682-441e-8934-cc77f885f275
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({ library(ggplot2); library(dplyr) })
d <- readRDS("pca_diag.rds"); pca <- d$pca; pcv <- d$pcv

load1 <- pca$rotation[,1]
top_pos <- sort(load1, decreasing=TRUE)[1:12]
top_neg <- sort(load1, decreasing=FALSE)[1:12]
ldf <- rbind(
  data.frame(gene=names(top_pos), loading=as.numeric(top_pos), end="Up in lesional (PP)"),
  data.frame(gene=names(top_neg), loading=as.numeric(top_neg), end="Up in normal (NN)")
)
ldf$gene <- factor(ldf$gene, levels=ldf$gene[order(ldf$loading)])
p4 <- ggplot(ldf, aes(loading, gene, fill=end)) +
  geom_col(width=0.75) +
  scale_fill_manual(values=c("Up in lesional (PP)"="#C44E52","Up in normal (NN)"="#4C72B0")) +
  labs(title="PC1 eigenvector reconstructs psoriasis biology",
       subtitle="Genes with the largest PC1 loadings: the axis is built from known psoriasis markers",
       x="PC1 loading (eigenvector weight)", y=NULL, fill=NULL) +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold"), legend.position="top")
ggsave("fig4_pc1_loadings.png", p4, width=8, height=6, dpi=300)