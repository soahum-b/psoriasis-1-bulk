# Auto-extracted generating script
# Produces: fig3_scree.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): pca_diag.rds
# Source artifact version: faba4a4b-022d-4455-add8-28d6d99690d0
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({ library(ggplot2); library(dplyr) })
d <- readRDS("pca_diag.rds"); pca <- d$pca; pcv <- d$pcv

scree <- data.frame(PC=factor(paste0("PC",1:10), levels=paste0("PC",1:10)), var=pcv[1:10])
p3 <- ggplot(scree, aes(PC, var)) +
  geom_col(fill="#4C72B0", width=0.7) +
  geom_text(aes(label=sprintf("%.1f%%", var)), vjust=-0.4, size=3.2) +
  labs(title="PC1 dominates: one axis carries the disease signal",
       subtitle="Scree plot: variance explained by each principal component (eigenvalues, %)",
       x=NULL, y="Variance explained (%)") +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold"),
        panel.grid.major.x=element_blank())
ggsave("fig3_scree.png", p3, width=8, height=5, dpi=300)