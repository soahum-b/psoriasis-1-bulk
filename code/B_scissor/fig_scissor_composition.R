# Auto-extracted generating script
# Produces: fig_scissor_composition.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: f28c9639-8126-4b13-82d7-017f2d29dbc8
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)

OUT <- "figures_full"; dir.create(OUT, showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))
pal <- c("Scissor-"="#2c7fb8","Background"="#bdbdbd","Scissor+"="#d7301f")
tierpal <- c("NN"="#4575b4","PN"="#fdae61","PP"="#d73027")
msg <- function(...) cat(sprintf(...), "\n")

so   <- readRDS("results_full/reference_scissor_full.rds")
enr  <- read.csv("results_full/celltype_enrichment_BH.csv", stringsAsFactors = FALSE)

enr$celltype <- factor(enr$celltype, levels=enr$celltype[order(enr$OR)])
pC1 <- ggplot(enr, aes(celltype, OR, fill=OR>1)) + geom_col() + coord_flip() +
  geom_hline(yintercept=1, linetype=2) +
  scale_fill_manual(values=c("TRUE"="#d7301f","FALSE"="#2c7fb8"), guide="none") +
  labs(x=NULL, y="odds ratio in Scissor+", title="Cell-type enrichment (full census)",
       subtitle="Endothelial OR 5.24, top lineage")
comp <- as.data.frame(prop.table(table(so$scissor, so$condition), 1))
names(comp) <- c("scissor","tier","frac")
pC2 <- ggplot(comp, aes(scissor, frac, fill=tier)) + geom_col() +
  scale_fill_manual(values=tierpal) + labs(x=NULL, y="fraction of cells", fill="tier",
       title="Tier composition by Scissor class", subtitle="Scissor+ shifts toward PN/PP")
ggsave(file.path(OUT,"fig_scissor_composition.png"), pC1|pC2, width=11, height=4.6, dpi=150)