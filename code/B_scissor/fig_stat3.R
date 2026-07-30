# Auto-extracted generating script
# Produces: fig_stat3.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 78780cb7-0355-48bd-bfc9-0523deb1c1b7
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
de <- read.csv("results_full/gradient_program_DE_BH.csv", stringsAsFactors = FALSE)

st3 <- FetchData(so, vars="STAT3", layer="data")[,1]
dfF <- data.frame(expr=st3, scissor=so$scissor)
s <- de[de$gene=="STAT3",]
sub <- if(nrow(s)) sprintf("log2FC %.2f, raw p=%.2f, BH q=%.2f -> n.s. at full census", s$avg_log2FC, s$p_val, s$fdr_BH) else "STAT3"
pF <- ggplot(dfF, aes(scissor, expr, fill=scissor)) +
  geom_violin(scale="width", trim=TRUE) +
  geom_boxplot(width=0.12, outlier.shape=NA, fill="white") +
  scale_fill_manual(values=pal, guide="none") +
  labs(x=NULL, y="STAT3 log-normalized expression",
       title="STAT3 across Scissor classes (full census)", subtitle=sub)
ggsave(file.path(OUT,"fig_stat3.png"), pF, width=6, height=5, dpi=150)