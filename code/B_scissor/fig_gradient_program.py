# Auto-extracted generating script
# Produces: fig_gradient_program.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 965a129d-a536-4b81-a9e7-461d9a00d7a1
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(ggrepel)

OUT <- "figures_full"
dir.create(OUT, showWarnings = FALSE)

theme_set(theme_bw(base_size = 12))

so <- readRDS("results_full/reference_scissor_full.rds")
de <- read.csv("results_full/gradient_program_DE_BH.csv", stringsAsFactors = FALSE)

de$sig <- de$fdr_BH < 0.05
de$neglog <- -log10(pmax(de$fdr_BH, 1e-300))
lab <- de[de$gene %in% c("CCL14","ACKR1","RAMP3","PLVAP","APLNR","CYTL1","SPNS2","PECAM1","VWF","STAT3"),]
pE <- ggplot(de, aes(avg_log2FC, neglog, color=sig)) +
  geom_point(size=0.5, alpha=0.5) +
  scale_color_manual(values=c("TRUE"="#d7301f","FALSE"="grey70"), guide="none") +
  ggrepel::geom_text_repel(data=lab, aes(label=gene), color="black", size=3, max.overlaps=20) +
  labs(x="avg log2FC (Scissor+ vs background)", y="-log10 BH-FDR",
       title="Gradient-tracking gene program (full census, BH-FDR)",
       subtitle=sprintf("%d/%d genes at BH q<0.05; vascular-led", sum(de$sig), nrow(de)))
ggsave(file.path(OUT,"fig_gradient_program.png"), pE, width=7, height=5.5, dpi=150)