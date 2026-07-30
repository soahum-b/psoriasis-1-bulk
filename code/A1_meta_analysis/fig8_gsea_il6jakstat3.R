# Auto-extracted generating script
# Produces: fig8_gsea_il6jakstat3.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds, gsea_ranks.rds
# Source artifact version: 7d9ab799-a300-49ab-912b-0f6d68fc1d7a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

ranks <- readRDS("gsea_ranks.rds")
de    <- readRDS("de_results_full.rds")

set.seed(42)
H <- msigdbr(species="human", collection="H")
H_list <- split(H$gene_symbol, H$gs_name)

fg_H <- fgsea(pathways=H_list, stats=ranks, minSize=10, maxSize=500, eps=0)
fg_H <- fg_H[order(fg_H$padj), ]

NNblue <- "#4C72B0"; PPred <- "#C44E52"

set.seed(42)
pw <- "HALLMARK_IL6_JAK_STAT3_SIGNALING"
gs <- H_list[[pw]]
row <- as.data.frame(fg_H[fg_H$pathway==pw, ])
nes <- row$NES; padj <- row$padj; sz <- row$size

r <- sort(ranks, decreasing=TRUE)
N <- length(r); hits <- names(r) %in% gs
Nh <- sum(hits)
Phit <- cumsum(ifelse(hits, abs(r), 0)) / sum(abs(r[hits]))
Pmiss <- cumsum(ifelse(!hits, 1, 0)) / (N - Nh)
res <- Phit - Pmiss
peak_i <- which.max(abs(res))

df_curve <- data.frame(rank=seq_len(N), RES=res)
tick <- data.frame(rank=which(hits))
le <- row$leadingEdge[[1]]

p1 <- ggplot(df_curve, aes(rank, RES)) +
  geom_hline(yintercept=0, colour="grey70", linewidth=0.3) +
  geom_line(colour=PPred, linewidth=0.9) +
  geom_vline(xintercept=peak_i, linetype="dashed", colour="grey40", linewidth=0.4) +
  geom_rug(data=tick, aes(x=rank), inherit.aes=FALSE, sides="b",
           colour=PPred, alpha=0.55, length=unit(0.06,"npc")) +
  annotate("text", x=N*0.62, y=max(res)*0.85,
           label=sprintf("NES = %.2f\nadj.P = %.1e\n%d genes", nes, padj, sz),
           hjust=0, size=4.1, colour="grey15") +
  annotate("text", x=N*0.02, y=min(res)-0.02, label="up in psoriasis (PP)", hjust=0, size=3.4, colour=PPred) +
  annotate("text", x=N*0.98, y=min(res)-0.02, label="down (NN)", hjust=1, size=3.4, colour=NNblue) +
  labs(title="GSEA enrichment - Hallmark IL6/JAK/STAT3 signaling",
       subtitle="Running enrichment score across all 24,528 ranked genes",
       x="Gene rank (high to low t-statistic)", y="Running enrichment score") +
  theme_bw(base_size=13) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"))

ggsave("fig8_gsea_il6jakstat3.png", p1, width=8.5, height=5.2, dpi=150)