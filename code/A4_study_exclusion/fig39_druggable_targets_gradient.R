# Auto-extracted generating script
# Produces: fig39_druggable_targets_gradient.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): clust_module_lfc.rds, gradient_gene_classes.csv, clust_module_genes.csv
# Source artifact version: 5ef2b532-b1b9-4465-b3e9-2ae35c0548ea
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(data.table)

theme_set(theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=12), plot.title.position="plot"))
BLUE<-"#4C72B0"; RED<-"#C44E52"; GREY<-"#8C8C8C"; GOLD<-"#DD8452"

gc <- read.csv("gradient_gene_classes.csv")
setDT(gc)

cm <- read.csv("clust_module_genes.csv")
mod_genes <- cm$gene
mg <- gc[gene %in% mod_genes]

targets <- c("STAT3","JAK1","JAK2","JAK3","TYK2","IL23A","IL12B","IL17A","IL17F","IL36G","IL36A",
             "TNF","IL6","IFNG","S100A7","S100A8","S100A9","DEFB4A","PI3","CXCL8","CCL20","SOCS3")
tg <- gc[gene %in% targets]
tg[, npn := (mu_PN-mu_NN)/(mu_PP-mu_NN)]
tg[, lfc := mu_PP - mu_NN]

act <- tg[, .(gene, cls, lfc, npn, FDR=trend_FDR)]
act <- act[lfc > 0.5 & FDR < 0.05]
act[, grp := ifelse(npn < 0.15, "late_PP (plaque-specific)", "progressive (rises through PN)")]
setorder(act, npn)
act[, gene_f := factor(gene, levels=rev(gene))]
act[, is_stat3 := gene=="STAT3"]

sig_nodes <- data.frame(
  gene=c("TYK2","JAK1","JAK2","JAK3"),
  lfc=c(0.123, -0.320, -0.169, tg[gene=="JAK3",lfc]),
  FDR=c(0.44, 6.1e-8, 0.33, tg[gene=="JAK3",trend_FDR]),
  drug=c("deucravacitinib","JAK1i (upadacitinib etc.)","JAK2i","JAK3/pan-JAKi (tofacitinib)"))
sig_nodes$sig <- sig_nodes$FDR < 0.05
sig_nodes$gene_f <- factor(sig_nodes$gene, levels=sig_nodes$gene[order(sig_nodes$lfc)])

pA <- ggplot(act, aes(npn, gene_f)) +
  geom_vline(xintercept=c(0,1), linetype=3, color=GREY) +
  geom_vline(xintercept=0.15, linetype=2, color="#999999") +
  geom_segment(aes(x=0, xend=npn, yend=gene_f, color=grp), linewidth=0.5, alpha=0.5) +
  geom_point(aes(color=grp, size=lfc)) +
  geom_point(data=act[is_stat3==TRUE], shape=21, size=5, stroke=1.2, color="black", fill=NA) +
  scale_color_manual(values=c("late_PP (plaque-specific)"=RED, "progressive (rises through PN)"="#4C72B0"), name=NULL) +
  scale_size_continuous(range=c(2,7), name="PP-NN log2FC") +
  labs(x="PN position on healthy(0) -> lesional(1) axis", y=NULL,
       title="A  Druggable targets on the psoriasis gradient",
       subtitle="Dashed line ~0.15 splits late/plaque-specific (left) from rises-through-PN (right); STAT3 circled") +
  coord_cartesian(xlim=c(-0.05,0.5)) +
  theme(legend.position="right")

pB <- ggplot(sig_nodes, aes(lfc, gene_f, fill=sig)) +
  geom_vline(xintercept=0, color=GREY) +
  geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%+.2f%s", lfc, ifelse(sig,"*"," n.s."))),
            hjust=ifelse(sig_nodes$lfc>0,-0.1,1.1), size=3) +
  geom_text(aes(x=0, label=paste0("  ",drug)), hjust=0, vjust=-1.6, size=2.6, color="#555555") +
  scale_fill_manual(values=c(`TRUE`="#55A868",`FALSE`="#CCCCCC"), name="FDR<0.05") +
  labs(x="STAT3-pathway kinase transcript, PP-NN log2FC", y=NULL,
       title="B  Most JAK/TYK2 drug nodes track disease poorly at the transcript level",
       subtitle="Exception JAK3 (+2.04); but TYK2 (deucravacitinib) & JAK2 barely move and JAK1 falls - these are inhibited at the ACTIVITY level, not abundance") +
  coord_cartesian(xlim=c(-1.1,3.0))

fig <- pA / pB + plot_layout(heights=c(1.3,0.8))
ggsave("fig39_druggable_targets_gradient.png", fig, width=11, height=9, dpi=150)