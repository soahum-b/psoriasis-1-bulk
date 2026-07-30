# Auto-extracted generating script
# Produces: stat3_headline.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): stat3_gradient_isoform_vs_gene.csv, trend_SRP165679.rds
# Source artifact version: 9b472df1-1556-4cdc-8d47-2944f80247e4
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(Gviz)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(grid)
  library(cowplot)
})

options(ucscChromosomeNames=FALSE)

ap <- function(vid) host$artifact_path(vid)
edb <- EnsDb.Hsapiens.v86

tr <- readRDS(ap("d61263b9-9355-4e92-97d8-4987bcb38e86")); setDT(tr)
psi <- read.csv(ap("087461f6-2a17-486a-bb9b-4f76963e8035"))
s3 <- tr[gene=="STAT3"]

gstart <- 42313324; gend <- 42388568
col_alpha <- "#B2182B"; col_beta <- "#2166AC"
grp_cols <- c(NN="#4DAF4A", PN="#FF7F00", PP="#984EA3")

exb <- exonsBy(edb, by="tx", filter=GeneNameFilter("STAT3"))
ex_df <- function(txid, lab){
  d <- as.data.frame(exb[[txid]], row.names=NULL)
  data.frame(chromosome="17", start=d$start, end=d$end, width=d$width, strand=as.character(d$strand),
             feature="protein_coding", gene="STAT3", exon=paste0(lab, seq_len(nrow(d))),
             transcript=lab, symbol=lab, stringsAsFactors=FALSE)
}
df_ab <- rbind(ex_df("ENST00000264657", "STAT3\u03b1"), ex_df("ENST00000585517", "STAT3\u03b2"))
mk <- function(lab, fill) GeneRegionTrack(df_ab[df_ab$transcript==lab,], chromosome="17", genome="hg38",
   name=lab, fill=fill, col="black", lwd=0.3, col.title="black", background.title="white", rotation.title=0)
grt_a <- mk("STAT3\u03b1", col_alpha); grt_b <- mk("STAT3\u03b2", col_beta)
win_s <- gstart - 1500; win_e <- gend + 1500

d1 <- data.frame(arm=factor(c("NN","PN","PP"), levels=c("NN","PN","PP")), expr=c(s3$mu_NN, s3$mu_PN, s3$mu_PP))
d2 <- data.frame(arm=factor(psi$arm, levels=c("NN","PN","PP")), psi=psi$PSI_beta_SRP165679_pct)
p1 <- ggplot(d1, aes(arm, expr, fill=arm)) + geom_col(width=.65) +
  geom_text(aes(label=sprintf("%.2f", expr)), vjust=-0.4, size=3.2) +
  scale_fill_manual(values=grp_cols, guide="none") + coord_cartesian(ylim=c(0, max(d1$expr)*1.12)) +
  labs(title="STAT3 gene expression rises with lesion stage", subtitle="mean log2 expression, SRP165679 gradient cohort", x=NULL, y="log2 expression") +
  theme_classic(base_size=11) + theme(plot.title=element_text(size=10.5, face="bold"), plot.subtitle=element_text(size=8.5), axis.text=element_text(color="black"))
p2 <- ggplot(d2, aes(arm, psi, fill=arm)) + geom_col(width=.65) +
  geom_text(aes(label=sprintf("%.1f%%", psi)), vjust=-0.4, size=3.2) +
  scale_fill_manual(values=grp_cols, guide="none") + coord_cartesian(ylim=c(0, max(d2$psi)*1.18)) +
  labs(title="STAT3\u03b2 isoform fraction does NOT track the gene", subtitle="PSI (percent-spliced-in) of \u03b2, per group", x=NULL, y="STAT3\u03b2 PSI (%)") +
  theme_classic(base_size=11) + theme(plot.title=element_text(size=10.5, face="bold"), plot.subtitle=element_text(size=8.5), axis.text=element_text(color="black"))

gax <- GenomeAxisTrack(fontcolor="black", col="black", fontsize=11)
full_grob <- grid.grabExpr(plotTracks(list(gax, grt_a, grt_b), from=win_s, to=win_e, chromosome="17",
  reverseStrand=TRUE, sizes=c(1.2, 1, 1), main="STAT3 locus (chr17, \u2212 strand): \u03b1 (full TAD) vs \u03b2 (truncated C-terminus)", cex.main=0.85, margin=14, add=TRUE))

cds <- cdsBy(edb, by="tx", filter=GeneNameFilter("STAT3"))
cds_df <- function(txid, lab){
  d <- as.data.frame(cds[[txid]], row.names=NULL)
  data.frame(chromosome="17", start=d$start, end=d$end, width=d$width,
             strand=as.character(d$strand), feature="CDS",
             gene="STAT3", exon=paste0(lab, "_", seq_len(nrow(d))),
             transcript=lab, symbol=lab, stringsAsFactors=FALSE)
}
df_c <- rbind(cds_df("ENST00000264657", "STAT3\u03b1"), cds_df("ENST00000585517", "STAT3\u03b2"))
grt_a_c <- GeneRegionTrack(df_c[df_c$transcript=="STAT3\u03b1",], chromosome="17", genome="hg38", name="STAT3\u03b1", fill=col_alpha, col="black", lwd=0.5, col.title="black", background.title="white", rotation.title=0)
grt_b_c <- GeneRegionTrack(df_c[df_c$transcript=="STAT3\u03b2",], chromosome="17", genome="hg38", name="STAT3\u03b2", fill=col_beta, col="black", lwd=0.5, col.title="black", background.title="white", rotation.title=0)
gax_z <- GenomeAxisTrack(fontcolor="black", col="black", fontsize=9)
ht <- HighlightTrack(trackList=list(grt_a_c, grt_b_c), chromosome="17", start=42315745, end=42316827, col="#F4A582", fill="#FDE0D2", inBackground=TRUE)
zoom_grob <- grid.grabExpr(plotTracks(list(gax_z, ht), from=42315400, to=42317350, chromosome="17", reverseStrand=TRUE,
  main="C-terminus (coding): \u03b1 extends ~1.1 kb further (shaded) \u2192 full transactivation domain; \u03b2 stops early", cex.main=0.72, margin=12, add=TRUE))

top <- plot_grid(full_grob, zoom_grob, ncol=1, rel_heights=c(1, 0.9), labels=c("a","b"), label_size=14)
bottom <- plot_grid(p1, p2, ncol=2, labels=c("c","d"), label_size=14)
lead <- plot_grid(top, bottom, ncol=1, rel_heights=c(1.9, 1))
ggsave("stat3_headline.png", lead, width=11, height=9, dpi=170)