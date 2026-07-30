# Auto-extracted generating script
# Produces: stat3_lead.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): trend_SRP165679.rds, stat3_gradient_isoform_vs_gene.csv
# Source artifact version: d65222b8-5425-4907-af65-e5d5aad75ca0
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(Gviz)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(data.table)
  library(GenomicRanges)
})
options(ucscChromosomeNames=FALSE)

edb <- EnsDb.Hsapiens.v86

col_alpha <- "#B2182B"
col_beta <- "#2166AC"

gstart <- 42313324
gend <- 42388568
chr <- "17"

exb <- exonsBy(edb, by="tx", filter=GeneNameFilter("STAT3"))

ex_df <- function(txid, lab) {
  d <- as.data.frame(exb[[txid]], row.names=NULL)
  data.frame(chromosome="17", start=d$start, end=d$end, width=d$width,
             strand=as.character(d$strand),
             feature="protein_coding", gene="STAT3",
             exon=paste0(lab, seq_len(nrow(d))),
             transcript=lab, symbol=lab, stringsAsFactors=FALSE)
}

df_ab <- rbind(ex_df("ENST00000264657", "STAT3\u03b1"),
               ex_df("ENST00000585517", "STAT3\u03b2"))

mk <- function(lab, fill) {
  GeneRegionTrack(df_ab[df_ab$transcript==lab,], chromosome="17", genome="hg38",
                  name=lab, fill=fill, col="black", lwd=0.3,
                  col.title="black", background.title="white", rotation.title=0)
}

grt_a <- mk("STAT3\u03b1", col_alpha)
grt_b <- mk("STAT3\u03b2", col_beta)

gax <- GenomeAxisTrack(fontcolor="black", col="black", fontsize=11)

win_s <- gstart - 1500
win_e <- gend + 1500

png("stat3_lead.png", width=1500, height=430, res=150)
plotTracks(list(gax, grt_a, grt_b), from=win_s, to=win_e, chromosome="17",
           reverseStrand=TRUE, sizes=c(1.2, 1, 1),
           main="STAT3 locus (chr17, \u2212 strand): \u03b1 (full TAD) vs \u03b2 (truncated C-terminus)",
           cex.main=0.85, margin=20)
dev.off()