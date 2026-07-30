# Auto-extracted generating script
# Produces: stat3_exon23_zoom.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: a1322047-4c8f-4049-9fd9-595753f3907d
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(Gviz); library(ensembldb); library(EnsDb.Hsapiens.v86)
  library(data.table); library(GenomicRanges)
})
options(ucscChromosomeNames=FALSE)

edb <- EnsDb.Hsapiens.v86

gr <- genes(edb, filter=GeneNameFilter("STAT3"))["ENSG00000168610"]
chr <- as.character(seqnames(gr)); gstart <- start(gr); gend <- end(gr)

col_alpha <- "#B2182B"; col_beta <- "#2166AC"

ex <- exonsBy(edb, by="tx", filter=GeneNameFilter("STAT3"))

build_df <- function(txids, symbols){
  do.call(rbind, lapply(seq_along(txids), function(i){
    g <- ex[[txids[i]]]
    if (is.null(g)) return(NULL)
    data.frame(chromosome=paste0("chr",as.character(seqnames(g))),
               start=start(g), end=end(g), width=width(g),
               strand=as.character(strand(g)),
               feature="protein_coding",
               gene="STAT3", exon=g$exon_id,
               transcript=symbols[i], symbol=symbols[i],
               stringsAsFactors=FALSE)
  }))
}
df_ab <- build_df(c("ENST00000264657","ENST00000585517"), c("STAT3\u03b1","STAT3\u03b2"))
df_ab$chromosome <- chr

grt_a_z <- GeneRegionTrack(df_ab[df_ab$transcript=="STAT3\u03b1",], chromosome=chr, genome="hg38",
             name="STAT3\u03b1", fill=col_alpha, col="black", lwd=0.5,
             col.title="black", background.title="white", rotation.title=0)
grt_b_z <- GeneRegionTrack(df_ab[df_ab$transcript=="STAT3\u03b2",], chromosome=chr, genome="hg38",
             name="STAT3\u03b2", fill=col_beta, col="black", lwd=0.5,
             col.title="black", background.title="white", rotation.title=0)
gax_z <- GenomeAxisTrack(fontcolor="black", col="black", fontsize=10)

ht <- HighlightTrack(trackList=list(grt_a_z, grt_b_z), chromosome=chr,
                     start=42348519, end=42348540, col="#F4A582", fill="#FDE0D2", inBackground=TRUE)
png("stat3_exon23_zoom.png", width=1400, height=380, res=150)
plotTracks(list(gax_z, ht), from=42348250, to=42348650, chromosome=chr, reverseStrand=TRUE,
           main="Exon 23: \u03b1 extends 19 bp (shaded) \u2192 full transactivation domain; \u03b2 truncated",
           cex.main=0.78, margin=12)
dev.off()
cat("zoom rendered\n")