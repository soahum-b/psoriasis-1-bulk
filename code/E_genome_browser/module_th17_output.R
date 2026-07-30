# Auto-extracted generating script
# Produces: module_th17_output.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): trend_SRP165679.rds
# Source artifact version: 25993623-d4e0-4119-9a07-ea71bc6db5b6
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(Gviz); library(ensembldb); library(EnsDb.Hsapiens.v86)
  library(data.table); library(GenomicRanges); library(org.Hs.eg.db)
  library(ggplot2); library(grid); library(cowplot)
})
options(ucscChromosomeNames=FALSE)
ap <- function(vid) vid

edb <- EnsDb.Hsapiens.v86

tr <- readRDS("trend_SRP165679.rds"); setDT(tr)

grp_cols <- c(NN="#4DAF4A", PN="#FF7F00", PP="#984EA3")

model_grob <- function(sym){
  g <- genes(edb, filter=GeneNameFilter(sym)); g <- g[g$gene_biotype=="protein_coding"]
  if(length(g)==0) return(NULL)
  g <- g[which.max(width(g))]
  cc<-as.character(seqnames(g)); s<-start(g); e<-end(g); st<-as.character(strand(g))
  pad <- max(1000, round(width(g)*0.03))
  grt <- GeneRegionTrack(edb, chromosome=cc, start=s-pad, end=e+pad,
           filter=GeneNameFilter(sym), name=sym, genome="hg38",
           collapseTranscripts="meta", transcriptAnnotation="none",
           fill="grey55", col="black", lwd=0.4, col.title="black",
           background.title="white", rotation.title=0, cex.title=0.8, fontface.title=3)
  gax <- GenomeAxisTrack(fontcolor="black", col="black", fontsize=8, cex=0.7)
  grid.grabExpr(plotTracks(list(gax, grt), from=s-pad, to=e+pad, chromosome=cc,
           reverseStrand=(st=="-"), sizes=c(0.9,1),
           title.width=1.1, margin=6, add=TRUE))
}

bar_grob <- function(sym){
  row <- tr[gene==sym]; if(nrow(row)==0) return(NULL)
  d <- data.frame(arm=factor(c("NN","PN","PP"),levels=c("NN","PN","PP")),
                  val=c(0, row$mu_PN-row$mu_NN, row$mu_PP-row$mu_NN))
  ggplot(d, aes(arm,val,fill=arm))+geom_col(width=.7)+
    geom_hline(yintercept=0,color="grey40",linewidth=.3)+
    scale_fill_manual(values=grp_cols,guide="none")+
    labs(x=NULL,y=NULL)+
    theme_classic(base_size=9)+
    theme(axis.text=element_text(color="black",size=7.5),
          plot.margin=margin(2,4,2,2))
}

module_fig <- function(genes_v, module_name, outfile){
  genes_v <- genes_v[genes_v %in% tr$gene]
  rows <- lapply(genes_v, function(sym){
    mg <- model_grob(sym); bg <- bar_grob(sym)
    if(is.null(mg)||is.null(bg)) return(NULL)
    plot_grid(mg, bg, nrow=1, rel_widths=c(2.4,1))
  })
  rows <- rows[!sapply(rows,is.null)]
  fig <- plot_grid(plotlist=rows, ncol=1)
  title <- ggdraw()+draw_label(
    sprintf("%s module: locus models (left) + per-group expression change vs NN, log2 (right)", module_name),
    fontface="bold", size=12, x=0.01, hjust=0)
  full <- plot_grid(title, fig, ncol=1, rel_heights=c(0.05, 1))
  ggsave(outfile, full, width=10, height=1.5*length(rows)+0.5, dpi=150, limitsize=FALSE)
  length(rows)
}

n <- module_fig(c("IL17A","IL23A","RORC","S100A7","S100A8","S100A9","DEFB4A","CCL20"),
                "Th17/output", "module_th17_output.png")
cat("Th17_output :", n, "rows\n")