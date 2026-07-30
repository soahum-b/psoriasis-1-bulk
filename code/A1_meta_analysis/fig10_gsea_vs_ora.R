# Auto-extracted generating script
# Produces: fig10_gsea_vs_ora.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds, gsea_ranks.rds
# Source artifact version: ee1391ed-7338-4fa7-8fc4-6ba07e291a18
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(data.table)
})

ranks <- readRDS("gsea_ranks.rds")
de    <- readRDS("de_results_full.rds")

sig_genes  <- de$gene[de$adj.P.Val < 0.05 & abs(de$logFC) > 1]
up_genes   <- de$gene[de$adj.P.Val < 0.05 & de$logFC >  1]
down_genes <- de$gene[de$adj.P.Val < 0.05 & de$logFC < -1]
universe   <- de$gene

build_list <- function(coll, sub=NULL) {
  df <- if (is.null(sub)) msigdbr(species="human", collection=coll)
        else msigdbr(species="human", collection=coll, subcollection=sub)
  split(df$gene_symbol, df$gs_name)
}

H_list     <- build_list("H")
React_list <- build_list("C2", "CP:REACTOME")
GOBP_list  <- build_list("C5", "GO:BP")
KEGG_list  <- build_list("C2", "CP:KEGG_LEGACY")

set.seed(42)
run_gsea <- function(lst) {
  r <- fgsea(pathways=lst, stats=ranks, minSize=10, maxSize=500, eps=0)
  r[order(r$padj), ]
}
fg_H     <- run_gsea(H_list)
fg_React <- run_gsea(React_list)
fg_GOBP  <- run_gsea(GOBP_list)
fg_KEGG  <- run_gsea(KEGG_list)

t2g <- function(lst) {
  do.call(rbind, lapply(names(lst), function(n) data.frame(term=n, gene=lst[[n]])))
}
T_H <- t2g(H_list); T_R <- t2g(React_list); T_G <- t2g(GOBP_list); T_K <- t2g(KEGG_list)

run_ora <- function(t2g_df, genes, uni) {
  as.data.frame(enricher(gene=genes, universe=uni, TERM2GENE=t2g_df,
                         pvalueCutoff=1, qvalueCutoff=1, minGSSize=10, maxGSSize=500))
}
ora_list <- list(
  Hallmark = run_ora(T_H, sig_genes, universe),
  Reactome = run_ora(T_R, sig_genes, universe),
  `GO:BP`  = run_ora(T_G, sig_genes, universe),
  KEGG     = run_ora(T_K, sig_genes, universe)
)

ora_std <- function(d, coll) {
  if (nrow(d)==0) return(data.table())
  gr <- as.numeric(sub("/.*","",d$GeneRatio)); gd <- as.numeric(sub(".*/","",d$GeneRatio))
  br <- as.numeric(sub("/.*","",d$BgRatio));   bd <- as.numeric(sub(".*/","",d$BgRatio))
  data.table(pathway=d$ID, collection=coll,
             ora_padj=d$p.adjust, ora_fold=(gr/gd)/(br/bd), ora_k=gr, ora_setsize=br)
}
ORA <- rbindlist(lapply(names(ora_list), function(n) ora_std(ora_list[[n]], n)), fill=TRUE)

gsea_std <- function(dt, coll) as.data.table(dt)[, .(pathway, collection=coll,
                              gsea_nes=NES, gsea_padj=padj, gsea_size=size)]
GSEA <- rbindlist(list(gsea_std(fg_H,"Hallmark"), gsea_std(fg_React,"Reactome"),
                       gsea_std(fg_GOBP,"GO:BP"), gsea_std(fg_KEGG,"KEGG")))

M <- merge(GSEA, ORA, by=c("pathway","collection"), all=TRUE)
M[, gsea_sig := !is.na(gsea_padj) & gsea_padj<0.05 & gsea_nes>0]
M[, ora_sig  := !is.na(ora_padj)  & ora_padj<0.05]

NNblue <- "#4C72B0"; PPred <- "#C44E52"; grey <- "#B0B0B0"

D <- M[!is.na(gsea_padj) & !is.na(ora_padj) & gsea_nes>0]
D[, glp := -log10(pmax(gsea_padj,1e-50))]
D[, olp := -log10(pmax(ora_padj,1e-50))]
D[, cls := fifelse(gsea_sig & ora_sig,"both",
            fifelse(gsea_sig & !ora_sig,"GSEA only",
            fifelse(!gsea_sig & ora_sig,"ORA only","neither")))]

focus_sets <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_INTERLEUKIN_17_SIGNALING","REACTOME_STAT3_NUCLEAR_EVENTS_DOWNSTREAM_OF_ALK_SIGNALING_AND_PTK6",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY","GOBP_KERATINIZATION","GOBP_KERATINOCYTE_DIFFERENTIATION",
  "GOBP_RESPONSE_TO_INTERLEUKIN_17")
lab <- D[pathway %in% focus_sets]
lab[, short := gsub("^(HALLMARK|REACTOME|KEGG|GOBP)_","",pathway)]
lab[, short := gsub("_"," ",short)]

p <- ggplot(D, aes(olp, glp)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour=grey, linewidth=0.4) +
  geom_vline(xintercept=-log10(0.05), linetype="dashed", colour=grey, linewidth=0.4) +
  geom_point(aes(colour=cls), alpha=0.55, size=1.6) +
  geom_point(data=lab, colour="black", size=2.6, shape=21, fill=PPred, stroke=0.6) +
  ggrepel::geom_text_repel(data=lab, aes(label=short), size=3.0, max.overlaps=20,
                           box.padding=0.5, min.segment.length=0, colour="grey15", seed=1) +
  scale_colour_manual(values=c(both=PPred,`GSEA only`=NNblue,`ORA only`="#DD8452",neither=grey),
                      name=NULL) +
  annotate("text", x=0.4, y=45, label="GSEA-significant\n(NES>0), ORA-null", hjust=0, size=3.2, colour=NNblue) +
  labs(title="GSEA vs ORA - two views of the same pathways",
       subtitle="Each point = one up-enriched gene set tested by both methods. Dashed = padj 0.05.",
       x=expression(ORA~-log[10]~adj.P~(Fisher)),
       y=expression(GSEA~-log[10]~adj.P~(fgsea))) +
  theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
       plot.title=element_text(face="bold"), legend.position="top")

ggsave("fig10_gsea_vs_ora.png", p, width=8.5, height=6.5, dpi=150)