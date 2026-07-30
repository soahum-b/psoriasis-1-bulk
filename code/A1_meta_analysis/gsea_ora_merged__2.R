# Auto-extracted generating script
# Produces: gsea_ora_merged.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds, gsea_ranks.rds, dge_filt_norm.rds
# Source artifact version: 9f68c49d-2d51-4a63-83e3-7990b6aedcf5
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(data.table)
})

de <- readRDS("de_results_full.rds")
ranks <- readRDS("gsea_ranks.rds")

sig_genes <- de$gene[de$adj.P.Val < 0.05 & abs(de$logFC) > 1]
universe  <- de$gene

build_list <- function(coll, sub=NULL) {
  df <- if (is.null(sub)) msigdbr(species="human", collection=coll)
        else msigdbr(species="human", collection=coll, subcollection=sub)
  split(df$gene_symbol, df$gs_name)
}

H_list     <- build_list("H")
React_list <- build_list("C2", "CP:REACTOME")
GOBP_list  <- build_list("C5", "GO:BP")
KEGG_list  <- build_list("C2", "CP:KEGG_LEGACY")

run_gsea <- function(lst) {
  r <- fgsea(pathways=lst, stats=ranks, minSize=10, maxSize=500, eps=0)
  r[order(r$padj), ]
}

set.seed(42)
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

saveRDS(M, "gsea_ora_merged.rds")
fwrite(M, "gsea_ora_merged.csv")