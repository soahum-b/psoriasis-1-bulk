# Auto-extracted generating script
# Produces: gsea_focus_stat3_axis.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_ranks.rds
# Source artifact version: 21ba01b2-4e12-4a96-8ed2-0bee6e56ea15
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
})

ranks <- readRDS("gsea_ranks.rds")

set.seed(42)
build_list <- function(coll, sub=NULL) {
  df <- if (is.null(sub)) msigdbr(species="human", collection=coll)
        else msigdbr(species="human", collection=coll, subcollection=sub)
  split(df$gene_symbol, df$gs_name)
}
run_gsea <- function(lst) {
  r <- fgsea(pathways=lst, stats=ranks, minSize=10, maxSize=500, eps=0)
  r[order(r$padj), ]
}

H_list     <- build_list("H")
React_list <- build_list("C2", "CP:REACTOME")
GOBP_list  <- build_list("C5", "GO:BP")
KEGG_list  <- build_list("C2", "CP:KEGG_LEGACY")

fg_H     <- run_gsea(H_list)
fg_React <- run_gsea(React_list)
fg_GOBP  <- run_gsea(GOBP_list)
fg_KEGG  <- run_gsea(KEGG_list)

focus <- function(dt, coll) {
  d <- as.data.table(dt)[, .(pathway, NES, padj, size)]
  d[, collection := coll]
  d
}
allp <- rbindlist(list(focus(fg_H,"Hallmark"), focus(fg_React,"Reactome"),
                       focus(fg_GOBP,"GO:BP"), focus(fg_KEGG,"KEGG")))

pat <- "STAT|IL17|IL_17|INTERLEUKIN_17|JAK|IL6|IL23|KERATINOCYTE|KERATINIZATION|CORNIF|NF_?KAPPA|NFKB|TNF"
hits <- allp[grepl(pat, pathway, ignore.case=TRUE) & padj<0.05][order(-NES)]

write.csv(as.data.frame(hits), "gsea_focus_stat3_axis.csv", row.names=FALSE)