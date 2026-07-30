# Auto-extracted generating script
# Produces: gsea_results_all.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_ranks.rds
# Source artifact version: 6701f53b-1c6a-402a-a351-ec00f7baef8a
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

H_list     <- split(msigdbr(species="human", collection="H")$gene_symbol,
                    msigdbr(species="human", collection="H")$gs_name)
React_list <- build_list("C2", "CP:REACTOME")
GOBP_list  <- build_list("C5", "GO:BP")
KEGG_list  <- build_list("C2", "CP:KEGG_LEGACY")

set.seed(42)
fg_H     <- fgsea(pathways=H_list,     stats=ranks, minSize=10, maxSize=500, eps=0)
fg_H     <- fg_H[order(fg_H$padj), ]

set.seed(42)
fg_React <- run_gsea(React_list)
fg_GOBP  <- run_gsea(GOBP_list)
fg_KEGG  <- run_gsea(KEGG_list)

flat <- function(dt, coll) {
  d <- as.data.frame(dt)
  d$leadingEdge <- sapply(d$leadingEdge, function(x) paste(head(x,30), collapse=";"))
  d$collection <- coll
  d[, c("collection","pathway","NES","ES","padj","pval","size","leadingEdge")]
}
full <- do.call(rbind, list(flat(fg_H,"Hallmark"), flat(fg_React,"Reactome"),
                            flat(fg_GOBP,"GO:BP"), flat(fg_KEGG,"KEGG")))
write.csv(full, "gsea_results_all.csv", row.names=FALSE)