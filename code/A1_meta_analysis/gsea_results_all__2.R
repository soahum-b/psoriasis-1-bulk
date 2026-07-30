# Auto-extracted generating script
# Produces: gsea_results_all.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_ranks.rds
# Source artifact version: c62a4b3f-8264-4b9b-8e81-29b537a930a6
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(fgsea)
  library(msigdbr)
})

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")

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

H <- msigdbr(species="human", collection="H")
H_list <- split(H$gene_symbol, H$gs_name)

set.seed(42)
fg_H <- fgsea(pathways=H_list, stats=ranks, minSize=10, maxSize=500, eps=0)
fg_H <- fg_H[order(fg_H$padj), ]

React_list <- build_list("C2", "CP:REACTOME")
GOBP_list  <- build_list("C5", "GO:BP")
KEGG_list  <- build_list("C2", "CP:KEGG_LEGACY")

fg_React <- run_gsea(React_list)
fg_GOBP  <- run_gsea(GOBP_list)
fg_KEGG  <- run_gsea(KEGG_list)

saveRDS(list(H=fg_H, Reactome=fg_React, GOBP=fg_GOBP, KEGG=fg_KEGG),
        "gsea_results_all.rds")