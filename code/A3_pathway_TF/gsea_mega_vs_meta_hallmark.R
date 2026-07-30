# Auto-extracted generating script
# Produces: gsea_mega_vs_meta_hallmark.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): paths.json
# Source artifact version: 4ab384c7-c2f5-480d-a288-32efe13c5883
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(fgsea); library(msigdbr); library(data.table); library(jsonlite)})

p <- fromJSON("paths.json")

build_list <- function(coll, sub=NULL){
  df <- if(is.null(sub)) msigdbr(species="Homo sapiens", collection=coll)
        else msigdbr(species="Homo sapiens", collection=coll, subcollection=sub)
  split(df$gene_symbol, df$gs_name)
}
colls <- list(
  Hallmark = build_list("H"),
  Reactome = build_list("C2","CP:REACTOME"),
  `GO:BP`  = build_list("C5","GO:BP"),
  KEGG     = build_list("C2","CP:KEGG_LEGACY")
)

make_ranks <- function(dt){
  dt <- dt[is.finite(logFC) & is.finite(SE) & SE>0 & !is.na(gene)]
  dt <- dt[!duplicated(gene)]
  z <- dt$logFC/dt$SE
  names(z) <- dt$gene
  sort(z, decreasing=TRUE)
}
cohorts <- list(
  k4_primary  = fread(p[["meta_de_PPvsNN_4study.csv"]]),
  k4_extended = fread(p[["meta_de_PPvsNN_7study.csv"]])
)
ranks <- lapply(cohorts, make_ranks)

run_all <- function(rk){
  set.seed(42)
  rbindlist(lapply(names(colls), function(cn){
    r <- fgsea(pathways=colls[[cn]], stats=rk, minSize=10, maxSize=500, eps=0)
    as.data.table(r)[, collection:=cn]
  }))
}
gsea_meta <- lapply(ranks, run_all)
for(nm in names(gsea_meta)){
  g <- gsea_meta[[nm]]
  g[, leadingEdge := sapply(leadingEdge, function(x) paste(head(x,30), collapse=";"))]
  g <- g[order(padj, -abs(NES))]
  fwrite(g, paste0("gsea_meta_", nm, ".csv"))
}

old <- fread(p[["gsea_results_all.csv"]])
new <- fread("gsea_meta_k4_primary.csv")

oh <- old[collection=="Hallmark", .(pathway, NES_mega=NES, padj_mega=padj)]
nh <- new[collection=="Hallmark", .(pathway, NES_meta=NES, padj_meta=padj)]
cmp <- merge(oh, nh, by="pathway")
cmp[, `:=`(nlp_mega=-log10(pmax(padj_mega,1e-320)), nlp_meta=-log10(pmax(padj_meta,1e-320)))]
cmp <- cmp[order(-nlp_mega)]
fwrite(cmp, "gsea_mega_vs_meta_hallmark.csv")