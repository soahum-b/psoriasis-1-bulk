# Auto-extracted generating script
# Produces: gsea_perstudy_raw_k4.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: b23eb2af-d867-42ae-946e-e56962bc3444
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(recount3); library(edgeR); library(limma); library(msigdbr); library(data.table); library(jsonlite)})

studies <- c("SRP035988","SRP165679","SRP126422","SRP065812")

classify <- function(srp, md){
  a <- md$sra.sample_attributes; st <- md$sra.sample_title; ext <- md$external_id
  cls <- rep("DROP", nrow(md))
  if(srp=="SRP035988"){ cls[grepl("lesional psoriatic",a)]<-"PP"; cls[grepl("normal skin",a)]<-"NN" }
  else if(srp=="SRP165679"){
    dis<-ifelse(grepl("^PSO",st),"PSO",ifelse(grepl("^AD",st),"AD","HC"))
    skin<-ifelse(grepl("non-lesional",st),"PN",ifelse(grepl("lesional",st),"PP","NN"))
    cls[dis=="PSO"&skin=="PP"]<-"PP"; cls[dis=="HC"]<-"NN"
  } else if(srp=="SRP126422"){
    sb<-grepl("skin biopsy",a,ignore.case=TRUE)
    cls[sb&grepl("\\(PP\\)",a)]<-"PP"; cls[sb&grepl("\\(NN\\)",a)]<-"NN"
  } else if(srp=="SRP065812"){
    cls[grepl("before adalimumab",a,ignore.case=TRUE) & grepl("lesional",a,ignore.case=TRUE)]<-"PP"
    cls[grepl("healthy control",a,ignore.case=TRUE) | grepl("normal skin",a,ignore.case=TRUE)]<-"NN"
    cls[grepl("after|post",a,ignore.case=TRUE)]<-"DROP"
  }
  data.table(external_id=ext, srp=srp, class=cls)
}
load_md <- function(srp) read_metadata(file_retrieve(locate_url(srp,"data_sources/sra","metadata",organism="human")))
get_rse <- function(srp){
  rse <- create_rse(subset(available_projects(organism="human"),
                    project==srp & project_type=="data_sources"), type="gene")
  assay(rse,"counts") <- transform_counts(rse); rse
}

de_one <- function(srp){
  md <- load_md(srp); cl <- classify(srp, md); rse <- get_rse(srp)
  cvec <- cl$class[match(colnames(rse), cl$external_id)]
  keep_s <- which(cvec %in% c("NN","PP"))
  grp <- factor(cvec[keep_s], levels=c("NN","PP"))
  dge <- DGEList(counts=assay(rse,"counts")[,keep_s], genes=data.frame(gene_name=rowData(rse)$gene_name))
  keep <- filterByExpr(dge, group=grp); dge<-dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge,"TMM")
  design <- model.matrix(~grp); v <- voom(dge, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
  dt <- as.data.table(tt); dt[, gene := dge$genes$gene_name]
  dt <- dt[order(-abs(t))][!duplicated(gene)]
  dt[, .(gene, logFC, t, P.Value, study=srp)]
}

set.seed(1)
per_study <- rbindlist(lapply(studies, de_one))
saveRDS(per_study, "per_study_de_k4_primary.rds")

build_list <- function(coll, sub=NULL){
  df <- if(is.null(sub)) msigdbr(species="Homo sapiens", collection=coll)
        else msigdbr(species="Homo sapiens", collection=coll, subcollection=sub)
  split(df$gene_symbol, df$gs_name)
}
allsets <- c(build_list("H"), build_list("C2","CP:REACTOME"),
             build_list("C5","GO:BP"), build_list("C2","CP:KEGG_LEGACY"))
allsets <- allsets[!duplicated(names(allsets))]

studies <- unique(per_study$study)
gsea_per <- rbindlist(lapply(studies, function(s){
  d <- per_study[study==s & is.finite(t) & !is.na(gene)][!duplicated(gene)]
  rk <- setNames(d$t, d$gene); rk <- sort(rk, decreasing=TRUE)
  set.seed(42)
  r <- fgsea(pathways=allsets, stats=rk, minSize=10, maxSize=500, eps=0)
  as.data.table(r)[, .(pathway, study=s, NES, pval, size)]
}))

fwrite(gsea_per[, .(pathway, study, NES, pval, size)], "gsea_perstudy_raw_k4.csv")