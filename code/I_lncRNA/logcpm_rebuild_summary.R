# Auto-extracted generating script
# Produces: logcpm_rebuild_summary.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 936a8b09-646d-4baa-b0d2-05e236c748a2
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(recount3); library(edgeR); library(limma); library(data.table)})
options(recount3_url = "http://duffel.rail.bio/recount3")
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

logcpm_one <- function(srp){
  md <- load_md(srp); cl <- classify(srp, md); rse <- get_rse(srp)
  cvec <- cl$class[match(colnames(rse), cl$external_id)]
  keep_s <- which(cvec %in% c("NN","PP"))
  grp <- factor(cvec[keep_s], levels=c("NN","PP"))
  dge <- DGEList(counts=assay(rse,"counts")[,keep_s],
                 genes=data.frame(gene_name=rowData(rse)$gene_name))
  keep <- filterByExpr(dge, group=grp); dge<-dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge,"TMM")
  lcpm <- cpm(dge, log=TRUE, prior.count=1)
  rownames(lcpm) <- dge$genes$gene_name
  o <- order(-rowMeans(lcpm)); lcpm <- lcpm[o,]; lcpm <- lcpm[!duplicated(rownames(lcpm)),]
  colnames(lcpm) <- colnames(dge)
  saveRDS(list(logcpm=lcpm, group=grp, study=srp), sprintf("logcpm_%s.rds", srp))
  data.table(study=srp, n_samples=ncol(lcpm), n_genes=nrow(lcpm),
             n_NN=sum(grp=="NN"), n_PP=sum(grp=="PP"))
}

set.seed(1)
summ <- rbindlist(lapply(studies, function(s){ cat("== ", s, " ==\n"); logcpm_one(s) }))
fwrite(summ, "logcpm_rebuild_summary.csv")
print(summ)
cat("DONE\n")