# Auto-extracted generating script
# Produces: gsea_perstudy_combined_k4.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 14ab83e2-2b01-4adc-ab49-46822f785e57
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(recount3); library(edgeR); library(limma); library(msigdbr); library(data.table)})

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

suppressMessages({library(fgsea)})

build_list <- function(coll, sub=NULL){
  df <- if(is.null(sub)) msigdbr(species="Homo sapiens", collection=coll)
        else msigdbr(species="Homo sapiens", collection=coll, subcollection=sub)
  split(df$gene_symbol, df$gs_name)
}
allsets <- c(build_list("H"), build_list("C2","CP:REACTOME"),
             build_list("C5","GO:BP"), build_list("C2","CP:KEGG_LEGACY"))
allsets <- allsets[!duplicated(names(allsets))]
coll_of <- function(pw) fifelse(grepl("^HALLMARK_",pw),"Hallmark",
              fifelse(grepl("^REACTOME_",pw),"Reactome",
              fifelse(grepl("^GOBP_",pw),"GO:BP",fifelse(grepl("^KEGG_",pw),"KEGG","other"))))

studies <- unique(per_study$study)
gsea_per <- rbindlist(lapply(studies, function(s){
  d <- per_study[study==s & is.finite(t) & !is.na(gene)][!duplicated(gene)]
  rk <- setNames(d$t, d$gene); rk <- sort(rk, decreasing=TRUE)
  set.seed(42)
  r <- fgsea(pathways=allsets, stats=rk, minSize=10, maxSize=500, eps=0)
  as.data.table(r)[, .(pathway, study=s, NES, pval, size)]
}))

gsea_per[, sgn := sign(NES)]
gsea_per[, p_up := ifelse(sgn>0, pval/2, 1 - pval/2)]
gsea_per[, z_up := qnorm(p_up, lower.tail=FALSE)]
gsea_per[!is.finite(z_up), z_up := ifelse(sgn>0, qnorm(pmax(p_up,1e-300),lower.tail=FALSE), -38)]

comb <- gsea_per[, {
  k <- .N
  Zc_up <- sum(z_up)/sqrt(k)
  p_stouffer_up <- pnorm(Zc_up, lower.tail=FALSE)
  Xf <- -2*sum(log(pmax(pval,1e-300))); pF <- pchisq(Xf, df=2*k, lower.tail=FALSE)
  .(k=k, meanNES=mean(NES), nUp=sum(sgn>0), Zc_up=Zc_up,
    p_stouffer_up=p_stouffer_up, p_fisher=pF, size=size[1])
}, by=pathway]
comb[, collection := coll_of(pathway)]
comb[, FDR_stouffer_up := p.adjust(p_stouffer_up,"BH")]
comb[, FDR_fisher := p.adjust(p_fisher,"BH")]
comb[, direction_consistent := (nUp==k | nUp==0)]
comb <- comb[order(p_stouffer_up)]

fwrite(comb, "gsea_perstudy_combined_k4.csv")