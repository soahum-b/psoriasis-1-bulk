# Auto-extracted generating script
# Produces: camera_old_vs_corrected.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): meta_pathway_camera_PPvsNN.csv
# Source artifact version: 3fe47dcc-50b9-4798-a37b-55ebc614933d
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(edgeR)
  library(limma)
  library(msigdbr)
  library(data.table)
})

## --- classify samples per study ---
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
    cls[grepl("lesional",a,ignore.case=TRUE) & !grepl("non-lesional|uninvolved",a,ignore.case=TRUE)]<-"PP"
    cls[grepl("normal|healthy|control",a,ignore.case=TRUE)]<-"NN"
    cls[grepl("after|post",a,ignore.case=TRUE)] <- "DROP"
  }
  data.table(external_id=ext, srp=srp, class=cls)
}

get_rse <- function(srp){
  rse <- create_rse(subset(available_projects(organism="human"),
                    project==srp & project_type=="data_sources"), type="gene")
  assay(rse,"counts") <- transform_counts(rse); rse
}
load_md <- function(srp) read_metadata(file_retrieve(locate_url(srp,"data_sources/sra","metadata",organism="human")))

## Hallmark
H <- msigdbr(species="Homo sapiens", collection="H")
Hlist <- split(H$gene_symbol, H$gs_name)

cam_one <- function(srp){
  md <- load_md(srp); cl <- classify(srp, md)
  rse <- get_rse(srp)
  cvec <- cl$class[match(colnames(rse), cl$external_id)]
  keep_s <- which(cvec %in% c("NN","PP"))
  grp <- factor(cvec[keep_s], levels=c("NN","PP"))
  if(length(unique(grp))<2 || min(table(grp))<3) return(NULL)
  dge <- DGEList(counts=assay(rse,"counts")[,keep_s],
                 genes=data.frame(gene_name=rowData(rse)$gene_name))
  keep <- filterByExpr(dge, group=grp); dge<-dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge,"TMM")
  design <- model.matrix(~grp); v <- voom(dge, design)
  idx <- ids2indices(Hlist, dge$genes$gene_name)
  cam <- camera(v, idx, design, contrast=2)
  dt <- as.data.table(cam, keep.rownames="pathway")
  dt[, `:=`(study=srp, n=length(grp), nNN=sum(grp=="NN"), nPP=sum(grp=="PP"))]
  dt
}

studies <- c("SRP035988","SRP165679","SRP126422","SRP065812")
set.seed(1)
cam_all <- rbindlist(lapply(studies, cam_one), fill=TRUE)

## Corrected classify for SRP065812
md_065812 <- load_md("SRP065812")
classify_065812 <- function(md){
  a <- md$sra.sample_attributes; st <- md$sra.sample_title; ext <- md$external_id
  cls <- rep("DROP", nrow(md))
  cls[grepl("before adalimumab",a,ignore.case=TRUE) & grepl("lesional",a,ignore.case=TRUE)] <- "PP"
  cls[grepl("healthy control",a,ignore.case=TRUE) | grepl("normal skin",a,ignore.case=TRUE)] <- "NN"
  cls[grepl("after|post",a,ignore.case=TRUE)] <- "DROP"
  data.table(external_id=ext, srp="SRP065812", class=cls)
}
cl <- classify_065812(md_065812)
rse_065812 <- get_rse("SRP065812")
cvec <- cl$class[match(colnames(rse_065812), cl$external_id)]
keep_s <- which(cvec %in% c("NN","PP"))
grp <- factor(cvec[keep_s], levels=c("NN","PP"))
dge <- DGEList(counts=assay(rse_065812,"counts")[,keep_s], genes=data.frame(gene_name=rowData(rse_065812)$gene_name))
keep <- filterByExpr(dge, group=grp); dge<-dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge,"TMM")
design <- model.matrix(~grp); v <- voom(dge, design)
idx <- ids2indices(Hlist, dge$genes$gene_name)
cam <- camera(v, idx, design, contrast=2)
dt <- as.data.table(cam, keep.rownames="pathway")
dt[, `:=`(study="SRP065812", n=length(grp), nNN=sum(grp=="NN"), nPP=sum(grp=="PP"))]

cam_all <- cam_all[study!="SRP065812"]
cam_all <- rbindlist(list(cam_all, dt), fill=TRUE)

## Compute corrected pooled statistics
setDT(cam_all)
cam_all[, p_use := pmin(pmax(PValue, .Machine$double.xmin), 1)]
cam_all[, sgn := ifelse(Direction=="Up", 1, -1)]
cam_all[, logp := log(pmax(PValue, .Machine$double.xmin))]
cam_all[, z := sgn * qnorm(logp - log(2), lower.tail=FALSE, log.p=TRUE)]
cam_all[, w := sqrt(n)]

pooled <- cam_all[, {
  Zc  <- sum(w*z)/sqrt(sum(w^2))
  logpS <- log(2) + pnorm(-abs(Zc), log.p=TRUE)
  Xf  <- -2*sum(logp); dfF <- 2*.N
  logpF <- pchisq(Xf, df=dfF, lower.tail=FALSE, log.p=TRUE)
  .(k=.N, nUp=sum(Direction=="Up"), nDown=sum(Direction=="Down"),
    Zc=Zc, meanNGenes=round(mean(NGenes)),
    log_p_stouffer=logpS, log_p_fisher=logpF,
    minStudyP=min(PValue), maxStudyP=max(PValue))
}, by=pathway]

pooled[, p_stouffer := exp(log_p_stouffer)]
pooled[, p_fisher   := exp(log_p_fisher)]
pooled[, FDR_stouffer := p.adjust(exp(pmax(log_p_stouffer,-745)),"BH")]
pooled[, FDR_fisher   := p.adjust(exp(pmax(log_p_fisher,-745)),"BH")]
pooled[, consistent := (nUp==k | nDown==k)]
pooled <- pooled[order(-abs(Zc))]

out <- pooled[, .(pathway, k, nUp, nDown, direction_consistent=consistent,
                  Zc=round(Zc,3), meanNGenes,
                  p_stouffer, FDR_stouffer, p_fisher, FDR_fisher,
                  minStudyP, maxStudyP)]

old_cam <- fread("meta_pathway_camera_PPvsNN.csv")
cmp <- merge(old_cam[, .(pathway, k_old=k, Zc_old=Zc, p_old=p_comb, FDR_old=FDR)],
             out[, .(pathway, k_new=k, Zc_new=Zc, p_new=p_stouffer, FDR_new=FDR_stouffer, direction_consistent)],
             by="pathway", all=TRUE)
cmp <- cmp[order(-abs(Zc_new))]
fwrite(cmp, "camera_old_vs_corrected.csv")