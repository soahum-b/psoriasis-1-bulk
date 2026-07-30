# Auto-extracted generating script
# Produces: key_genes_7study_PPvsNN.csv
# Conda env: psoriasis-r
# Source artifact version: 5fd6338b-66d8-45e2-bcc4-6f931f2b358d
#------------------------------------------------------------

library(recount3)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(metafor)

# Load saved per-study DE and meta results
psd_old <- readRDS("{{artifact:MISSING:5e70760d-61d9-4b40-b0e5-79d94a01f78e}}")

# Get RSEs
ap <- available_projects(organism="human")
get_rse <- function(srp){
  pj <- subset(ap, project==srp & project_type=="data_sources")
  rse <- create_rse(pj); assay(rse,"counts") <- transform_counts(rse); rse
}
rse_s3 <- get_rse("SRP154474")
rse_s4 <- get_rse("SRP057087")

# de_one function
de_one <- function(rse, cls, cA, cB){
  keep_s <- which(cls %in% c(cA,cB))
  if(length(unique(cls[keep_s]))<2) return(NULL)
  grp <- factor(cls[keep_s], levels=c(cB,cA))
  if(min(table(grp))<2) return(NULL)
  counts <- assay(rse,"counts")[,keep_s]
  gn <- rowData(rse)$gene_name
  dge <- DGEList(counts=counts, genes=data.frame(gene_name=gn))
  keep <- filterByExpr(dge, group=grp); dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge,"TMM")
  des <- model.matrix(~grp); v <- voom(dge,des); fit <- eBayes(lmFit(v,des))
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
  data.frame(gene=dge$genes$gene_name, logFC=tt$logFC, SE=tt$logFC/tt$t,
             t=tt$t, P=tt$P.Value, adjP=p.adjust(tt$P.Value,"BH"),
             AveExpr=tt$AveExpr, n1=sum(grp==cA), n2=sum(grp==cB))
}

# S3 classification
at3 <- tolower(as.character(colData(rse_s3)$sra.sample_attributes))
cls_s3 <- rep(NA_character_, length(at3))
cls_s3[grepl("healthy control", at3)]              <- "NN"
cls_s3[grepl("conventional psoriatic skin", at3)]  <- "PP"
cls_s3[grepl("scalp psoriatic skin", at3)]         <- "SCALP"
cls_s3[grepl("palmoplantar", at3)]                 <- "PALMO"

# S4 classification
at4 <- tolower(as.character(colData(rse_s4)$sra.sample_attributes))
cls_s4 <- rep(NA_character_, length(at4))
cls_s4[grepl("group;;pp", at4)] <- "PP"
cls_s4[grepl("group;;pn", at4)] <- "PN"

# Per-study DE for new cohorts
de_s3_PPvsNN <- de_one(rse_s3, cls_s3, "PP","NN")
de_s4_PNvsPP <- de_one(rse_s4, cls_s4, "PN","PP")

# meta_dl function
meta_dl <- function(de_list){
  genes <- sort(unique(unlist(lapply(de_list, function(d) d$gene))))
  Y <- sapply(de_list, function(d) d$logFC[match(genes,d$gene)])
  V <- sapply(de_list, function(d) (d$SE[match(genes,d$gene)])^2)
  rownames(Y)<-rownames(V)<-genes
  out <- t(apply(cbind(Y,V),1,function(row){
    k<-ncol(Y); yi<-row[1:k]; vi<-row[(k+1):(2*k)]
    ok<-is.finite(yi)&is.finite(vi)&vi>0
    if(sum(ok)<2) return(c(logFC=NA,SE=NA,P=NA,I2=NA,k=sum(ok)))
    yi<-yi[ok]; vi<-vi[ok]; wi<-1/vi; ybar<-sum(wi*yi)/sum(wi)
    Q<-sum(wi*(yi-ybar)^2); df<-length(yi)-1; C<-sum(wi)-sum(wi^2)/sum(wi)
    tau2<-max(0,(Q-df)/C); ws<-1/(vi+tau2); mu<-sum(ws*yi)/sum(ws); se<-sqrt(1/sum(ws))
    I2<-max(0,(Q-df)/Q)*100
    c(logFC=mu,SE=se,P=2*pnorm(-abs(mu/se)),I2=I2,k=length(yi))
  }))
  d<-data.frame(gene=genes,out); d$FDR<-p.adjust(d$P,"BH"); d
}

# Build new 7-study pool
psd_new <- psd_old
psd_new$PPvsNN$SRP154474 <- de_s3_PPvsNN
psd_new$PNvsPP$SRP057087 <- de_s4_PNvsPP

meta_new <- lapply(psd_new, meta_dl)

# Key genes table
key <- c("STAT3","IL17A","IL17F","IL23A","IL12B","SOCS3","S100A7","S100A8","S100A9",
         "DEFB4A","PI3","IL36G","CCL20","KRT16","IL19","IL36A","LCN2","CXCL8")
kt <- merge(data.frame(gene=key), meta_new$PPvsNN[,c("gene","logFC","FDR","I2","k")], by="gene", all.x=TRUE)
kt <- kt[match(key, kt$gene),]

write.csv(kt, "key_genes_7study_PPvsNN.csv", row.names=FALSE)