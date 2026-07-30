# Auto-extracted generating script
# Produces: stat3_5v7study_comparison.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 21da016f-8ba2-47cd-a104-af1a6bd8ceb4
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(recount3); library(SummarizedExperiment)})
suppressMessages({library(edgeR); library(limma); library(metafor)})
suppressMessages({library(ggplot2); library(patchwork); library(data.table)})

theme_set(theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=12), plot.title.position="plot"))

ap <- available_projects(organism="human")
get_rse <- function(srp){
  pj <- subset(ap, project==srp & project_type=="data_sources")
  rse <- create_rse(pj); assay(rse,"counts") <- transform_counts(rse); rse
}
rse_s3 <- get_rse("SRP154474")
rse_s4 <- get_rse("SRP057087")

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

at3 <- tolower(as.character(colData(rse_s3)$sra.sample_attributes))
cls_s3 <- rep(NA_character_, length(at3))
cls_s3[grepl("healthy control", at3)]              <- "NN"
cls_s3[grepl("conventional psoriatic skin", at3)]  <- "PP"
cls_s3[grepl("scalp psoriatic skin", at3)]         <- "SCALP"
cls_s3[grepl("palmoplantar", at3)]                 <- "PALMO"

at4 <- tolower(as.character(colData(rse_s4)$sra.sample_attributes))
cls_s4 <- rep(NA_character_, length(at4))
cls_s4[grepl("group;;pp", at4)] <- "PP"
cls_s4[grepl("group;;pn", at4)] <- "PN"

de_s3_PPvsNN <- de_one(rse_s3, cls_s3, "PP","NN")
de_s4_PNvsPP <- de_one(rse_s4, cls_s4, "PN","PP")

psd_old <- readRDS("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/per_study_de.rds")
meta_old <- readRDS("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/meta_de_results.rds")

psd_new <- psd_old
psd_new$PPvsNN$SRP154474 <- de_s3_PPvsNN
psd_new$PNvsPP$SRP057087 <- de_s4_PNvsPP

stat3_rows <- function(psd){
  do.call(rbind, lapply(names(psd$PPvsNN), function(s){
    d<-psd$PPvsNN[[s]]; r<-d[d$gene=="STAT3",][1,]
    data.frame(study=s, yi=r$logFC, sei=r$SE, n1=r$n1, n2=r$n2)
  }))
}
s3_old <- stat3_rows(psd_old); s3_new <- stat3_rows(psd_new)
r_old <- rma(yi=yi, sei=sei, data=s3_old, method="DL")
r_new <- rma(yi=yi, sei=sei, data=s3_new, method="DL")
rk_old <- rma(yi=yi, sei=sei, data=s3_old, method="REML", test="knha")
rk_new <- rma(yi=yi, sei=sei, data=s3_new, method="REML", test="knha")

stat3_cmp <- data.frame(
  pooling=c("DL random-effects","HKSJ (REML, few-study honest)"),
  study_set_5=c(sprintf("%+.2f [%.2f, %.2f], p=%.1e", r_old$b,r_old$ci.lb,r_old$ci.ub,r_old$pval),
                sprintf("%+.2f [%.2f, %.2f], p=%.3f", rk_old$b,rk_old$ci.lb,rk_old$ci.ub,rk_old$pval)),
  study_set_7=c(sprintf("%+.2f [%.2f, %.2f], p=%.1e", r_new$b,r_new$ci.lb,r_new$ci.ub,r_new$pval),
                sprintf("%+.2f [%.2f, %.2f], p=%.3f", rk_new$b,rk_new$ci.lb,rk_new$ci.ub,rk_new$pval)))
write.csv(stat3_cmp, "stat3_5v7study_comparison.csv", row.names=FALSE)