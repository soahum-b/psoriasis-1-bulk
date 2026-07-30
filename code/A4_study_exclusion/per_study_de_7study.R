# Auto-extracted generating script
# Produces: per_study_de_7study.rds
# Conda env: psoriasis-r
# Source artifact version: 1a70e129-8c63-482b-91e1-218eb73ef156
#------------------------------------------------------------

library(recount3)
library(SummarizedExperiment)
library(edgeR)
library(limma)

# Load existing per-study DE results
psd_old <- readRDS("{{artifact:MISSING:5e70760d-61d9-4b40-b0e5-79d94a01f78e}}")

# Get RSE objects from recount3
ap <- available_projects(organism="human")
get_rse <- function(srp){
  pj <- subset(ap, project==srp & project_type=="data_sources")
  rse <- create_rse(pj); assay(rse,"counts") <- transform_counts(rse); rse
}
rse_s3 <- get_rse("SRP154474")
rse_s4 <- get_rse("SRP057087")

# DE function
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

# Run DE for new cohorts
de_s3_PPvsNN <- de_one(rse_s3, cls_s3, "PP","NN")
de_s4_PNvsPP <- de_one(rse_s4, cls_s4, "PN","PP")

# Build new 7-study per-study DE list
psd_new <- psd_old
psd_new$PPvsNN$SRP154474 <- de_s3_PPvsNN
psd_new$PNvsPP$SRP057087 <- de_s4_PNvsPP

saveRDS(psd_new, "per_study_de_7study.rds")