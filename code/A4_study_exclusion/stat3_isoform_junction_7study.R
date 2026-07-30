# Auto-extracted generating script
# Produces: stat3_isoform_junction_7study.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 65d1e364-b7a1-43e0-b923-e8d7bb0e1b1d
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(recount3)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(metafor)
library(data.table)
library(ggplot2)
library(patchwork)

suppressMessages({library(recount3); library(SummarizedExperiment)})
ap <- available_projects(organism="human")
get_rse <- function(srp){
  pj <- subset(ap, project==srp & project_type=="data_sources")
  rse <- create_rse(pj); assay(rse,"counts") <- transform_counts(rse); rse
}
rse_s3 <- get_rse("SRP154474")
rse_s4 <- get_rse("SRP057087")

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

acc <- 42317181L; a_don <- 42316902L; b_don <- 42316852L; tol <- 2L
psi_one <- function(rse_jxn, cls){
  rj <- rse_jxn; gr <- rowRanges(rj)
  onchr <- as.character(seqnames(gr)) %in% c("chr17","17")
  is_a <- onchr & abs(start(gr)-a_don)<=tol & abs(end(gr)-acc)<=tol
  is_b <- onchr & abs(start(gr)-b_don)<=tol & abs(end(gr)-acc)<=tol
  cat("  alpha junctions matched:", sum(is_a), " beta junctions matched:", sum(is_b), "\n")
  if(sum(is_a)<1 || sum(is_b)<1) return(NULL)
  a_cnt <- colSums(matrix(assay(rj,"counts")[is_a,,drop=FALSE], ncol=ncol(rj)))
  b_cnt <- colSums(matrix(assay(rj,"counts")[is_b,,drop=FALSE], ncol=ncol(rj)))
  data.frame(sample=colnames(rj), class=cls, a=a_cnt, b=b_cnt,
             depth=a_cnt+b_cnt, psi_beta=b_cnt/(a_cnt+b_cnt))
}
get_jxn <- function(srp){
  pj <- subset(ap, project==srp & project_type=="data_sources")
  create_rse(pj, type="jxn")
}
cat("S3 SRP154474 junction RSE ...\n"); rj_s3 <- get_jxn("SRP154474")
cat("S4 SRP057087 junction RSE ...\n"); rj_s4 <- get_jxn("SRP057087")
cat("\n== S3 junction match ==\n"); psi_s3 <- psi_one(rj_s3, cls_s3)
cat("== S4 junction match ==\n"); psi_s4 <- psi_one(rj_s4, cls_s4)

summ_psi <- function(d, srp, groups){
  cat(sprintf("\n===== %s =====\n", srp))
  d2 <- d[d$class %in% groups & !is.na(d$class),]
  for(g in groups){
    x <- d2[d2$class==g,]
    cat(sprintf("  %-6s n=%2d | median depth=%.0f | n(depth>=20)=%2d | mean PSI_beta(all)=%.3f  mean PSI_beta(depth>=20)=%.3f\n",
        g, nrow(x), median(x$depth), sum(x$depth>=20),
        mean(x$psi_beta,na.rm=TRUE),
        ifelse(sum(x$depth>=20)>0, mean(x$psi_beta[x$depth>=20]),NA)))
  }
  d2
}
d3 <- summ_psi(psi_s3, "S3 SRP154474 (conv. plaque PP vs control NN)", c("PP","NN"))
d4 <- summ_psi(psi_s4, "S4 SRP057087 (lesional PP vs peri-lesional PN)", c("PP","PN"))

test_shift <- function(d2, gA, gB, srp){
  d2 <- d2[d2$depth>=20,]
  a <- d2$psi_beta[d2$class==gA]; b <- d2$psi_beta[d2$class==gB]
  if(length(a)<3 || length(b)<3){ cat(sprintf("  %s: too few depth>=20 samples (%s=%d, %s=%d) - not testable\n",srp,gA,length(a),gB,length(b))); return(NULL)}
  w <- suppressWarnings(wilcox.test(a,b))
  delta <- mean(a)-mean(b); se <- sqrt(var(a)/length(a)+var(b)/length(b))
  cat(sprintf("  %s: PSI_beta %s - %s = %+.4f pp*100=%.2f  Wilcoxon p=%.3f  (n %s=%d, %s=%d, depth>=20)\n",
      srp, gA, gB, delta, delta*100, w$p.value, gA, length(a), gB, length(b)))
  data.frame(study=srp, contrast=paste0(gA,"-",gB), delta=delta, se=se, p=w$p.value, nA=length(a), nB=length(b))
}
cat("\n=== Exon-exon junction PSI_beta shift tests (depth>=20) ===\n")
t3 <- test_shift(d3,"PP","NN","SRP154474")
t4 <- test_shift(d4,"PP","PN","SRP057087")

psi_beta_all <- tryCatch(readRDS("psi_beta_allstudies.rds"), error=function(e) NULL)
cat("psi_beta_allstudies.rds class:", class(psi_beta_all), "\n")
if(!is.null(psi_beta_all)) print(psi_beta_all)

setDT(psi_beta_all)
eff_one <- function(dt){
  d <- dt[depth>=20 & class %in% c("PP","NN")]
  pp <- d$psi_beta[d$class=="PP"]; nn <- d$psi_beta[d$class=="NN"]
  if(length(pp)<3 || length(nn)<3) return(NULL)
  data.table(delta=mean(pp)-mean(nn), se=sqrt(var(pp)/length(pp)+var(nn)/length(nn)),
             nPP=length(pp), nNN=length(nn))
}
orig <- psi_beta_all[, eff_one(.SD), by=srp]

d3f <- d3[d3$depth>=20,]
pp3 <- d3f$psi_beta[d3f$class=="PP"]; nn3 <- d3f$psi_beta[d3f$class=="NN"]
s3_eff <- data.table(srp="SRP154474", delta=mean(pp3)-mean(nn3),
                     se=sqrt(var(pp3)/length(pp3)+var(nn3)/length(nn3)), nPP=length(pp3), nNN=length(nn3))

old_pool <- rma(yi=delta, sei=se, data=orig, method="DL")
new_dt <- rbind(orig, s3_eff)
new_pool <- rma(yi=delta, sei=se, data=new_dt, method="DL")
new_hk <- rma(yi=delta, sei=se, data=new_dt, method="REML", test="knha")

jxn_tab <- data.frame(
  study=c("SRP035988","SRP165679","SRP126422","SRP154474 (NEW)","SRP057087 (NEW)"),
  contrast=c("PP-NN","PP-NN","PP-NN","PP-NN","PP-PN"),
  median_depth_PP=c(NA,NA,NA,round(median(d3$depth[d3$class=="PP"])),round(median(d4$depth[d4$class=="PP"]))),
  psi_beta_shift_pp=round(c(orig$delta[orig$srp=="SRP035988"], orig$delta[orig$srp=="SRP165679"],
                            orig$delta[orig$srp=="SRP126422"], s3_eff$delta, t4$delta)*100,3),
  wilcox_p=c(0.017,0.70,0.86, round(t3$p,3), round(t4$p,3)))
write.csv(jxn_tab, "stat3_isoform_junction_7study.csv", row.names=FALSE)