# Auto-extracted generating script
# Produces: stat3_isoform_pool_5v7.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 509f324d-ee44-4c76-af12-c4488e1e2cd6
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
cat("=== Original per-study PP-vs-NN PSI_beta effects (depth>=20, >=3/arm) ===\n")
print(orig)

d3f <- d3[d3$depth>=20,]
pp3 <- d3f$psi_beta[d3f$class=="PP"]; nn3 <- d3f$psi_beta[d3f$class=="NN"]
s3_eff <- data.table(srp="SRP154474", delta=mean(pp3)-mean(nn3),
                     se=sqrt(var(pp3)/length(pp3)+var(nn3)/length(nn3)), nPP=length(pp3), nNN=length(nn3))
cat("\nS3 SRP154474 PP-vs-NN PSI_beta effect:\n"); print(s3_eff)

old_pool <- rma(yi=delta, sei=se, data=orig, method="DL")
new_dt <- rbind(orig, s3_eff)
new_pool <- rma(yi=delta, sei=se, data=new_dt, method="DL")
new_hk <- rma(yi=delta, sei=se, data=new_dt, method="REML", test="knha")

cat(sprintf("\nOLD pool (%d studies): shift %+.4f [%.4f, %.4f] I2=%.0f%% p=%.3f -> %s\n",
    nrow(orig), old_pool$b, old_pool$ci.lb, old_pool$ci.ub, old_pool$I2, old_pool$pval,
    ifelse(old_pool$pval<0.05,"replicates","does NOT replicate")))
cat(sprintf("NEW pool (%d studies, +S3): shift %+.4f [%.4f, %.4f] I2=%.0f%% p=%.3f -> %s\n",
    nrow(new_dt), new_pool$b, new_pool$ci.lb, new_pool$ci.ub, new_pool$I2, new_pool$pval,
    ifelse(new_pool$pval<0.05,"replicates","does NOT replicate")))
cat(sprintf("NEW pool HKSJ: shift %+.4f [%.4f, %.4f] p=%.3f\n", new_hk$b, new_hk$ci.lb, new_hk$ci.ub, new_hk$pval))

pool_tab <- data.frame(
  pool=c("3-study (original)","4-study (+S3)"),
  DL=c(sprintf("%+.3f pp [%.3f, %.3f] p=%.3f", old_pool$b*100, old_pool$ci.lb*100, old_pool$ci.ub*100, old_pool$pval),
       sprintf("%+.3f pp [%.3f, %.3f] p=%.3f", new_pool$b*100, new_pool$ci.lb*100, new_pool$ci.ub*100, new_pool$pval)),
  HKSJ=c("+0.79 pp p~0.11 (reported in deep-dive)",
         sprintf("%+.3f pp [%.3f, %.3f] p=%.3f", new_hk$b*100, new_hk$ci.lb*100, new_hk$ci.ub*100, new_hk$pval)),
  verdict=c("does NOT replicate","DL crosses 0.05; HKSJ still n.s."))
write.csv(pool_tab, "stat3_isoform_pool_5v7.csv", row.names=FALSE)
cat("tables saved\n"); print(pool_tab)