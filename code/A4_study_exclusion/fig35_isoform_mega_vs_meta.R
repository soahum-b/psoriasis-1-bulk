# Auto-extracted generating script
# Produces: fig35_isoform_mega_vs_meta.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: b63969ea-970f-40ba-b28c-76367886940a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(ggplot2); library(patchwork); library(data.table)})
theme_set(theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=12), plot.title.position="plot"))
BLUE<-"#4C72B0"; RED<-"#C44E52"; GREY<-"#8C8C8C"; GOLD<-"#DD8452"

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
t3 <- test_shift(d3,"PP","NN","SRP154474")
t4 <- test_shift(d4,"PP","PN","SRP057087")

psi_beta_all <- tryCatch(readRDS("psi_beta_allstudies.rds"), error=function(e) NULL)

suppressMessages({library(metafor)})
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

setDT(psi_beta_all)
comb <- psi_beta_all[class %in% c("PP","NN"), .(srp, sample, class, a, b, depth)]
s3add <- as.data.table(d3)[class %in% c("PP","NN"), .(srp="SRP154474", sample, class, a, b, depth)]
comb <- rbind(comb, s3add)

comb <- comb[depth>=20]
has_both <- comb[, .(nc=uniqueN(class)), by=srp][nc==2, srp]
comb <- comb[srp %in% has_both]
comb[, class := factor(class, levels=c("NN","PP"))]
comb[, srp := factor(srp)]

m_mega  <- glm(cbind(b, a) ~ class + srp, family=quasibinomial, data=comb)

co_m <- summary(m_mega)$coef["classPP",]

studies <- levels(comb$srp)
pred_grid <- expand.grid(class=factor(c("NN","PP"),levels=c("NN","PP")), srp=factor(studies,levels=studies))
pred_grid$p <- predict(m_mega, pred_grid, type="response")
mm <- aggregate(p~class, pred_grid, mean)
pp_diff <- (mm$p[mm$class=="PP"]-mm$p[mm$class=="NN"])*100

comb[, psi := b/(a+b)]
m_lm <- lm(psi ~ class + srp, data=comb, weights=depth)
co_lm <- summary(m_lm)$coef["classPP",]

cmpf <- data.frame(
  method=c("Meta: DL random-effects","Meta: HKSJ","Mega: quasibinomial GLM\n(+ study)","Mega: depth-weighted LM\n(+ study)"),
  type=c("meta","meta","mega","mega"),
  est=c(new_pool$b*100, new_hk$b*100, pp_diff, co_lm[1]*100),
  lo=c(new_pool$ci.lb*100, new_hk$ci.lb*100,
       (mm$p[2]-mm$p[1])*100 - 1.96*summary(m_mega)$coef["classPP",2]*100*0.25,
       co_lm[1]*100-1.96*co_lm[2]*100),
  p=c(new_pool$pval, new_hk$pval, co_m[4], co_lm[4]))
cmpf$lo[3] <- NA; cmpf$hi <- c(new_pool$ci.ub*100, new_hk$ci.ub*100, NA, co_lm[1]*100+1.96*co_lm[2]*100)
cmpf$sig <- cmpf$p<0.05
cmpf$method <- factor(cmpf$method, levels=rev(cmpf$method))
cmpf$lab <- sprintf("%+.2f pp, p=%.3f", cmpf$est, cmpf$p)

k_pp_per_logit <- pp_diff / co_m[1]
glm_lo <- (co_m[1]-1.96*co_m[2]) * k_pp_per_logit
glm_hi <- (co_m[1]+1.96*co_m[2]) * k_pp_per_logit

cmpf$lo[3] <- glm_lo; cmpf$hi[3] <- glm_hi
pC <- ggplot(cmpf, aes(est, method)) +
  geom_vline(xintercept=0, linetype=2, color=GREY) +
  geom_errorbarh(aes(xmin=lo, xmax=hi, color=type), height=0.15, linewidth=0.8) +
  geom_point(aes(color=type, shape=sig), size=4) +
  geom_text(aes(label=lab), vjust=-1.1, size=3.1) +
  scale_color_manual(values=c(meta=RED, mega="#4C72B0"), guide="none") +
  scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1), name="p<0.05") +
  labs(x="STAT3-beta PSI shift, PP - NN (percentage points)", y=NULL,
       title="Combining 268 samples directly (mega) does not confirm the isoform shift",
       subtitle="Only the per-study random-effects mean clears p<0.05; sample-level pooling (higher power) says n.s.") +
  coord_cartesian(xlim=c(-1.2,2.6)) +
  theme(legend.position=c(0.88,0.15), legend.background=element_blank(),
        axis.text.y=element_text(size=8.5))
ggsave("fig35_isoform_mega_vs_meta.png", pC, width=10, height=4.6, dpi=150)