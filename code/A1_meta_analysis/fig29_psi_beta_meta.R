# Auto-extracted generating script
# Produces: fig29_psi_beta_meta.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 972457b5-8222-4141-b3c6-6c3163a066ba
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(data.table)
  library(metafor)
  library(ggplot2)
  library(patchwork)
})

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/")

# Sample classification
cls_all <- readRDS("/tmp/cls_all.rds")
setDT(cls_all)

ap <- available_projects(organism="human")

acc <- 42317181L; a_don <- 42316902L; b_don <- 42316852L; tol <- 2L

psi_for_study <- function(srp) {
  proj <- subset(ap, project==srp & project_type=="data_sources")
  rj <- create_rse(proj, type="jxn")
  gr <- rowRanges(rj)
  onchr <- as.character(seqnames(gr)) %in% c("chr17","17")
  is_a <- onchr & abs(start(gr)-a_don)<=tol & abs(end(gr)-acc)<=tol
  is_b <- onchr & abs(start(gr)-b_don)<=tol & abs(end(gr)-acc)<=tol
  if (sum(is_a)<1 || sum(is_b)<1) return(data.table(srp=srp, note="junction not found",
                                                    n_a=sum(is_a), n_b=sum(is_b)))
  a_cnt <- colSums(matrix(assay(rj,"counts")[is_a,,drop=FALSE], ncol=ncol(rj)))
  b_cnt <- colSums(matrix(assay(rj,"counts")[is_b,,drop=FALSE], ncol=ncol(rj)))
  cmap <- cls_all[cls_all$srp==srp]
  cl <- cmap$class[match(colnames(rj), cmap$external_id)]
  dt <- data.table(srp=srp, sample=colnames(rj), class=cl,
                   a=a_cnt, b=b_cnt, depth=a_cnt+b_cnt)
  dt[, psi_beta := b/depth]
  dt[class %in% c("PP","PN","NN")]
}

studies <- c("SRP035988","SRP165679","SRP076982","SRP126422","SRP016583")
psi_list <- list()
for (s in studies) { message("jxn ", s); psi_list[[s]] <- tryCatch(psi_for_study(s), error=function(e) data.table(srp=s, note=conditionMessage(e))) }
psi_all <- rbindlist(psi_list, fill=TRUE)

MIN_DEPTH <- 20
psi_f <- psi_all[!is.na(psi_beta) & depth>=MIN_DEPTH]

pp_nn_effect <- function(d) {
  pp <- d[class=="PP", psi_beta]; nn <- d[class=="NN", psi_beta]
  if (length(pp)<3 || length(nn)<3) return(NULL)
  m <- mean(pp)-mean(nn)
  se <- sqrt(var(pp)/length(pp) + var(nn)/length(nn))
  data.table(delta=m, se=se, n_pp=length(pp), n_nn=length(nn),
             mean_pp=mean(pp), mean_nn=mean(nn),
             w_p=wilcox.test(pp,nn)$p.value)
}
eff <- psi_f[, pp_nn_effect(.SD), by=srp]

res <- rma(yi=delta, sei=se, data=eff, method="DL")

NNblue<-"#4C72B0"; PPred<-"#C44E52"; grey<-"#B0B0B0"

eff2 <- copy(eff); eff2[, `:=`(d=delta*100, lo=(delta-1.96*se)*100, hi=(delta+1.96*se)*100)]
eff2[, lab := sprintf("%s (%d PP vs %d NN)", srp, n_pp, n_nn)]
pooled_row <- data.table(lab="RE pooled", d=as.numeric(res$b)*100,
                         lo=as.numeric(res$ci.lb)*100, hi=as.numeric(res$ci.ub)*100, srp="POOL")
fa <- rbind(eff2[,.(lab,d,lo,hi,srp)], pooled_row)
fa[, lab := factor(lab, levels=rev(lab))]
fa[, ispool := srp=="POOL"]
pA <- ggplot(fa, aes(d, lab)) +
  geom_vline(xintercept=0, linetype="dashed", color=grey) +
  geom_errorbarh(aes(xmin=lo, xmax=hi, color=ispool), height=0.2) +
  geom_point(aes(color=ispool, shape=ispool, size=ispool)) +
  scale_color_manual(values=c("FALSE"="black","TRUE"=PPred), guide="none") +
  scale_shape_manual(values=c("FALSE"=15,"TRUE"=18), guide="none") +
  scale_size_manual(values=c("FALSE"=2.4,"TRUE"=4), guide="none") +
  labs(x="PSI_beta shift, lesional - healthy (percentage points)", y=NULL,
       title="A. STAT3 beta-inclusion shift does NOT replicate",
       subtitle=sprintf("pooled +%.2f pp [%.2f, %.2f], I2=%.0f%%, p=%.3f (n.s.)",
                        as.numeric(res$b)*100,as.numeric(res$ci.lb)*100,as.numeric(res$ci.ub)*100,res$I2,res$pval)) +
  theme_bw(base_size=9) + theme(plot.subtitle=element_text(size=7))
pf <- psi_all[!is.na(psi_beta) & depth>=20 & class %in% c("PP","NN")]
pf[, cls := factor(class, levels=c("NN","PP"))]
pB <- ggplot(pf, aes(srp, psi_beta*100, color=cls)) +
  geom_boxplot(outlier.size=0.4, width=0.6, position=position_dodge(0.7)) +
  scale_color_manual(values=c("NN"=NNblue,"PP"=PPred), name=NULL, labels=c("healthy","lesional")) +
  labs(x=NULL, y="PSI_beta (%)", title="B. Per-study distributions (depth>=20)",
       subtitle="anchor shows PP>NN; deep replication (SRP165679) shows no gap") +
  theme_bw(base_size=9) +
  theme(plot.subtitle=element_text(size=7), axis.text.x=element_text(angle=30,hjust=1,size=6.5),
        legend.position=c(0.85,0.9), legend.background=element_blank())
p <- pA + pB + plot_layout(widths=c(1,1.1))
ggsave("fig29_psi_beta_meta.png", p, width=11, height=4.2, dpi=150)