# Auto-extracted generating script
# Produces: fig_s3_pathway_timing.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): sample_classification.rds
# Source artifact version: 8f34ddb2-907f-4f69-a92a-f07c9599da1b
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(msigdbr)
library(ggplot2)

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")

# base-R ssGSEA (Barbie 2009)
ssgsea_scores <- function(X, gsets, alpha=0.25){
  R <- apply(X, 2, function(x) rank(x, ties.method="average"))
  N <- nrow(X); scores <- matrix(NA, length(gsets), ncol(X),
        dimnames=list(names(gsets), colnames(X)))
  for(s in names(gsets)){
    inset <- rownames(X) %in% gsets[[s]]
    for(j in seq_len(ncol(X))){
      ord <- order(R[,j], decreasing=TRUE); ind <- inset[ord]
      rw <- abs(R[ord,j])^alpha
      Pin <- cumsum(ifelse(ind, rw, 0)); Pin <- Pin/Pin[N]
      Pout <- cumsum(ifelse(!ind, 1, 0)); Pout <- Pout/Pout[N]
      scores[s,j] <- sum(Pin - Pout)
    }
  }
  # range-normalize per set
  t(apply(scores,1,function(z)(z-min(z))/(max(z)-min(z))))
}

mat <- fread("clust_input/SRP165679.tsv"); X<-as.matrix(mat[,-1]); rownames(X)<-mat$Genes
cls <- readRDS("sample_classification.rds"); cls<-as.data.table(cls)
lab <- cls[srp=="SRP165679"][class %in% c("NN","PN","PP")]
X <- X[, lab$external_id]; grp <- factor(lab$class, levels=c("NN","PN","PP"))

H <- as.data.table(msigdbr(species="Homo sapiens", category="H"))
keysets <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS","HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MTORC1_SIGNALING","HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_ADIPOGENESIS")
gsets <- lapply(keysets, function(s) H[gs_name==s, gene_symbol]); names(gsets)<-sub("HALLMARK_","",keysets)

sc <- ssgsea_scores(X, gsets)

stage <- as.integer(grp)-1L

res <- rbindlist(lapply(rownames(sc), function(p){
  y <- sc[p,]
  lp <- summary(lm(y~stage))$coefficients["stage",4]
  pn_p <- wilcox.test(y[grp=="PN"], y[grp=="NN"])$p.value
  data.table(pathway=p, trend_p=lp, pn_vs_nn_p=pn_p,
             NN=mean(y[grp=="NN"]), PN=mean(y[grp=="PN"]), PP=mean(y[grp=="PP"]))
}))
res[, pn_frac := (PN-NN)/(PP-NN)]
res[, trend_FDR := p.adjust(trend_p,"BH")]
res[, pn_FDR := p.adjust(pn_vs_nn_p,"BH")]

plotd <- rbindlist(lapply(rownames(sc), function(p){
  data.table(pathway=p, stage=factor(c("NN","PN","PP"),levels=c("NN","PN","PP")),
    mean=tapply(sc[p,],grp,mean),
    se=tapply(sc[p,],grp,function(z)sd(z)/sqrt(length(z))))
}))
prolif <- c("E2F_TARGETS","G2M_CHECKPOINT","MYC_TARGETS_V1","MTORC1_SIGNALING")
immune <- c("INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","IL6_JAK_STAT3_SIGNALING",
            "INFLAMMATORY_RESPONSE","TNFA_SIGNALING_VIA_NFKB","ALLOGRAFT_REJECTION")
plotd <- plotd[pathway %in% c(prolif,immune)]
plotd[, program := ifelse(pathway %in% immune, "Immune / interferon (early)", "Proliferation (late)")]
plotd[, lbl := gsub("_"," ", pathway)]

th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold",size=13), strip.text=element_text(face="bold"))
p <- ggplot(plotd, aes(stage, mean, group=lbl, colour=program))+
  geom_line(linewidth=.8, alpha=.8)+geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width=.12, alpha=.6)+
  facet_wrap(~program)+
  scale_colour_manual(values=c("Immune / interferon (early)"="#55A868","Proliferation (late)"="#8172B3"))+
  labs(title="Pathway activation timing along the staging axis (SRP165679)",
       subtitle="Immune/interferon programs already rise at peri-lesional (PN); proliferation stays flat until lesional (PP)",
       x="Staging axis", y="ssGSEA score (mean +/- SE)")+
  th+theme(legend.position="none", plot.subtitle=element_text(size=9.5))
ggsave("fig_s3_pathway_timing.png", p, width=10, height=4.8, dpi=150)