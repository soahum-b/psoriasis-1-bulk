# Auto-extracted generating script
# Produces: pathway_timing_stats.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: e1704b85-4740-460f-a2e7-f33362ea1f9a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(msigdbr)})

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
fwrite(res[order(-pn_frac)], "pathway_timing_stats.csv")