# Auto-extracted generating script
# Produces: ssgsea_3group.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): sample_classification.rds
# Source artifact version: 39d80619-4f97-4e84-95da-804bc824851a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(msigdbr)

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

cls <- as.data.table(readRDS("sample_classification.rds"))

mat <- fread("clust_input/SRP165679.tsv"); X<-as.matrix(mat[,-1]); rownames(X)<-mat$Genes
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
saveRDS(list(scores=sc, grp=grp), "ssgsea_3group.rds")