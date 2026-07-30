# Auto-extracted generating script
# Produces: meta_summary_stats.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds
# Source artifact version: 7b8d327d-cae3-46d0-8a79-6dbfac547125
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)

de <- readRDS("per_study_de.rds")

build_mats <- function(delist) {
  studies <- names(delist)
  genes <- sort(unique(unlist(lapply(delist, function(x) x$gene))))
  Y <- matrix(NA, length(genes), length(studies), dimnames=list(genes, studies))
  V <- Y
  for (s in studies) {
    x <- delist[[s]]
    idx <- match(x$gene, genes)
    Y[idx, s] <- x$logFC
    V[idx, s] <- x$SE^2
  }
  list(Y=Y, V=V, studies=studies, genes=genes)
}

dl_meta <- function(Y, V) {
  W <- 1/V
  W[is.na(Y)] <- NA
  k <- rowSums(!is.na(Y))
  sumW  <- rowSums(W, na.rm=TRUE)
  sumWY <- rowSums(W*Y, na.rm=TRUE)
  mu_fe <- sumWY / sumW
  Q <- rowSums(W*(Y-mu_fe)^2, na.rm=TRUE)
  sumW2 <- rowSums(W^2, na.rm=TRUE)
  C <- sumW - sumW2/sumW
  tau2 <- pmax(0, (Q-(k-1))/C)
  tau2[k<2] <- NA
  Wr <- 1/(V + matrix(tau2, nrow(Y), ncol(Y)))
  Wr[is.na(Y)] <- NA
  sumWr  <- rowSums(Wr, na.rm=TRUE)
  mu_re  <- rowSums(Wr*Y, na.rm=TRUE)/sumWr
  se_re  <- sqrt(1/sumWr)
  z <- mu_re/se_re
  p <- 2*pnorm(-abs(z))
  I2 <- ifelse(Q>0 & k>1, pmax(0, 100*(Q-(k-1))/Q), 0)
  data.table(gene=rownames(Y), k=k, logFC=mu_re, SE=se_re, z=z, P=p,
             tau2=tau2, I2=I2, Q=Q)
}

meta_res <- list()
for (cn in names(de)) {
  m <- build_mats(de[[cn]])
  res <- dl_meta(m$Y, m$V)
  res <- res[k>=2]
  res[, FDR := p.adjust(P, "BH")]
  res <- res[order(P)]
  meta_res[[cn]] <- res
}

a <- as.data.table(de$PPvsNN$SRP035988)
single_strict <- a[adjP<0.05 & abs(logFC)>1, gene]
mp <- meta_res$PPvsNN
meta_strict <- mp[FDR<0.05 & abs(logFC)>1, gene]
gained <- setdiff(meta_strict, single_strict)

summary_stats <- list(
  single_strict_n = length(single_strict),
  meta_strict_n   = length(meta_strict),
  gained_n        = length(gained),
  gained_genes    = gained,
  shared_n        = length(intersect(meta_strict, single_strict)),
  meta_fdr_only_n = nrow(mp[FDR<0.05]),
  stat3           = as.list(mp[gene=="STAT3", .(k,logFC,SE,z,P,FDR,I2,tau2)])
)

saveRDS(summary_stats, "meta_summary_stats.rds")