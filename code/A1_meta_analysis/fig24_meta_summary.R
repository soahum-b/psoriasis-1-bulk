# Auto-extracted generating script
# Produces: fig24_meta_summary.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds
# Source artifact version: e7b3b668-3956-4b6d-a5ad-88246338ebd8
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(ggplot2)
library(patchwork)
library(ggrepel)

de <- readRDS("per_study_de.rds")

NNblue <- "#4C72B0"; PPred <- "#C44E52"; grey <- "#B0B0B0"; green <- "#55A868"

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

mp <- copy(meta_res$PPvsNN)
a <- as.data.table(de$PPvsNN$SRP035988)
single_strict <- a[adjP<0.05 & abs(logFC)>1, gene]
meta_strict <- mp[FDR<0.05 & abs(logFC)>1, gene]
gained <- setdiff(meta_strict, single_strict)

mp[, cat := "ns"]
mp[FDR<0.05 & abs(logFC)>1, cat := "meta-sig"]
mp[gene %in% gained, cat := "gained"]
mp[, nlp := pmin(-log10(P), 60)]
hi <- c("STAT3","SOCS3","STAT1","S100A7","DEFB4A","CCR1","C1QB","PTPN7")
mp[, lab := ifelse(gene %in% hi, gene, "")]

pA <- ggplot(mp, aes(logFC, nlp)) +
  geom_point(data=mp[cat=="ns"], color=grey, size=0.4, alpha=0.3) +
  geom_point(data=mp[cat=="meta-sig"], color=PPred, size=0.5, alpha=0.5) +
  geom_point(data=mp[cat=="gained"], color=green, size=0.9, alpha=0.9) +
  geom_point(data=mp[gene %in% hi], color="black", size=1.5) +
  geom_text_repel(aes(label=lab), size=2.6, max.overlaps=20, min.segment.length=0) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color=grey, linewidth=0.3) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color=grey, linewidth=0.3) +
  labs(x="pooled log2 fold-change (random-effects)", y="-log10 P (capped 60)",
       title="A. Meta-analysis PP vs NN",
       subtitle=sprintf("green = %d genes gained over single study (cross-study replicated)", length(gained))) +
  theme_bw(base_size=9) + theme(plot.subtitle=element_text(size=7,color=green))

sigg <- mp[FDR<0.05 & abs(logFC)>1]
lowI <- sum(sigg$I2<25); hiI <- sum(sigg$I2>=75)
pB <- ggplot(sigg, aes(I2)) +
  geom_histogram(binwidth=5, fill=NNblue, color="white", linewidth=0.2) +
  geom_vline(xintercept=median(sigg$I2), linetype="dashed", color=PPred) +
  annotate("text", x=median(sigg$I2)-3, y=Inf, vjust=2, hjust=1,
           label=sprintf("median I2 = %.0f%%", median(sigg$I2)), size=2.8, color=PPred) +
  labs(x="I-squared (between-study heterogeneity, %)", y="genes",
       title="B. Heterogeneity among meta-significant genes",
       subtitle=sprintf("bimodal: %d genes I2<25%% (highly consistent), %d genes I2>75%% (effect-size varies across cohorts)", lowI, hiI)) +
  theme_bw(base_size=9) + theme(plot.subtitle=element_text(size=6.6))

p <- pA + pB + plot_layout(widths=c(1.4,1))
ggsave("fig24_meta_summary.png", p, width=11, height=4.4, dpi=150)