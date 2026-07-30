# Auto-extracted generating script
# Produces: ssgsea_scores.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): stat3_isoform_psi.rds, gsea_results_all.rds, de_results_full.rds, dge_filt_norm.rds
# Source artifact version: f973b1cf-90e8-4199-9e68-83f1440fd368
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR); library(limma); library(msigdbr)})

dge <- readRDS("dge_filt_norm.rds")
de_tab <- readRDS("de_results_full.rds")
psi <- readRDS("stat3_isoform_psi.rds")

grp <- factor(dge$samples$group, levels=c("NN","PP"))
design <- model.matrix(~grp)
colnames(design) <- c("Intercept","PPvsNN")
v <- voom(dge, design, plot=FALSE)

h_df <- msigdbr(species="Homo sapiens", category="H")
H <- split(h_df$gene_symbol, h_df$gs_name)

ssgsea_scores <- function(expr, gene_sets, alpha=0.25, norm=TRUE, min_size=5) {
  genes <- rownames(expr); n <- nrow(expr); ns <- ncol(expr)
  R <- apply(expr, 2, function(x) rank(x, ties.method="average"))
  gene_sets <- lapply(gene_sets, function(g) intersect(g, genes))
  gene_sets <- gene_sets[sapply(gene_sets, length) >= min_size]
  ES <- matrix(NA, length(gene_sets), ns,
               dimnames=list(names(gene_sets), colnames(expr)))
  for (j in seq_len(ns)) {
    ord <- order(R[,j], decreasing=TRUE)
    rnk <- R[ord, j]
    g_ord <- genes[ord]
    for (s in seq_along(gene_sets)) {
      inset <- g_ord %in% gene_sets[[s]]
      w <- abs(rnk)^alpha
      Pin  <- cumsum(ifelse(inset, w, 0)); Pin  <- Pin / Pin[n]
      Pout <- cumsum(ifelse(inset, 0, 1)); Pout <- Pout / Pout[n]
      ES[s, j] <- sum(Pin - Pout)
    }
  }
  if (norm) ES <- ES / (max(ES) - min(ES))
  ES
}

sets_use <- c(H["HALLMARK_IL6_JAK_STAT3_SIGNALING"],
              H["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
              H["HALLMARK_INFLAMMATORY_RESPONSE"],
              H["HALLMARK_E2F_TARGETS"],
              H["HALLMARK_TNFA_SIGNALING_VIA_NFKB"],
              H["HALLMARK_MYC_TARGETS_V1"],
              H["HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
              H["HALLMARK_MYOGENESIS"])
scores <- ssgsea_scores(v$E, sets_use)

logcpm <- cpm(dge, log=TRUE)
stat3_expr <- logcpm["STAT3", ]
il6 <- scores["HALLMARK_IL6_JAK_STAT3_SIGNALING", ]

sid <- rownames(dge$samples)
df <- data.frame(sample=sid, il6=il6, stat3=stat3_expr, grp=grp)
m <- merge(df, as.data.frame(psi)[,c("sample","psi_beta")], by="sample")

c_expr <- cor.test(m$il6, m$stat3, method="spearman")
c_psi  <- cor.test(m$il6, m$psi_beta, method="spearman")
c_psi_pp <- cor.test(m$il6[m$grp=="PP"], m$psi_beta[m$grp=="PP"], method="spearman")
c_psi_nn <- cor.test(m$il6[m$grp=="NN"], m$psi_beta[m$grp=="NN"], method="spearman")

saveRDS(list(scores=scores, merged=m, cor_expr=c_expr, cor_psi=c_psi,
             cor_psi_pp=c_psi_pp, cor_psi_nn=c_psi_nn), "ssgsea_scores.rds")