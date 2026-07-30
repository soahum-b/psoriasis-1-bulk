# Auto-extracted generating script
# Produces: ssgsea_panel.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_results_all.rds, dge_filt_norm.rds
# Source artifact version: 213f77e5-2c69-4a51-9411-2fbd3d1fe630
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(edgeR)
  library(limma)
  library(msigdbr)
})

dge <- readRDS("dge_filt_norm.rds")

grp <- factor(dge$samples$group, levels = c("NN", "PP"))
design <- model.matrix(~grp)
colnames(design) <- c("Intercept", "PPvsNN")

v <- voom(dge, design, plot = FALSE)

h_df <- msigdbr(species = "Homo sapiens", category = "H")
H <- split(h_df$gene_symbol, h_df$gs_name)

ssgsea_scores <- function(expr, gene_sets, alpha = 0.25, norm = TRUE, min_size = 5) {
  genes <- rownames(expr)
  n <- nrow(expr)
  ns <- ncol(expr)
  R <- apply(expr, 2, function(x) rank(x, ties.method = "average"))
  gene_sets <- lapply(gene_sets, function(g) intersect(g, genes))
  gene_sets <- gene_sets[sapply(gene_sets, length) >= min_size]
  ES <- matrix(NA, length(gene_sets), ns,
               dimnames = list(names(gene_sets), colnames(expr)))
  for (j in seq_len(ns)) {
    ord <- order(R[, j], decreasing = TRUE)
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

panel_sets <- H[c("HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MTORC1_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYOGENESIS", "HALLMARK_ADIPOGENESIS", "HALLMARK_FATTY_ACID_METABOLISM")]
sc <- ssgsea_scores(v$E, panel_sets)
rownames(sc) <- gsub("HALLMARK_", "", rownames(sc))

saveRDS(list(panel_scores = sc), "ssgsea_panel.rds")