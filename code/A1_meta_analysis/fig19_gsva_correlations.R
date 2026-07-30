# Auto-extracted generating script
# Produces: fig19_gsva_correlations.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): stat3_isoform_psi.rds, gsea_results_all.rds, de_results_full.rds, dge_filt_norm.rds
# Source artifact version: df12b72c-0eac-4fc6-850b-34baf5a37249
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(edgeR)
  library(limma)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
})

NNblue <- "#4C72B0"
PPred <- "#C44E52"

dge <- readRDS("dge_filt_norm.rds")
de_tab <- readRDS("de_results_full.rds")

grp <- factor(dge$samples$group, levels = c("NN", "PP"))
design <- model.matrix(~ grp)
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
      Pin <- cumsum(ifelse(inset, w, 0)); Pin <- Pin / Pin[n]
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

logcpm <- cpm(dge, log = TRUE)
stat3_expr <- logcpm["STAT3", ]

psi <- readRDS("stat3_isoform_psi.rds")

il6 <- scores["HALLMARK_IL6_JAK_STAT3_SIGNALING", ]
sid <- rownames(dge$samples)

df <- data.frame(sample = sid, il6 = il6, stat3 = stat3_expr, grp = grp)
m <- merge(df, as.data.frame(psi)[, c("sample", "psi_beta")], by = "sample")

c_expr <- cor.test(m$il6, m$stat3, method = "spearman")
c_psi <- cor.test(m$il6, m$psi_beta, method = "spearman")
c_psi_pp <- cor.test(m$il6[m$grp == "PP"], m$psi_beta[m$grp == "PP"], method = "spearman")
c_psi_nn <- cor.test(m$il6[m$grp == "NN"], m$psi_beta[m$grp == "NN"], method = "spearman")

pA <- ggplot(m, aes(stat3, il6, color = grp)) +
  geom_point(alpha = 0.8, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  scale_color_manual(values = c(NN = NNblue, PP = PPred), name = NULL) +
  labs(x = "STAT3 expression (log2-CPM)", y = "IL6/JAK/STAT3 ssGSEA score",
       title = "A. Per-sample pathway score tracks STAT3 expression",
       subtitle = sprintf("Spearman rho=%.2f, p=%.1e", c_expr$estimate, c_expr$p.value)) +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 8))

pB <- ggplot(m, aes(psi_beta * 100, il6, color = grp)) +
  geom_point(alpha = 0.8, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, aes(group = grp), linewidth = 0.6) +
  scale_color_manual(values = c(NN = NNblue, PP = PPred), name = NULL) +
  labs(x = "PSI_beta (% STAT3-beta)", y = "IL6/JAK/STAT3 ssGSEA score",
       title = "B. Pathway score vs beta-isoform fraction",
       subtitle = sprintf("all rho=%.2f (p=%.2f); within PP rho=%.2f (p=%.3f), NN rho=%.2f (p=%.3f)",
         c_psi$estimate, c_psi$p.value, c_psi_pp$estimate, c_psi_pp$p.value,
         c_psi_nn$estimate, c_psi_nn$p.value)) +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 7))

fig19 <- pA / pB
ggsave("fig19_gsva_correlations.png", fig19, width = 8, height = 8, dpi = 150)