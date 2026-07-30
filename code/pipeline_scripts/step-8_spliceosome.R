# 8. Spliceosome Dysregulation in Psoriasis ----
# Purpose: test whether the core spliceosome (and U5 snRNP component EFTUD2
# in particular) is differentially expressed in psoriatic skin, and whether
# its expression covaries with STAT3 and downstream inflammatory programs.
# This sets up the hypothesis tested in Steps 9-10.

# 1. Setup & Load Packages -----
library(tidyverse)
library(limma)
library(edgeR)
library(GSVA)
library(msigdbr)
library(pheatmap)
library(ggplot2)

# Requires: myTopHits.df (Step 4), v.DEGList (Step 4), log2.cpm.filtered.norm
# (Step 3), design and contrast.matrix (Step 4), targets (Step 1).

# 2. Define Core Spliceosome Components -----
# U5 snRNP catalytic core + auxiliary factors implicated in disease splicing
core_factors <- c("EFTUD2", "PRPF8", "SNRNP200", "PRPF4", "PRPF3", "PRPF31",
                  "RAVER1", "PRPF40A", "RBM22")

# 3. Pull Differential Expression Stats for Spliceosome -----
machinery_stats <- myTopHits.df %>%
  dplyr::filter(geneID %in% core_factors) %>%
  dplyr::select(geneID, logFC, adj.P.Val, AveExpr) %>%
  dplyr::arrange(adj.P.Val)

cat("=== Spliceosome DEGs (PP vs NN) ===\n")
print(machinery_stats)

write_tsv(machinery_stats, "spliceosome_DEG_stats.tsv")

# 4. Boxplot: Spliceosome Expression by Group -----
targets_clean <- as.data.frame(targets) %>%
  tibble::rownames_to_column("sample_id")

splice_data <- v.DEGList$E[rownames(v.DEGList$E) %in% core_factors, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "sample_id", values_to = "Expression") %>%
  left_join(targets_clean, by = "sample_id")

fig_spliceosome_box <- ggplot(splice_data, aes(x = Gene, y = Expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1),
              size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c("NN" = "#2C467A", "PP" = "#BE684D")) +
  theme_bw(base_size = 14) +
  labs(title = "Core Spliceosome Components are Upregulated in Psoriasis",
       subtitle = "Psoriasis (PP) vs Normal Skin (NN)",
       y = "Expression (log2 CPM)", x = "") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("fig_spliceosome_boxplot.png", fig_spliceosome_box,
       width = 8, height = 5, dpi = 300)

# 5. Heatmap: Spliceosome Expression Across Samples -----
# 1. Subset the matrix for core factors
heat_mat <- v.DEGList$E[rownames(v.DEGList$E) %in% core_factors, , drop = FALSE]

# 2. Ensure matrix column names match the targets external_id
# If you previously renamed columns (e.g., in Step 5/7), reset them here
colnames(heat_mat) <- targets$external_id 

# 3. Prepare annotation data frame using the same IDs
anno_df <- data.frame(
  group = targets$group,
  row.names = targets$external_id  # Use external_id as row names to match matrix
)

# 4. Alignment Check: Ensure columns of heat_mat and rows of anno_df match
# This prevents the "subscript out of bounds" error
common_ids <- intersect(colnames(heat_mat), rownames(anno_df))
heat_mat <- heat_mat[, common_ids]
anno_df <- anno_df[common_ids, , drop = FALSE]

# 5. Define Colors and Generate Heatmap
ann_colors <- list(group = c(NN = "#2C467A", PP = "#BE684D"))

pheatmap(heat_mat,
         annotation_col    = anno_df,
         annotation_colors = ann_colors,
         scale             = "row",
         show_colnames     = FALSE,
         border_color      = NA,
         color             = colorRampPalette(c("#2C467A", "white", "#BE684D"))(100),
         main              = "Spliceosome Expression Signature",
         filename          = "fig_spliceosome_heatmap.png",
         width             = 8, height = 4)

# 6. Direct Correlation: EFTUD2 vs STAT3 -----
# This is the central observation linking spliceosome upregulation
# to the STAT3 axis tested in Step 9.
expr_matrix_t <- as.data.frame(t(log2.cpm.filtered.norm))

cor_eftud2_stat3 <- cor.test(expr_matrix_t$EFTUD2,
                             expr_matrix_t$STAT3,
                             method = "pearson")

cat("\n=== EFTUD2 ~ STAT3 (Pearson) ===\n")
print(cor_eftud2_stat3)

# Scatter plot
plot_df <- data.frame(
  EFTUD2 = expr_matrix_t$EFTUD2,
  STAT3  = expr_matrix_t$STAT3,
  Group  = targets$group
)

r_val <- round(cor_eftud2_stat3$estimate, 2)
p_val <- format.pval(cor_eftud2_stat3$p.value, digits = 2)

fig_eftud2_stat3 <- ggplot(plot_df, aes(x = EFTUD2, y = STAT3)) +
  geom_point(aes(color = Group), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = "grey80", linewidth = 1) +
  scale_color_manual(values = c("NN" = "#2C467A", "PP" = "#BE684D")) +
  annotate("text",
           x = min(plot_df$EFTUD2),
           y = max(plot_df$STAT3),
           label = paste0("Pearson r = ", r_val, "\np = ", p_val),
           hjust = 0, vjust = 1, size = 5, fontface = "italic") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        axis.line = element_line(linewidth = 0.8)) +
  labs(title = "EFTUD2 Expression Tracks with STAT3",
       subtitle = "Suggests coupling between spliceosome and STAT3 program",
       x = "EFTUD2 Expression (log2 CPM)",
       y = "STAT3 Expression (log2 CPM)")

ggsave("fig_eftud2_stat3_correlation.png", fig_eftud2_stat3,
       width = 7, height = 5, dpi = 300)

# 7. EFTUD2 ~ Inflammatory Pathway Activity (Hallmark GSVA) -----
# Tests whether EFTUD2 expression covaries with broader inflammatory signatures
# beyond STAT3 alone.

hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

GSVA_hallmark_matrix <- gsva(gsvaParam(exprData = v.DEGList$E,
                                       geneSets = hallmark_list,
                                       minSize  = 5,
                                       maxSize  = 500))

eftud2_expr <- log2.cpm.filtered.norm["EFTUD2", ]

target_hallmarks <- c(
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
)

cor_results <- lapply(target_hallmarks, function(pw) {
  if (!pw %in% rownames(GSVA_hallmark_matrix)) return(NULL)
  test <- cor.test(eftud2_expr, GSVA_hallmark_matrix[pw, ])
  data.frame(
    Pathway = pw,
    R       = round(test$estimate, 3),
    P       = format.pval(test$p.value, digits = 2)
  )
}) %>% bind_rows()

cat("\n=== EFTUD2 vs Hallmark Pathway Activity ===\n")
print(cor_results)
write_tsv(cor_results, "eftud2_hallmark_correlations.tsv")

# Save GSVA matrix for reuse in Step 11
save(GSVA_hallmark_matrix, eftud2_expr, file = "gsva_hallmark.RData")

cat("\nStep 8 complete. Outputs:\n")
cat("  - spliceosome_DEG_stats.tsv\n")
cat("  - fig_spliceosome_boxplot.png\n")
cat("  - fig_spliceosome_heatmap.png\n")
cat("  - fig_eftud2_stat3_correlation.png\n")
cat("  - eftud2_hallmark_correlations.tsv\n")
cat("  - gsva_hallmark.RData (for Step 11)\n")
