# 11. NF-kB and IL-17 Pathway Engagement ----
# Purpose: situate the spliceosome-STAT3 finding within the canonical
# psoriasis inflammatory framework. Tests whether NF-kB and IL-17 signaling
# are enriched in psoriasis (expected) and whether EFTUD2 expression
# covaries with these pathways (mechanistically suggestive).

# 1. Setup & Load Packages -----
library(tidyverse)
library(limma)
library(GSVA)
library(msigdbr)

# Requires: v.DEGList (Step 4), design and contrast.matrix (Step 4),
# log2.cpm.filtered.norm (Step 3), GSVA_hallmark_matrix and eftud2_expr (Step 8).
# If you've restarted R since Step 8, load:
# load("gsva_hallmark.RData")

# 2. Define Targeted Pathway Set (IL-17 + NF-kB) -----
all_pathways <- msigdbr(species = "Homo sapiens")

target_set_names <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_INTERLEUKIN_17_SIGNALING",
  "KEGG_IL_17_SIGNALING_PATHWAY",
  "GOBP_IKK_NFKB_SIGNALING"
)

target_sets <- all_pathways %>%
  dplyr::filter(gs_name %in% target_set_names)

cat("=== Pathways found in MSigDB ===\n")
print(unique(target_sets$gs_name))

target_list <- split(target_sets$gene_symbol, target_sets$gs_name)

# 3. Run GSVA on the Targeted Set -----
GSVA_inflam_matrix <- gsva(gsvaParam(exprData = v.DEGList$E,
                                     geneSets = target_list,
                                     minSize  = 5,
                                     maxSize  = 500))

# 4. Test Differential Pathway Activity (PP vs NN) -----
fit_inflam          <- lmFit(GSVA_inflam_matrix, design)
fit_inflam_contrast <- contrasts.fit(fit_inflam, contrast.matrix)
ebFit_inflam        <- eBayes(fit_inflam_contrast)

inflam_enrichment <- topTable(ebFit_inflam, coef = 1, number = Inf) %>%
  tibble::rownames_to_column("Pathway")

cat("\n=== IL-17 / NF-kB enrichment in Psoriasis ===\n")
print(inflam_enrichment)
write_tsv(inflam_enrichment, "il17_nfkb_enrichment.tsv")

# 5. Correlate EFTUD2 with Each Pathway Score -----
# Tests whether spliceosome upregulation tracks with downstream inflammation.
cor_inflam <- lapply(rownames(GSVA_inflam_matrix), function(pw) {
  test <- cor.test(eftud2_expr, GSVA_inflam_matrix[pw, ])
  data.frame(
    Pathway = pw,
    R       = round(test$estimate, 3),
    P       = format.pval(test$p.value, digits = 2)
  )
}) %>% bind_rows()

cat("\n=== EFTUD2 vs IL-17 / NF-kB pathway activity ===\n")
print(cor_inflam)
write_tsv(cor_inflam, "eftud2_il17_nfkb_correlations.tsv")

# 6. Visualize: Pathway Activity by Group -----
inflam_long <- as.data.frame(GSVA_inflam_matrix) %>%
  tibble::rownames_to_column("Pathway") %>%
  pivot_longer(-Pathway, names_to = "sample_id", values_to = "Score") %>%
  left_join(
    targets %>% as.data.frame() %>% tibble::rownames_to_column("sample_id"),
    by = "sample_id"
  )

fig_inflam_box <- ggplot(inflam_long, aes(x = Pathway, y = Score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15),
              size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c("NN" = "#2C467A", "PP" = "#BE684D")) +
  theme_bw(base_size = 12) +
  labs(title = "IL-17 and NF-kB Activity in Psoriasis",
       y     = "GSVA Enrichment Score",
       x     = "") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("fig_il17_nfkb_activity.png", fig_inflam_box,
       width = 9, height = 5, dpi = 300)

# 7. EFTUD2 vs Pathway Scatter Plots -----
scatter_df <- data.frame(
  EFTUD2 = eftud2_expr,
  Group  = targets$group,
  t(GSVA_inflam_matrix)
) %>%
  tibble::rownames_to_column("sample_id") %>%
  pivot_longer(cols = -c(sample_id, EFTUD2, Group),
               names_to = "Pathway", values_to = "Score")

fig_eftud2_scatter <- ggplot(scatter_df, aes(x = EFTUD2, y = Score)) +
  geom_point(aes(color = Group), alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", color = "black", fill = "grey80") +
  scale_color_manual(values = c("NN" = "#2C467A", "PP" = "#BE684D")) +
  facet_wrap(~ Pathway, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(title = "EFTUD2 Expression vs Inflammatory Pathway Activity",
       x     = "EFTUD2 Expression (log2 CPM)",
       y     = "GSVA Score")

ggsave("fig_eftud2_vs_inflam.png", fig_eftud2_scatter,
       width = 9, height = 7, dpi = 300)

cat("\nStep 11 complete. Outputs:\n")
cat("  - il17_nfkb_enrichment.tsv\n")
cat("  - eftud2_il17_nfkb_correlations.tsv\n")
cat("  - fig_il17_nfkb_activity.png\n")
cat("  - fig_eftud2_vs_inflam.png\n")
