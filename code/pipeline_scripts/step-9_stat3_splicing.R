# 9. STAT3 Alternative Splicing: Exon and Junction Analysis ----
# Purpose: test whether STAT3 itself shows differential exon usage and
# differential junction usage between psoriatic and normal skin. The 3'
# end of STAT3 is the locus that determines the alpha (full-length) vs
# beta (truncated) isoform switch via alternative splicing of exon 23.

# 1. Setup & Load Packages -----
library(tidyverse)
library(limma)
library(edgeR)
library(GenomicRanges)
library(SummarizedExperiment)

# Requires: rse_gene (Step 1), rse_exon and rse_jxn (Step 7),
# design and contrast.matrix (Step 4), targets (Step 1).

# 2. Subset to STAT3 -----
# Get STAT3's gene-level GRanges
stat3_gr <- rowRanges(rse_gene)[rowData(rse_gene)$gene_name == "STAT3"]

# Pull the Ensembl ID (cleaning NAs) for exon-level subsetting
stat3_gene_id <- rowData(rse_gene)$gene_id[rowData(rse_gene)$gene_name == "STAT3"]
stat3_gene_id <- stat3_gene_id[!is.na(stat3_gene_id)][1]

# 3. EXON-LEVEL: Differential Exon Usage -----
# Subset the exon RSE to STAT3 exons only
stat3_exons_idx  <- which(rowData(rse_exon)$gene_id == stat3_gene_id)
rse_exon_stat3   <- rse_exon[stat3_exons_idx, ]
stat3_exon_counts <- assay(rse_exon_stat3, "counts")

# Use exon_id as row labels if available
if ("exon_id" %in% colnames(rowData(rse_exon_stat3))) {
  rownames(stat3_exon_counts) <- rowData(rse_exon_stat3)$exon_id
}

cat(sprintf("STAT3 exon bins extracted: %d\n", nrow(stat3_exon_counts)))

# Build exon-level DGEList and run voom
dge_exon <- DGEList(counts = stat3_exon_counts)
dge_exon <- calcNormFactors(dge_exon)
v_exon   <- voom(dge_exon, design, plot = FALSE)

# Fit linear model, apply contrast, then run diffSplice
# (diffSplice does its own eBayes shrinkage — do NOT call eBayes() before it)
fit_exon          <- lmFit(v_exon, design)
fit_exon_contrast <- contrasts.fit(fit_exon, contrast.matrix)
fit_splice_exon   <- diffSplice(
  fit_exon_contrast,
  geneid = rep("STAT3", nrow(stat3_exon_counts))
)

# Pull all exon stats
stat3_exon_stats <- topSplice(fit_splice_exon, coef = 1, test = "t", number = Inf)
cat("\n=== Top differentially used STAT3 exons ===\n")
print(head(stat3_exon_stats, 10))

# 4. JUNCTION-LEVEL: Differential Junction Usage -----
# Find junctions that overlap the STAT3 gene body
stat3_overlaps <- findOverlaps(rowRanges(rse_jxn), stat3_gr)
stat3_jxn_idx  <- queryHits(stat3_overlaps)
rse_jxn_stat3  <- rse_jxn[stat3_jxn_idx, ]

# Extract counts and label junctions with their genomic coordinates
stat3_jxn_counts <- assay(rse_jxn_stat3, "counts")

# Filter rare/noise junctions: >= 5 reads in at least 3 samples
keep_jxn <- rowSums(stat3_jxn_counts >= 5) >= 3
stat3_jxn_counts <- stat3_jxn_counts[keep_jxn, ]
rse_jxn_stat3    <- rse_jxn_stat3[keep_jxn, ]
rownames(stat3_jxn_counts) <- as.character(rowRanges(rse_jxn_stat3))

cat(sprintf("\nSTAT3 junctions kept after filtering: %d\n", nrow(stat3_jxn_counts)))

# Build junction-level DGEList and run voom
dge_jxn <- DGEList(counts = stat3_jxn_counts)
dge_jxn <- calcNormFactors(dge_jxn)
v_jxn   <- voom(dge_jxn, design, plot = FALSE)

fit_jxn          <- lmFit(v_jxn, design)
fit_jxn_contrast <- contrasts.fit(fit_jxn, contrast.matrix)
fit_splice_jxn   <- diffSplice(
  fit_jxn_contrast,
  geneid = rep("STAT3", nrow(stat3_jxn_counts))
)

# Pull all junction stats
top_jxn <- topSplice(fit_splice_jxn, coef = 1, test = "t", number = Inf)
cat("\n=== Top differentially used STAT3 junctions ===\n")
print(head(top_jxn, 10))

# 5. Map Significant Exons to Genomic Coordinates -----
# Ensure sig_exons exists and has data
sig_exons <- stat3_exon_stats %>% dplyr::filter(FDR < 0.05)

if (nrow(sig_exons) > 0) {
  # Robust Indexing: topSplice's ExonID usually refers to the row index
  # of the input matrix if row names were not perfectly propagated.
  exon_indices <- as.numeric(as.character(sig_exons$ExonID))
  
  # Filter out any NAs from the index conversion just in case
  valid_idx <- !is.na(exon_indices)
  sig_exons <- sig_exons[valid_idx, ]
  exon_indices <- exon_indices[valid_idx]

  # Pull ranges directly from the subsetted RSE using indices
  sig_exon_ranges <- rowRanges(rse_exon_stat3)[exon_indices]

  sig_exons_with_coords <- sig_exons %>%
    mutate(
      Chromosome = as.character(seqnames(sig_exon_ranges)),
      Start      = start(sig_exon_ranges),
      End        = end(sig_exon_ranges)
    ) %>%
    dplyr::select(ExonID, Chromosome, Start, End, logFC, t, FDR)

  cat("\n=== Significant STAT3 Exons with coordinates ===\n")
  print(sig_exons_with_coords)
  write_tsv(sig_exons_with_coords, "stat3_significant_exons.tsv")
} else {
  cat("\nNo exons reached the FDR < 0.05 threshold.\n")
}

# 6. Map Top Junctions to Genomic Coordinates -----
# Use a similar index-based approach for junctions
jxn_indices <- as.numeric(as.character(top_jxn$ExonID))
valid_jxn_idx <- !is.na(jxn_indices)

top_jxn_ranges <- rowRanges(rse_jxn_stat3)[jxn_indices[valid_jxn_idx]]

top_jxn_with_coords <- top_jxn[valid_jxn_idx, ] %>%
  mutate(
    Chromosome = as.character(seqnames(top_jxn_ranges)),
    Start      = start(top_jxn_ranges),
    End        = end(top_jxn_ranges)
  ) %>%
  dplyr::select(ExonID, Chromosome, Start, End, logFC, t, FDR)

cat("\n=== Top STAT3 junctions with coordinates ===\n")
print(head(top_jxn_with_coords, 10))
write_tsv(top_jxn_with_coords, "stat3_junctions_with_coords.tsv")

# 7. Isolate the Alpha/Beta Locus (chr17:42,313,000-42,316,000, hg38) -----
# This is the region containing the alternative 5' splice site in exon 23
# that determines STAT3-alpha vs STAT3-beta.
stat3_ab_locus <- top_jxn_with_coords %>%
  dplyr::filter((Start >= 42313000 & Start <= 42316000) |
                (End   >= 42313000 & End   <= 42316000)) %>%
  arrange(FDR)

cat("\n=== Junctions at the STAT3 alpha/beta locus ===\n")
print(stat3_ab_locus)
write_tsv(stat3_ab_locus, "stat3_alpha_beta_junctions.tsv")

# 8. Visualize: Differential Exon Usage Across STAT3 -----
stat3_exon_stats_plot <- stat3_exon_stats %>%
  mutate(
    ExonID_num  = suppressWarnings(as.numeric(as.character(ExonID))),
    Significant = ifelse(FDR < 0.05, "Significant", "Not Significant")
  ) %>%
  arrange(ExonID_num)

fig_exon_usage <- ggplot(stat3_exon_stats_plot,
                         aes(x = factor(ExonID, levels = ExonID),
                             y = logFC, fill = Significant)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Not Significant" = "grey70",
                               "Significant" = "#BE684D")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(title    = "Differential Exon Usage in STAT3",
       subtitle = "Psoriasis (PP) vs Normal Skin (NN)",
       x        = "Exon Bin (5' -> 3')",
       y        = "Log2 Fold Change (Exon vs Gene baseline)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("fig_stat3_exon_usage.png", fig_exon_usage,
       width = 10, height = 5, dpi = 300)

# 9. Visualize: Top Junction Usage Forest Plot -----
forest_data <- top_jxn %>%
  head(20) %>%
  mutate(
    SE          = abs(logFC / t),
    CI_low      = logFC - 1.96 * SE,
    CI_high     = logFC + 1.96 * SE,
    Significant = ifelse(FDR < 0.05, "Significant", "Not Significant"),
    Junction    = reorder(factor(ExonID), logFC)
  )

fig_jxn_forest <- ggplot(forest_data,
                         aes(x = logFC, y = Junction, color = Significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high), size = 0.7, linewidth = 1) +
  scale_color_manual(values = c("Not Significant" = "grey70",
                                "Significant" = "#BE684D")) +
  labs(title    = "Differential Junction Usage in STAT3",
       subtitle = "Psoriasis vs Normal (Log2 FC ± 95% CI), top 20",
       x        = "Log2 Fold Change (Junction vs Gene baseline)",
       y        = "Junction (genomic coordinates)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave("fig_stat3_junction_forest.png", fig_jxn_forest,
       width = 9, height = 7, dpi = 300)

# 10. Save Objects for Step 10 -----
save(rse_exon_stat3, rse_jxn_stat3,
     stat3_exon_stats, top_jxn, top_jxn_with_coords,
     stat3_ab_locus, stat3_jxn_counts,
     file = "stat3_splicing_objects.RData")

cat("\nStep 9 complete. Outputs:\n")
cat("  - stat3_significant_exons.tsv\n")
cat("  - stat3_junctions_with_coords.tsv\n")
cat("  - stat3_alpha_beta_junctions.tsv\n")
cat("  - fig_stat3_exon_usage.png\n")
cat("  - fig_stat3_junction_forest.png\n")
cat("  - stat3_splicing_objects.RData (for Step 10)\n")

cat("\n>>> NEXT: inspect 'stat3_alpha_beta_junctions.tsv' against GENCODE\n")
cat(">>> STAT3-201 (alpha) and STAT3-203 (beta) before running Step 10.\n")
