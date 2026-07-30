# =============================================================================
# Step 10 (Revised): STAT3 Alpha/Beta Isoform Switch — Recount3 Data Only
# =============================================================================
# Confirmed findings from diagnostic run:
#
#   ALPHA junction  ExonID 2: chr17:42,316,852-42,317,181 (PP-enriched, FDR=0.59)
#   BETA  junction  ExonID 3: chr17:42,316,902-42,317,181 (NN-enriched, FDR=0.02)
#   Shared acceptor at chr17:42,317,181; donors differ by 50 bp.
#
#   ALPHA exon bin  idx 5:  chr17:42,313,324-42,315,800 (2477 bp, long terminal exon)
#   BETA  exon bin  idx 2:  chr17:42,315,559-42,315,800 (242 bp,  shared short exon)
#
# Requires: stat3_splicing_objects.RData (from Step 9)
#           targets dataframe (from Step 1)
# =============================================================================

library(tidyverse)
library(ggplot2)
library(GenomicRanges)
library(SummarizedExperiment)

load("stat3_splicing_objects.RData")


# =============================================================================
# PART 1: MAP THE EXON LANDSCAPE
# =============================================================================
cat("=== ALL STAT3 EXON BINS (sorted low to high genomic coord) ===\n")
cat("LOW coord = 3' end of transcript (terminal exon, minus strand)\n\n")

all_bins <- as.data.frame(rowRanges(rse_exon_stat3), row.names = NULL) %>%
  mutate(bin_idx = seq_len(n())) %>%
  arrange(start) %>%
  mutate(transcript_order = rev(seq_len(n())))

print(all_bins[, c("bin_idx","seqnames","start","end","width","strand","transcript_order")])

cat("\n=== SIGNIFICANT EXON BINS (FDR < 0.05) ===\n")
sig_exons <- read_tsv("stat3_significant_exons.tsv", show_col_types = FALSE)
sig_with_order <- sig_exons %>%
  left_join(all_bins %>% dplyr::select(bin_idx, start, end, transcript_order),
            by = c("ExonID" = "bin_idx")) %>%
  arrange(Start)
print(sig_with_order)


# =============================================================================
# PART 2: MAP THE JUNCTION LANDSCAPE
# =============================================================================
cat("\n=== ALL STAT3 JUNCTIONS (sorted low to high genomic coord) ===\n")
cat("Minus strand: Start (low) = 3' side of intron; End (high) = 5' side\n\n")

all_jxn <- read_tsv("stat3_junctions_with_coords.tsv", show_col_types = FALSE) %>%
  arrange(Start) %>%
  mutate(span_kb   = round((End - Start) / 1000, 1),
         sig       = ifelse(FDR < 0.05, "***", ""),
         direction = ifelse(logFC > 0, "PP-enriched", "NN-enriched"))

print(all_jxn %>% dplyr::select(ExonID, Start, End, span_kb, logFC, FDR, sig, direction))

cat("\n=== SIGNIFICANT JUNCTIONS ONLY ===\n")
print(all_jxn %>% filter(FDR < 0.05) %>%
        dplyr::select(ExonID, Start, End, span_kb, logFC, FDR, direction))


# =============================================================================
# PART 3: COMPETING JUNCTION PAIRS (shared splice site)
# =============================================================================
cat("\n=== PAIRS SHARING THE SAME END (same upstream donor) ===\n")
shared_end <- all_jxn %>%
  group_by(End) %>% filter(n() > 1) %>%
  arrange(End, Start) %>%
  mutate(donor_diff_bp = End - Start)
print(shared_end %>% dplyr::select(ExonID, Start, End, span_kb, logFC, FDR, sig, direction))

cat("\n=== PAIRS SHARING THE SAME START (same downstream acceptor) ===\n")
shared_start <- all_jxn %>%
  group_by(Start) %>% filter(n() > 1) %>%
  arrange(Start, End)
print(shared_start %>% dplyr::select(ExonID, Start, End, span_kb, logFC, FDR, sig, direction))


# =============================================================================
# PART 4: TERMINAL EXON REGION JUNCTIONS
# =============================================================================
terminal_lower <- min(all_bins$start)
terminal_upper <- terminal_lower + 5000

cat(sprintf("\n=== JUNCTIONS IN TERMINAL EXON REGION chr17:%d-%d ===\n",
            terminal_lower, terminal_upper))
terminal_jxn <- all_jxn %>%
  filter(Start <= terminal_upper | End <= terminal_upper) %>%
  arrange(Start)
print(terminal_jxn %>% dplyr::select(ExonID, Start, End, span_kb, logFC, FDR, sig, direction))


# =============================================================================
# PART 5: CONFIRMED IDs
# =============================================================================

# CONFIRMED alpha/beta junction IDs from Part 3 shared-End analysis:
alpha_junction_ExonID <- 2   # chr17:42,316,852-42,317,181 (PP-enriched, FDR=0.59)
beta_junction_ExonID  <- 3   # chr17:42,316,902-42,317,181 (NN-enriched, FDR=0.02)

# CONFIRMED exon bin indices from Part 1 terminal_bins output:
alpha_long_idx <- 5   # chr17:42,313,324-42,315,800  (2477 bp, full alpha terminal exon)
beta_short_idx <- 2   # chr17:42,315,559-42,315,800  (242 bp,  shared short exon)

cat("\n=== CONFIRMED IDENTITIES ===\n")
cat(sprintf("Alpha junction  ExonID %d: chr17:%d-%d\n",
            alpha_junction_ExonID,
            all_jxn$Start[all_jxn$ExonID == alpha_junction_ExonID],
            all_jxn$End[all_jxn$ExonID   == alpha_junction_ExonID]))
cat(sprintf("Beta  junction  ExonID %d: chr17:%d-%d  (FDR=%.3f, NN-enriched)\n",
            beta_junction_ExonID,
            all_jxn$Start[all_jxn$ExonID == beta_junction_ExonID],
            all_jxn$End[all_jxn$ExonID   == beta_junction_ExonID],
            all_jxn$FDR[all_jxn$ExonID   == beta_junction_ExonID]))

alpha_coords <- as.data.frame(rowRanges(rse_exon_stat3)[alpha_long_idx])
beta_coords  <- as.data.frame(rowRanges(rse_exon_stat3)[beta_short_idx])
cat(sprintf("Alpha exon bin  idx %d:    chr17:%d-%d  (%d bp)\n",
            alpha_long_idx, alpha_coords$start, alpha_coords$end, alpha_coords$width))
cat(sprintf("Beta  exon bin  idx %d:    chr17:%d-%d  (%d bp)\n",
            beta_short_idx, beta_coords$start,  beta_coords$end,  beta_coords$width))


# =============================================================================
# PART 6: JUNCTION-BASED PSI
# Uses confirmed junction read counts directly
# =============================================================================
cat("\n=== PART 6: JUNCTION-BASED PSI ===\n")

alpha_jxn_counts <- as.numeric(stat3_jxn_counts[alpha_junction_ExonID, ])
beta_jxn_counts  <- as.numeric(stat3_jxn_counts[beta_junction_ExonID,  ])

# Beta ratio: fraction of reads supporting the beta-producing junction
# Pseudocount of 1 stabilises low-coverage samples
beta_ratio <- (beta_jxn_counts + 1) / (alpha_jxn_counts + beta_jxn_counts + 2)

psi_jxn_df <- data.frame(
  Sample      = colnames(stat3_jxn_counts),
  Group       = targets$group,
  Alpha_reads = alpha_jxn_counts,
  Beta_reads  = beta_jxn_counts,
  Total_reads = alpha_jxn_counts + beta_jxn_counts,
  Beta_Ratio  = beta_ratio
)

# Keep only samples with >= 10 reads at this locus for robustness
psi_jxn_robust <- psi_jxn_df %>% filter(Total_reads >= 10)
cat(sprintf("Samples with >= 10 reads at locus: %d / %d\n",
            nrow(psi_jxn_robust), nrow(psi_jxn_df)))

psi_jxn_test <- wilcox.test(Beta_Ratio ~ Group, data = psi_jxn_robust)
cat("\nWilcoxon test (junction-based Beta Ratio):\n")
print(psi_jxn_test)

psi_jxn_summary <- psi_jxn_robust %>%
  group_by(Group) %>%
  summarise(n           = n(),
            median_beta = round(median(Beta_Ratio), 3),
            mean_beta   = round(mean(Beta_Ratio), 3),
            .groups     = "drop")
cat("\nGroup summary (higher beta ratio = more STAT3-beta):\n")
print(psi_jxn_summary)

write_tsv(psi_jxn_df, "stat3_psi_junction_confirmed.tsv")
cat("Saved: stat3_psi_junction_confirmed.tsv\n")


# =============================================================================
# PART 7: EXON-BIN PSI (second independent method)
# Uses read density on the alpha terminal exon vs shared short exon
# =============================================================================
cat("\n=== PART 7: EXON-BIN PSI ===\n")

alpha_exon_counts <- as.numeric(assay(rse_exon_stat3[alpha_long_idx, ], "counts"))
beta_exon_counts  <- as.numeric(assay(rse_exon_stat3[beta_short_idx,  ], "counts"))

width_alpha <- alpha_coords$width   # 2477 bp
width_beta  <- beta_coords$width    #  242 bp

# Reads per bp — corrects for the 10x size difference between the two bins
alpha_density <- (alpha_exon_counts + 1) / width_alpha
beta_density  <- (beta_exon_counts  + 1) / width_beta

# Alpha-Inclusion Ratio: higher = more STAT3-alpha
alpha_ratio <- alpha_density / beta_density

psi_exon_df <- data.frame(
  Sample       = colnames(assay(rse_exon_stat3, "counts")),
  Group        = targets$group,
  Alpha_counts = alpha_exon_counts,
  Beta_counts  = beta_exon_counts,
  Alpha_Ratio  = alpha_ratio
)

psi_exon_test <- wilcox.test(Alpha_Ratio ~ Group, data = psi_exon_df)
cat("\nWilcoxon test (exon-bin Alpha-Inclusion Ratio):\n")
print(psi_exon_test)

psi_exon_summary <- psi_exon_df %>%
  group_by(Group) %>%
  summarise(n            = n(),
            median_ratio = round(median(Alpha_Ratio), 3),
            mean_ratio   = round(mean(Alpha_Ratio), 3),
            .groups      = "drop")
cat("\nGroup summary (lower ratio in PP = more STAT3-beta in psoriasis):\n")
print(psi_exon_summary)

write_tsv(psi_exon_df, "stat3_psi_exon_confirmed.tsv")
cat("Saved: stat3_psi_exon_confirmed.tsv\n")


# =============================================================================
# PART 8: VISUALISE — side-by-side comparison of both methods
# =============================================================================
cat("\n=== PART 8: FIGURES ===\n")

# Plot 1: Junction-based Beta Ratio
fig_jxn <- ggplot(psi_jxn_robust, aes(x = Group, y = Beta_Ratio, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("NN" = "#2C467A", "PP" = "#BE684D")) +
  labs(
    title    = "STAT3 Isoform Switch: Junction Evidence",
    subtitle = paste0("Beta Junction Ratio  |  Wilcoxon p = ",
                      signif(psi_jxn_test$p.value, 3),
                      "  (n = ", nrow(psi_jxn_robust), " samples, >= 10 reads)"),
    y        = "Beta Junction Ratio  [ beta reads / (alpha + beta reads) ]",
    x        = "Tissue Status",
    caption  = paste0("Alpha junction ExonID ", alpha_junction_ExonID,
                      ": chr17:42,316,852-42,317,181\n",
                      "Beta  junction ExonID ", beta_junction_ExonID,
                      ": chr17:42,316,902-42,317,181  (shared acceptor, 50 bp donor diff)")
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggsave("fig_stat3_psi_junction.png", fig_jxn, width = 6, height = 5, dpi = 300)
cat("Saved: fig_stat3_psi_junction.png\n")

# Plot 2: Exon-bin Alpha-Inclusion Ratio
fig_exon <- ggplot(psi_exon_df, aes(x = Group, y = Alpha_Ratio, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("NN" = "#2C467A", "PP" = "#BE684D")) +
  labs(
    title    = "STAT3 Isoform Switch: Exon Bin Evidence",
    subtitle = paste0("Alpha-Inclusion Ratio  |  Wilcoxon p = ",
                      signif(psi_exon_test$p.value, 3)),
    y        = "Alpha-Inclusion Ratio  [ alpha density / beta density ]",
    x        = "Tissue Status",
    caption  = paste0("Alpha bin idx ", alpha_long_idx,
                      ": chr17:42,313,324-42,315,800  (", width_alpha, " bp)\n",
                      "Beta  bin idx ", beta_short_idx,
                      ": chr17:42,315,559-42,315,800  (", width_beta,  " bp)")
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggsave("fig_stat3_psi_exon.png", fig_exon, width = 6, height = 5, dpi = 300)
cat("Saved: fig_stat3_psi_exon.png\n")

# Plot 3: Raw read count sanity check
counts_long <- psi_exon_df %>%
  dplyr::select(Sample, Group, Alpha_counts, Beta_counts) %>%
  pivot_longer(cols = c(Alpha_counts, Beta_counts),
               names_to  = "Exon", values_to = "Counts") %>%
  mutate(Exon = recode(Exon,
                       Alpha_counts = "Alpha exon bin (2477 bp)",
                       Beta_counts  = "Beta/shared exon bin (242 bp)"))

fig_check <- ggplot(counts_long, aes(x = Group, y = log2(Counts + 1), fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~Exon) +
  scale_fill_manual(values = c("NN" = "#2C467A", "PP" = "#BE684D")) +
  labs(title    = "Sanity Check: Raw Counts at Alpha and Beta Exon Bins",
       subtitle = "Both bins should have detectable reads in all samples",
       y = "log2(counts + 1)", x = "Tissue Status") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("fig_stat3_exon_counts_check.png", fig_check, width = 8, height = 5, dpi = 300)
cat("Saved: fig_stat3_exon_counts_check.png\n")

# =============================================================================
cat("\n=== Step 10 complete ===\n")
cat("Outputs:\n")
cat("  stat3_psi_junction_confirmed.tsv  — junction-based beta ratio per sample\n")
cat("  stat3_psi_exon_confirmed.tsv      — exon-bin alpha-inclusion ratio per sample\n")
cat("  fig_stat3_psi_junction.png        — junction-based isoform switch plot\n")
cat("  fig_stat3_psi_exon.png            — exon-bin isoform switch plot\n")
cat("  fig_stat3_exon_counts_check.png   — sanity check raw counts\n")
