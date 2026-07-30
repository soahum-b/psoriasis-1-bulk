# 1. Setup & Load Packages -----

library(tidyverse)    # Core data science toolkit
library(edgeR)        # DGEList object and normalization methods
library(matrixStats)  # Summary stats
library(cowplot)      # Combine multiple ggplots into one figure
library(recount3)     # Access to recount3 RSE objects

# Extract sample labels early for easy access later
sampleLabels <- colData(rse_gene)$external_id


# 2. Gene Annotation & Matching-----
# # Extract, filter, and keep only unique gene symbols
gene_anno_unique <- as_tibble(rowData(rse_gene)) %>%
  mutate(gene_id = rownames(rse_gene)) %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::filter(!is.na(gene_name), gene_name != "") %>%
  dplyr::distinct(gene_name, .keep_all = TRUE)

# Align the annotation rows to match the exact order of myCounts
gene_anno_matched <- gene_anno_unique[match(rownames(myCounts), gene_anno_unique$gene_name), ]

# Quick validation checks (Both should be TRUE)
stopifnot(nrow(myCounts) == nrow(gene_anno_matched))
stopifnot(all(rownames(myCounts) == gene_anno_matched$gene_name))


# 3. Build & Save Initial DGEList-----
myDGEList <- DGEList(
    counts = myCounts,                 
    samples = targets,                 
    group = targets$group,             
    genes = gene_anno_matched          
)

save(myDGEList, file = "myDGElist.RData")


# 4. Helper Function for Plotting Log2 CPM -----

# This function handles all the reshaping and plotting so we don't repeat code.
plot_cpm_violin <- function(dge_obj, plot_title, plot_subtitle, labels) {
  
  # Calculate Log2 CPM
  log2_cpm <- edgeR::cpm(dge_obj, log = TRUE)
  
  # Convert to tibble and set column names
  df <- as_tibble(log2_cpm, rownames = "geneID")
  colnames(df) <- c("geneID", labels)
  
  # Pivot to long format
  df_pivot <- pivot_longer(df, cols = -1, 
                           names_to = "samples", 
                           values_to = "expression")
  
  # Generate Plot
  ggplot(df_pivot, aes(x = samples, y = expression, fill = samples)) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median", geom = "point", shape = 95, 
                 size = 10, color = "black", show.legend = FALSE) +
    labs(y = "log2 expression", x = "sample",
         title = plot_title,
         subtitle = plot_subtitle,
         caption = paste0("produced on ", Sys.Date())) +
    theme_bw()
}


# 5. Filtering Lowly Expressed Genes -----
## Keep genes with >1 CPM in at least 3 samples (adjust 3 based on smallest group) -----
keepers <- rowSums(edgeR::cpm(myDGEList) >= 1) >= 60 # to account for highly dynamic, disease-specific transcripts that suffer from technical dropout
myDGEList.filtered <- myDGEList[keepers, ]

# Summary of filtering
total_genes <- nrow(myDGEList)
deleted_genes <- rownames(myDGEList)[!keepers]

# Export deleted genes to CSV
if (!dir.exists("deleted_genes")) dir.create("deleted_genes")
write.csv(data.frame(DeletedGenes = deleted_genes), 
          file = "deleted_genes/filtered_out_genes.csv", 
          row.names = FALSE)

# Print filtering stats to console
cat(sprintf("Total genes: %d\n", total_genes))
cat(sprintf("Genes deleted: %d (%.2f%%)\n", length(deleted_genes), (length(deleted_genes) / total_genes) * 100))
cat(sprintf("Genes kept: %d (%.2f%%)\n", nrow(myDGEList.filtered), (nrow(myDGEList.filtered) / total_genes) * 100))


# 6. Normalization (TMM)-------
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")


# # 7. Visualization Output ------
# ## Generate the three plots using our custom helper function

# p1 <- plot_cpm_violin(myDGEList, "Log2 Counts per Million (CPM)", "Unfiltered, non-normalized", labels = sampleLabels)
# p2 <- plot_cpm_violin(myDGEList.filtered, "Log2 Counts per Million (CPM)", "Filtered, non-normalized", labels = sampleLabels)
# p3 <- plot_cpm_violin(myDGEList.filtered.norm, "Log2 Counts per Million (CPM)", "Filtered & Normalized", labels = sampleLabels)

# library(ggplot2)

# # Save each plot to your working directory
# # ggsave("Log2_CPM_.png", plot = p1, width = 90, height = 5, dpi = 300, limitsize = FALSE)
# ggsave("Log2CPM_filtered.png", plot = p2, width = 90, height = 5, dpi = 300, limitsize = FALSE)
# ggsave("Log2CPM_Filtered_Normalized.png", plot = p3, width = 90, height = 5, dpi = 300, limitsize = FALSE)

# ## Combine into a single grid ----
# combined_grid <- plot_grid(p1, p2, p3, 
#                            labels = c('A', 'B', 'C'), 
#                            label_size = 12, 
#                            ncol = 3) # Ensures they stay in one row

# ggsave("Combined_CPM_Plots.png", 
#        plot = combined_grid, 
#        scale = 1, 
#       #  height = 5, 
#       #  dpi = 300,
#       limitsize = FALSE)
