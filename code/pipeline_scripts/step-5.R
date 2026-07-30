# 5. Hierarchical Clustering and Heatmaps ----

# Load packages ----
library(tidyverse)
library(limma) 
library(gplots)       # For heatmap.2
library(RColorBrewer) # For color palettes

# 1. Prepare the Data ----
# We use the decideTests output from Code 4 to get our significant genes
# (Note: lfc = 1 means a 2-fold change. If you want stricter, change to 2)
results <- decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 0.01, lfc = 1)

# Ensure the column names of our expression matrix match our sample IDs
colnames(v.DEGList$E) <- paste(targets$group, targets$external_id, sep = "_")

# Subset the expression matrix to keep ONLY the significant genes
diffGenes <- v.DEGList$E[results[,1] != 0, ]

# 2. Cluster DEGs (Rows) and Samples (Columns) ----
# Cluster genes (rows) by Pearson correlation (measures trend similarity)
clustRows <- hclust(as.dist(1 - cor(t(diffGenes), method = "pearson")), method = "complete") 

# Cluster samples (columns) by Spearman correlation (robust to outliers)
clustColumns <- hclust(as.dist(1 - cor(diffGenes, method = "spearman")), method = "complete")

# Cut the resulting gene tree into 2 distinct clusters (Up-regulated vs Down-regulated modules)
module.assign <- cutree(clustRows, k = 2)

# Assign a color to each module for the heatmap sidebar
module.color <- rainbow(length(unique(module.assign)), start = 0.1, end = 0.9) 
module.color <- module.color[as.vector(module.assign)] 

# 3. Choose Color Palette ----
# Using a classic Red/Blue palette (Red = high expression, Blue = low expression)
myheatcolors <- rev(brewer.pal(name = "RdBu", n = 11))

# 4. Produce the Main Heatmap ----
# This plots all significant DEGs across all samples
heatmap.2(diffGenes, 
          Rowv = as.dendrogram(clustRows), 
          Colv = as.dendrogram(clustColumns),
          RowSideColors = module.color,
          col = myheatcolors, 
          scale = 'row',     # Scales expression by gene to highlight differences
          labRow = NA,       # Hides gene names (too messy with thousands of genes)
          density.info = "none", 
          trace = "none",  
          cexCol = 0.5,      # Makes sample labels smaller to fit
          margins = c(10, 5),
          main = "All DEGs: Psoriasis vs Normal") 

# 5. Extract and View Specific Modules ----

# Let's pull out Module 2 (For example, genes that might be strictly UP in Psoriasis)
modulePick <- 2 

# Subset the matrix for only genes in Module 2
myModule_2 <- diffGenes[names(module.assign[module.assign %in% modulePick]), ] 

# Re-cluster just this subset
hrsub_2 <- hclust(as.dist(1 - cor(t(myModule_2), method = "pearson")), method = "complete") 

# Heatmap for just Module 2
heatmap.2(myModule_2, 
          Rowv = as.dendrogram(hrsub_2), 
          Colv = as.dendrogram(clustColumns), 
          labRow = NA,
          col = myheatcolors, 
          scale = "row", 
          density.info = "none", 
          trace = "none", 
          RowSideColors = module.color[module.assign %in% modulePick], 
          margins = c(10, 5),
          main = "Module 2 Sub-cluster")

# 6. Export Modules for Downstream Analysis ----
# Export Module 2 to a text file for things like Pathway Analysis (e.g., Enrichr, STRING)
moduleSymbols <- tibble(geneID = rev(hrsub_2$labels[hrsub_2$order]))
moduleData <- diffGenes[moduleSymbols$geneID, ]
moduleData.df <- as_tibble(moduleData, rownames = "geneSymbol")

write_tsv(moduleData.df, "Module_2_Genes.tsv")



# 5. Focused Hierarchical Clustering and Heatmaps ----

# 1. Define and Prepare the Targeted Data ----
# Define the specific genes related to the requested pathways
target_genes <- c(
  "EFTUD2", "STAT3",                       # Core requested genes
  "RELA", "RELB", "REL", "NFKB1", "NFKB2",  # NF-kB transcription factor family
  "IL17A", "IL17F", "IL17RA", "IL17RC"     # IL-17 ligands and receptors
)

# Ensure sample names are formatted correctly
colnames(v.DEGList$E) <- paste(targets$group, targets$external_id, sep = "_")

# Subset the matrix: Only keep target genes that exist in your dataset
# This prevents errors if a specific gene was filtered out during preprocessing
available_targets <- intersect(target_genes, rownames(v.DEGList$E))
target_matrix <- v.DEGList$E[available_targets, ]

# 2. Cluster Targeted Genes and Samples ----
# Gene clustering (Rows)
clustRows_target <- hclust(as.dist(1 - cor(t(target_matrix), method = "pearson")), method = "complete") 

# Sample clustering (Columns)
clustCols_target <- hclust(as.dist(1 - cor(target_matrix, method = "spearman")), method = "complete")

# 3. Choose Color Palette ----
myheatcolors <- rev(brewer.pal(name = "RdBu", n = 11))

# 4. Produce the Targeted Heatmap ----
# Note: labRow is set to available_targets to display gene names
# heatmap_3 <- heatmap.2(target_matrix, 
#           Rowv = as.dendrogram(clustRows_target), 
#           Colv = as.dendrogram(clustCols_target),
#           col = myheatcolors, 
#           scale = 'row',           # Z-score scaling to show relative changes
#           labRow = rownames(target_matrix), # Show gene symbols
#           cexRow = 0.8,            # Adjust font size for gene names
#           cexCol = 0.7,            # Adjust font size for sample names
#           density.info = "none", 
#           trace = "none",  
#           margins = c(10, 8),      # Wider margins for labels
#           main = "Targeted Pathway Expression\n(EFTUD2, STAT3, NF-kB, IL-17)")

# ggsave('heatmap-target-genes.png',heatmap_3, scale = 1)

# Open the file connection
png("heatmap-target-genes.png", width = 10, height = 10, units = "in", res = 300)

# Run the heatmap code (do not assign it to a variable)
heatmap.2(target_matrix, 
          Rowv = as.dendrogram(clustRows_target), 
          Colv = as.dendrogram(clustCols_target),
          col = myheatcolors, 
          scale = 'row',
          labRow = rownames(target_matrix),
          cexRow = 0.8,
          cexCol = 0.7,
          density.info = "none", 
          trace = "none",  
          margins = c(10, 8),
          main = "Targeted Pathway Expression\n(EFTUD2, STAT3, NF-kB, IL-17)")

# Crucial: Close the device to save the file
dev.off()


# 1. Subset the expression matrix for the 9 found genes
# Use the 'available_targets' variable we verified in the diagnostic
target_matrix <- v.DEGList$E[available_targets, ]

# 2. Align Column Names and Annotation
# We must ensure colnames(target_matrix) matches rownames(anno_col)
# Let's use the external_id as the primary key
colnames(target_matrix) <- targets$external_id

anno_col <- data.frame(Group = targets$group)
rownames(anno_col) <- targets$external_id

# 3. Generate the Heatmap (Plotting to the Viewer)
library(pheatmap)
p_heat <- pheatmap(target_matrix, 
                   scale = "row",                      # Z-score normalization
                   clustering_distance_rows = "correlation", 
                   clustering_distance_cols = "euclidean",
                   annotation_col = anno_col,          # Adds the NN vs PP color bar
                   show_colnames = FALSE,              # Sample names are often too crowded
                   main = "Psoriasis Pathway Analysis: NN vs PP",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# 4. Save the plot using a dedicated device (avoiding ggsave errors)
# This creates a high-res PNG of the plot you just saw
png("Targeted_Heatmap_Final.png", width = 8, height = 6, units = "in", res = 300)
grid::grid.newpage()
grid::grid.draw(p_heat$gtable)
dev.off()