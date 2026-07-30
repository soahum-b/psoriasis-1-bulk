# 4. Differential Gene Expression (Psoriasis vs Normal) ----

# Load packages ----
library(tidyverse) 
library(limma) 
library(edgeR)
library(DT)
library(plotly)

# Set up your design matrix ----
# Ensure 'group' is a factor
group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Model mean-variance trend and fit linear model ----
# Use voom to transform TMM-normalized counts to log2-CPM and model the mean-variance relationship
v.DEGList <- voom(myDGEList.filtered.norm, design, plot = TRUE)

# Fit the linear model
fit <- lmFit(v.DEGList, design)

# Contrast matrix and Empirical Bayes ----
contrast.matrix <- makeContrasts(Psoriasis = PP - NN, levels = design)

# Extract the linear model fit and get Bayesian stats
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)

# TopTable to view DEGs ----
# Extract all genes (number = 40000) to ensure the volcano plot is completely populated
myTopHits <- topTable(ebFit, adjust ="BH", coef = 1, number = 40000, sort.by = "logFC")

# Convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

# Volcano Plots ----
# Create the base ggplot
vplot <- ggplot(myTopHits.df) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "longdash", colour = "grey", linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "longdash", colour = "#BE684D", linewidth = 0.5) +
  geom_vline(xintercept = -1, linetype = "longdash", colour = "#2C467A", linewidth = 0.5) +
  labs(title = "Volcano Plot: Psoriasis vs Normal Skin",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       caption = paste0("produced on ", Sys.Date())) +
  theme_bw()

# Save static plot to your working directory
ggsave("Volcano_Plot_PP_vs_NN.png", plot = vplot, width = 8, height = 6, dpi = 300)

# Make the volcano plot interactive in the Viewer pane
ggplotly(vplot)

# DecideTests to pull out significant DEGs ----
# Standard thresholds: FDR < 0.01 and an absolute Log2 Fold Change of 1 (which equals a 2-fold change)
results <- decideTests(ebFit, method = "global", adjust.method = "BH", p.value = 0.01, lfc = 1)

# Summary of up/down regulated genes
print(summary(results))

# Retrieve expression data for your significant DEGs ----
# Ensure column names map correctly to your target sample IDs
colnames(v.DEGList$E) <- targets$external_id

# Subset the expression matrix to keep only genes identified as significantly different
diffGenes <- v.DEGList$E[results[,1] != 0, ]

# Convert your DEGs to a dataframe
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# Create interactive tables to display your DEGs ----
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: Significant DEGs in Psoriasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c(10, 25, 50, 100))) %>%
  formatRound(columns = 2:ncol(diffGenes.df), digits = 2)

# Write your DEGs to a file for downstream network analysis (e.g., STRING)
write_tsv(diffGenes.df, "Significant_DEGs_PPvsNN.txt")

