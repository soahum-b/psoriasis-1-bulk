# Multivariate analysis ----
# Comparison: PP (lesional) vs normal NN (study 2)


# 1. Setup & Load Packages -----
library(tidyverse)   # Core data science toolkit
library(DT)          # Interactive tables
library(plotly)      # Interactive plots
library(gt)          # Publication-quality tables

## Ensure group is a factor with correct reference levels (PL as reference)-----
is.factor(targets$group)
group = targets$group
# Returns TRUE or FALSE
### returns true skip next part
### skip targets$group <- factor(targets$group, levels = c("NN", "PP")) -----
## group <- targets$group

# Update sample selection for averages (matching your Code 1 column names)
NN_samples <- targets$external_id[targets$group == "NN"] # Changed target$sample to target$external_id
PP_samples <- targets$external_id[targets$group == "PP"]

# 2. Principal Component Analysis (PCA) ----

### Generate Log2 CPM matrix for PCA and downstream stats -----
# This applies the TMM normalization factors and log2 transforms the data
log2.cpm.filtered.norm <- edgeR::cpm(myDGEList.filtered.norm, log = TRUE)
# Create the dataframe version needed multivariate analysis
log2.cpm.filtered.norm.df <- as.data.frame(log2.cpm.filtered.norm) %>%
  rownames_to_column(var = "geneID")
# Run PCA on the normalized log2 CPM, NOT raw counts
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale. = FALSE, retx = TRUE)

# Calculate percentage of variance explained by each PC
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var / sum(pc.var) * 100, 1)

# Prepare a single PCA dataframe for plotting
pca.res.df <- as_tibble(pca.res$x) %>%
  add_column(sample = paste(targets$group, targets$patient_num, sep = "_"),
             group  = group)

# Main PCA Plot (PC1 vs PC2)
pca_plot <- ggplot(pca.res.df, aes(x = PC1, y = PC2, label = sample, color = group)) +
  geom_point(size = 4) +
  geom_label(vjust = 1.5, show.legend = FALSE) +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(title = "PCA plot — Study-2 (NN vs PP)",
       caption = paste0("produced on ", Sys.Date())) +
  theme_bw()

print(pca_plot)

# PCA 'Small Multiples' Chart (PC1 to PC4)
pca.pivot <- pca.res.df %>%
  dplyr::select(sample, group, PC1:PC4) %>%
  pivot_longer(cols = PC1:PC4, names_to = "PC", values_to = "loadings")

ggplot(pca.pivot, aes(x = sample, y = loadings, fill = group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~PC) +
  labs(title = "PCA 'Small Multiples' Plot",
       caption = paste0("produced on ", Sys.Date())) +
  theme_bw() +
  coord_flip()


# 3. Summary Statistics (Averages & LogFC)-----
# Dynamically identify samples for each group
# Use external_id as the sample identifier
NN_samples <- targets$external_id[targets$group == "NN"]
PP_samples <- targets$external_id[targets$group == "PP"]

# Initialize a new column
targets$new_id <- NA

# Assign the new IDs to the correct rows
targets$new_id[targets$group == "NN"] <- paste0("NN_", 1:sum(targets$group == "NN"))
targets$new_id[targets$group == "PP"] <- paste0("PP_", 1:sum(targets$group == "PP"))

# Compute row means and LogFC (Positive = Up in Lesional/PP)
mydata.df <- log2.cpm.filtered.norm.df %>%
  mutate(
    healthy.AVG = rowMeans(log2.cpm.filtered.norm[, NN_samples, drop = FALSE]),
    disease.AVG = rowMeans(log2.cpm.filtered.norm[, PP_samples, drop = FALSE]),
    LogFC       = disease.AVG - healthy.AVG
  ) %>%
  mutate_if(is.numeric, round, 2)


# 4. Filtering and Tables ----
# # ------------------------------------------------------------------------------
# # this example is from the paper by Dr. Daniel. can skip filtering by genes of interest. 
# # Set up publication table using `gt`
# genes_of_interest <- c("MMP1", "GZMB", "IL1B", "GNLY", "IFNG", "CCL4", "PRF1", "APOBEC3A", "UNC13A")

# mydata.filter <- mydata.df %>%
#   dplyr::filter(geneID %in% genes_of_interest) %>%
#   dplyr::select(geneID, healthy.AVG, disease.AVG, LogFC) %>%
#   dplyr::arrange(desc(LogFC))


# mydata.filter %>%
#   gt() %>%
#   fmt_number(columns = 2:4, decimals = 1) %>%
#   tab_header(
#     title = md("**Regulators of skin pathogenesis**"),
#     subtitle = md("*during cutaneous leishmaniasis (SRP016583)*")
#   ) %>%
#   tab_source_note(source_note = md("Data: recount3 project SRP016583"))

# Set up interactive searchable table using `DT`
datatable(mydata.df,
          extensions = c('KeyTable', "FixedHeader"),
          filter = 'top',
          # server = TRUE,
          options = list(
            keys = TRUE,
            searchHighlight = TRUE,
            pageLength = 10,
            lengthMenu = c(10, 25, 50, 100)
          ))

mydata_interac <-datatable(mydata.df,
          options = list(
            dom = 't',       # Only show the 't'able (removes search, paging, etc.)
            ordering = FALSE, # Disables column sorting
            paging = FALSE    # Shows all rows at once
          ))

library(knitr)
kable(mydata.df)

# 5. Interactive Scatter Plot
# ------------------------------------------------------------------------------
scatter_plot <- ggplot(mydata.df, aes(x = healthy.AVG, y = disease.AVG, text = paste("Symbol:", geneID))) +
  geom_point(shape = 16, size = 1, alpha = 0.5) +
  ggtitle("Disease vs. Healthy Expression") +
  theme_bw()

# Convert to interactive plotly graphic
ggsave("scatter_plot_1.png", plot = scatter_plot, width = 90, height = 5, dpi = 300, limitsize = FALSE)



# 4. Filtering and Tables ----
# ------------------------------------------------------------------------------
# Define the specific genes requested: EFTUD2, STAT3, NF-kB, and IL-17 related genes
genes_of_interest <- c(
  "EFTUD2", "STAT3",                       # Requested specific genes
  "RELA", "RELB", "REL", "NFKB1", "NFKB2",  # NF-kB transcription factor family
  "IL17A", "IL17F", "IL17RA", "IL17RC"     # IL-17 ligands and receptors
)

# Filter the main dataframe for these specific genes
mydata.filter <- mydata.df %>%
  dplyr::filter(geneID %in% genes_of_interest) %>%
  dplyr::select(geneID, healthy.AVG, disease.AVG, LogFC) %>%
  dplyr::arrange(desc(LogFC))

# Set up publication-quality table using `gt`
mydata.filter %>%
  gt() %>%
  fmt_number(columns = 2:4, decimals = 2) %>%
  tab_header(
    title = md("**Expression of Targeted Regulatory Genes**"),
    subtitle = md("*Analysis of EFTUD2, STAT3, NF-κB, and IL-17 pathways*")
  ) %>%
  tab_source_note(source_note = md("Data: Study-2 (NN vs PP) normalized log2 CPM"))

# Set up interactive searchable table for the filtered set
datatable(mydata.filter,
          extensions = c('KeyTable', "FixedHeader"),
          filter = 'top',
          options = list(
            keys = TRUE,
            searchHighlight = TRUE,
            pageLength = 10
          ))

# 5. Interactive Scatter Plot
# ------------------------------------------------------------------------------
# Highlighting the genes of interest in the scatter plot
scatter_plot <- ggplot(mydata.df, aes(x = healthy.AVG, y = disease.AVG)) +
  geom_point(shape = 16, size = 1, alpha = 0.3, color = "grey") +
  geom_point(data = mydata.filter, aes(color = "Target Genes"), size = 2) +
  geom_text(data = mydata.filter, aes(label = geneID), vjust = -1, size = 3) +
  scale_color_manual(values = c("Target Genes" = "red")) +
  ggtitle("Disease vs. Healthy Expression: Pathway Targets") +
  theme_bw()

# Save the plot
ggsave("scatter_plot_targets.png", plot = scatter_plot, width = 10, height = 8, dpi = 300)