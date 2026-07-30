## 1. Load Packages
library(recount3)
library(edgeR)       
library(dplyr)       
library(rhdf5)

## 2. Find and Load Human Project
human_projects <- available_projects()

proj_info <- subset(
    human_projects,
    project == "SRP035988" & project_type == "data_sources"
)

# Create RSE at the gene level
rse_gene <- create_rse(proj_info)

## 3. Transform Counts (CRITICAL STEP)
# recount3 provides raw coverage; we must transform it to estimated read counts
# We assign it directly back to the RSE object's "counts" assay
assay(rse_gene, "counts") <- transform_counts(rse_gene)

## 4. Clean Metadata & Assign Groups
# Work directly with colData to avoid creating separate, disconnected dataframes
colData(rse_gene)$group <- NA

# Assign groups based on attributes (Ensure exact string match with your dataset)
is_normal <- grepl("normal_skin", colData(rse_gene)$sra.sample_attributes, ignore.case = TRUE)
is_psoriasis <- grepl("Psoriasis_skin", colData(rse_gene)$sra.sample_attributes, ignore.case = TRUE)

colData(rse_gene)$group[is_normal] <- "NN"
colData(rse_gene)$group[is_psoriasis] <- "PP"

# Convert to factor
colData(rse_gene)$group <- as.factor(colData(rse_gene)$group)

# Verify grouping
table(colData(rse_gene)$group)

## 5. Extract and Clean the Count Matrix
myCounts <- assay(rse_gene, "counts")
gene_symbols <- rowData(rse_gene)$gene_name

# Filter out missing or blank gene symbols
valid <- !is.na(gene_symbols) & gene_symbols != ""
myCounts <- myCounts[valid, ]
gene_symbols <- gene_symbols[valid]

# Collapse duplicated symbols by summing counts
myCounts <- rowsum(myCounts, group = gene_symbols)

## 6. Build the DGEList for edgeR
# edgeR requires integer counts, so we must round the transform_counts() output
myDGEList <- DGEList(
    counts = round(myCounts), 
    group = colData(rse_gene)$group
)

## 7. Build a Clean Targets Dataframe 
targets <- as.data.frame(colData(rse_gene)) %>%
  dplyr::select(
    external_id,              
    sra.sample_title,         
    group,    
    sra.experiment_attributes 
  )