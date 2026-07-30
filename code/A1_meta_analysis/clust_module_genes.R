# Auto-extracted generating script
# Produces: clust_module_genes.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): meta_de_results.rds
# Source artifact version: 8a10bd30-09f6-42a0-ad27-3a6a67c72bea
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(msigdbr)

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/")

meta <- readRDS("meta_de_results.rds")
mp <- as.data.table(meta$PPvsNN)
focus <- mp[FDR<0.05 & abs(logFC)>1, gene]

dir.create("clust_input2", showWarnings=FALSE)

lines <- readLines("clust_out3/Clusters_Objects.tsv")
tab <- fread("clust_out3/Clusters_Objects.tsv", skip=1, header=TRUE)
mod <- tab[[1]]; mod <- mod[mod!="" & !is.na(mod)]

fwrite(data.table(gene=mod, logFC=mp$logFC[match(mod,mp$gene)]), "clust_module_genes.csv")