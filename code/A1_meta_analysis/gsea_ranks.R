# Auto-extracted generating script
# Produces: gsea_ranks.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds
# Source artifact version: acdaa257-263c-43dc-9ac9-06cf76e0dc3f
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR); library(limma)})

de <- readRDS("de_results_full.rds")

ranks <- de$t
names(ranks) <- de$gene
ranks <- sort(ranks, decreasing = TRUE)
saveRDS(ranks, "gsea_ranks.rds")