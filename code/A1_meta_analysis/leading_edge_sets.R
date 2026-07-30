# Auto-extracted generating script
# Produces: leading_edge_sets.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_results_all.rds
# Source artifact version: c04465db-f557-4dc1-b802-76d653d0998e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(msigdbr); library(fgsea); library(data.table)})

# Load GSEA results
gsea_all <- readRDS("gsea_results_all.rds")

gh <- as.data.frame(gsea_all$H)

pick <- c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_COMPLEMENT")
le <- setNames(lapply(pick, function(p) {
  x <- gh$leadingEdge[gh$pathway==p][[1]]
  if (is.list(x)) x <- x[[1]]
  as.character(x)
}), gsub("HALLMARK_","",pick))

saveRDS(le, "leading_edge_sets.rds")