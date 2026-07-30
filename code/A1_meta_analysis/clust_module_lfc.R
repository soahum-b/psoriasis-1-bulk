# Auto-extracted generating script
# Produces: clust_module_lfc.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds
# Source artifact version: 9e4f2014-6045-4f96-a4e3-ca73d6991292
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)

de <- readRDS("per_study_de.rds")

lines <- readLines("clust_out3/Clusters_Objects.tsv")
tab <- fread("clust_out3/Clusters_Objects.tsv", skip=1, header=TRUE)
mod <- tab[[1]]; mod <- mod[mod!="" & !is.na(mod)]

studies_use <- names(de$PPvsNN)
lfc_mat <- matrix(NA_real_, length(mod), length(studies_use), dimnames=list(mod, studies_use))
for (s in studies_use) {
  x <- de$PPvsNN[[s]]; x <- x[!duplicated(x$gene),]
  lfc_mat[, s] <- x$logFC[match(mod, x$gene)]
}
lfc_mat <- lfc_mat[order(-rowMeans(lfc_mat, na.rm=TRUE)),,drop=FALSE]

saveRDS(lfc_mat, "clust_module_lfc.rds")