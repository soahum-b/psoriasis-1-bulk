# Auto-extracted generating script
# Produces: tf_activity_perstudy.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): collectri_regulon.rds, per_study_de.rds
# Source artifact version: 1fa473be-7455-47ed-ae1a-de8d85941ed3
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(decoupleR); library(data.table)})

de <- readRDS("per_study_de.rds")
reg <- readRDS("collectri_regulon.rds")
reg_df <- as.data.frame(reg)

tf_activity_study <- function(dt) {
  m <- matrix(dt$t, ncol=1, dimnames=list(dt$gene, "t"))
  res <- run_ulm(mat=m, network=reg_df, .source="source", .target="target", .mor="mor", minsize=5)
  as.data.table(res)
}

tf_by_contrast <- list()
for (cn in names(de)) {
  perstudy <- list()
  for (s in names(de[[cn]])) {
    dt <- as.data.table(de[[cn]][[s]])
    act <- tf_activity_study(dt)
    act[, study := s]
    perstudy[[s]] <- act
  }
  tf_by_contrast[[cn]] <- rbindlist(perstudy)
}
saveRDS(tf_by_contrast, "tf_activity_perstudy.rds")