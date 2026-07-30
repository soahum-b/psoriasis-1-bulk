# Auto-extracted generating script
# Produces: plan_rebuild-the-psoriasis-pathway-enrichment_afa8da2d.json
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 39a171c0-7492-47dd-9232-6d6eda41fc0c
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

# [lineage] reconstruction failed validation (prose) — raw producing cell below; see Execution Log for full trace
m4 <- data.table::fread(p[["meta_de_PPvsNN.csv"]])
cat("meta_de_PPvsNN.csv cols:", paste(colnames(m4), collapse=", "), "\n")
cat("nrow:", nrow(m4), "\n")
print(m4[gene=="STAT3" | Symbol=="STAT3" | genes=="STAT3", , drop=FALSE][1,])
cat("\n--- STAT3 row (grep any col) ---\n")
print(m4[apply(m4, 1, function(r) any(r=="STAT3", na.rm=TRUE)), ][1,])