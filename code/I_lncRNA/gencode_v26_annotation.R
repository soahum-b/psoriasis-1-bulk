# Auto-extracted generating script
# Produces: gencode_v26_annotation.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: d1260db6-d2f1-4f52-ae7f-49a06bb05132
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages(library(recount3))
options(recount3_url = "http://duffel.rail.bio/recount3")
hp <- available_projects()
proj <- subset(hp, project == "SRP035988" & file_source == "sra")
rse <- create_rse(proj)
rr <- as.data.frame(SummarizedExperiment::rowData(rse))
ann <- rr[, intersect(c("gene_id","gene_name","gene_type","bp_length"), colnames(rr))]
write.csv(ann, "gencode_v26_annotation.csv", row.names=FALSE)