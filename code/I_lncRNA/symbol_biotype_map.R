# Auto-extracted generating script
# Produces: symbol_biotype_map.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 29dba5e9-ff9d-4939-8814-6e53916d8d5d
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages(library(recount3))
library(data.table)
options(recount3_url = "http://duffel.rail.bio/recount3")
hp <- available_projects()
proj <- subset(hp, project == "SRP035988" & file_source == "sra")
rse <- create_rse(proj)
rr <- as.data.frame(SummarizedExperiment::rowData(rse))
ann <- rr[, intersect(c("gene_id","gene_name","gene_type","bp_length"), colnames(rr))]
write.csv(ann, "gencode_v26_annotation.csv", row.names=FALSE)
ann <- read.csv("gencode_v26_annotation.csv", stringsAsFactors=FALSE)

lnc_types <- c("lincRNA","antisense","sense_intronic","sense_overlapping",
               "processed_transcript","3prime_overlapping_ncRNA",
               "bidirectional_promoter_lncRNA","macro_lncRNA","non_coding")
ann$biotype_class <- ifelse(ann$gene_type %in% lnc_types, "lncRNA",
                     ifelse(ann$gene_type == "protein_coding", "protein_coding",
                            "other"))
dt <- as.data.table(ann)
sym_map <- dt[gene_name != "", .(
    n_ids = .N,
    any_pc  = any(biotype_class=="protein_coding"),
    any_lnc = any(biotype_class=="lncRNA"),
    gene_types = paste(sort(unique(gene_type)), collapse=";")
), by = gene_name]
sym_map[, biotype := fifelse(any_pc, "protein_coding",
                     fifelse(any_lnc, "lncRNA", "other"))]
fwrite(sym_map, "symbol_biotype_map.csv")