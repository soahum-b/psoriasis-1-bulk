# Auto-extracted generating script
# Produces: meta_de_PPvsNN_4study_biotyped.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): meta_de_PPvsNN_4study.csv
# Source artifact version: 8c104c0b-9864-4f53-a73c-60a1c07899be
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(recount3)

options(recount3_url = "http://duffel.rail.bio/recount3")
hp <- available_projects()
proj <- subset(hp, project == "SRP035988" & file_source == "sra")
rse <- create_rse(proj)
rr <- as.data.frame(SummarizedExperiment::rowData(rse))
ann <- rr[, intersect(c("gene_id","gene_name","gene_type","bp_length"), colnames(rr))]

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

de <- fread("meta_de_PPvsNN_4study.csv")

de_ann <- merge(de, sym_map[, .(gene=gene_name, biotype, gene_types)],
                by="gene", all.x=TRUE)
de_ann[is.na(biotype), biotype := "unmapped"]
de_ann[, is_lncRNA := biotype=="lncRNA"]

de_ann[, FDR_lncRNA := NA_real_]
lnc_idx <- de_ann$is_lncRNA & !is.na(de_ann$P)
de_ann[lnc_idx, FDR_lncRNA := p.adjust(P, method="BH")]

fwrite(de_ann, "meta_de_PPvsNN_4study_biotyped.csv")