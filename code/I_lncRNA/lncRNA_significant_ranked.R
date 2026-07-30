# Auto-extracted generating script
# Produces: lncRNA_significant_ranked.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de_k4_primary.rds, meta_de_PPvsNN_4study.csv
# Source artifact version: 76020b2b-62fa-4ab0-85ff-7a93eedc2f9b
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

lnc <- de_ann[is_lncRNA==TRUE & !is.na(FDR)]

lnc[, I2_flag := fifelse(I2>=75,"high", fifelse(I2>=50,"moderate","low"))]
lnc[, k_flag  := fifelse(k>=4,"all4", fifelse(k==3,"k3","k2_or_less"))]
lnc[, confidence := fifelse(FDR<0.05 & k>=3 & I2<50, "A_high",
                    fifelse(FDR<0.05 & k>=3 & I2<75, "B_moderate",
                    fifelse(FDR<0.05, "C_significant_but_flagged", "NS")))]
lnc[, direction := fifelse(logFC>0,"up_in_PP","down_in_PP")]
lnc[, absFC := abs(logFC)]

sig <- lnc[FDR<0.05]
setorder(sig, confidence, -absFC)
sig[, rank := .I]
cols <- c("rank","gene","gene_types","logFC","direction","SE","P","FDR","FDR_lncRNA","I2","I2_flag","k","k_flag","confidence")
fwrite(sig[, ..cols], "lncRNA_significant_ranked.csv")