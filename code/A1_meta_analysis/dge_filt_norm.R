# Auto-extracted generating script
# Produces: dge_filt_norm.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 7fe8fcf7-90ca-4a27-abd2-7fd5bd62810a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({ library(recount3); library(edgeR); library(dplyr) })

human_projects <- available_projects()
proj_info <- subset(human_projects, project == "SRP035988" & project_type == "data_sources")
rse_gene <- create_rse(proj_info)
assay(rse_gene, "counts") <- transform_counts(rse_gene)

colData(rse_gene)$group <- NA
colData(rse_gene)$group[grepl("normal_skin",   colData(rse_gene)$sra.sample_attributes, ignore.case=TRUE)] <- "NN"
colData(rse_gene)$group[grepl("Psoriasis_skin",colData(rse_gene)$sra.sample_attributes, ignore.case=TRUE)] <- "PP"
grp <- factor(colData(rse_gene)$group, levels=c("NN","PP"))

myCounts <- assay(rse_gene, "counts")
gs <- rowData(rse_gene)$gene_name
ok <- !is.na(gs) & gs != ""
myCounts <- rowsum(myCounts[ok,], group = gs[ok])

dge <- DGEList(counts = round(myCounts), group = grp)
keep_fbe <- filterByExpr(dge, group = grp)
dge.filt <- dge[keep_fbe, , keep.lib.sizes = FALSE]
dge.filt.norm <- calcNormFactors(dge.filt, method = "TMM")
saveRDS(dge.filt.norm, "dge_filt_norm.rds")