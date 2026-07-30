# Auto-extracted generating script
# Produces: de_results_significant.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 964d448f-8689-4661-9550-fbfd30c84570
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(edgeR)
  library(dplyr)
  library(limma)
  library(qvalue)
})

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

grp <- dge.filt.norm$samples$group

design <- model.matrix(~ grp)
colnames(design) <- c("Intercept", "PPvsNN")

v <- voom(dge.filt.norm, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "PPvsNN", number = Inf, sort.by = "P")
res$gene <- rownames(res)

qobj <- qvalue(p = res$P.Value)
res$qvalue <- qobj$qvalues

de <- res[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","qvalue","B")]

de_sig <- de[de$adj.P.Val < 0.05 & abs(de$logFC) > 1, ]
de_sig <- de_sig[order(de_sig$adj.P.Val), ]
write.csv(de_sig, "de_results_significant.csv", row.names = FALSE)