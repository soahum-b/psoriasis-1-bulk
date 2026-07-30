# Auto-extracted generating script
# Produces: de_results_full.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): dge_filt_norm.rds
# Source artifact version: 40f54a42-9257-485b-9795-a4c4c059d5e9
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR); library(limma); library(qvalue)})

dge <- readRDS("dge_filt_norm.rds")
grp <- dge$samples$group

design <- model.matrix(~ grp)
colnames(design) <- c("Intercept", "PPvsNN")

v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "PPvsNN", number = Inf, sort.by = "P")
res$gene <- rownames(res)

qobj <- qvalue(p = res$P.Value)
res$qvalue <- qobj$qvalues

de <- res[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","qvalue","B")]

saveRDS(de, "de_results_full.rds")