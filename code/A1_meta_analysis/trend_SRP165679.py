# Auto-extracted generating script
# Produces: trend_SRP165679.rds
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): sample_classification.rds
# Source artifact version: d61263b9-9355-4e92-97d8-4987bcb38e86
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(limma)})
setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")

cls <- readRDS("sample_classification.rds")
cls <- as.data.table(cls)

mat <- fread("clust_input/SRP165679.tsv")
genes <- mat$Genes
X <- as.matrix(mat[,-1]); rownames(X) <- genes
lab <- cls[srp=="SRP165679"]
keep <- lab[class %in% c("NN","PN","PP")]
X <- X[, keep$external_id]
grp <- factor(keep$class, levels=c("NN","PN","PP"))
stage <- as.integer(grp)-1L

des <- model.matrix(~ stage)
fit <- lmFit(X, des); fit <- eBayes(fit, trend=TRUE)
tt <- topTable(fit, coef="stage", number=Inf, sort.by="none")
trend <- data.table(gene=rownames(X), trend_slope=tt$logFC, trend_t=tt$t,
                    trend_P=tt$P.Value, trend_FDR=tt$adj.P.Val)

mu <- sapply(c("NN","PN","PP"), function(g) rowMeans(X[, grp==g, drop=FALSE]))
colnames(mu) <- c("mu_NN","mu_PN","mu_PP")
trend <- cbind(trend, mu)
trend[, mono_up   := mu_NN < mu_PN & mu_PN < mu_PP]
trend[, mono_down := mu_NN > mu_PN & mu_PN > mu_PP]
trend[, monotonic := mono_up | mono_down]
trend[, pn_frac := round((mu_PN-mu_NN)/(mu_PP-mu_NN), 3)]

saveRDS(trend, "trend_SRP165679.rds")