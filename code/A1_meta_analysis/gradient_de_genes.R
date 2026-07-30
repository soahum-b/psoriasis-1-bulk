# Auto-extracted generating script
# Produces: gradient_de_genes.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds, sample_classification.rds
# Source artifact version: a27965d8-adbb-4661-9b0b-bcd01c3bfb76
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({library(limma); library(data.table)})

cls <- readRDS("sample_classification.rds")
cls <- as.data.table(cls)

de <- readRDS("per_study_de.rds")

# Build expression matrix + group vector for SRP165679 (3 groups only)
mat <- fread("clust_input/SRP165679.tsv")
genes <- mat$Genes
X <- as.matrix(mat[,-1]); rownames(X) <- genes
lab <- cls[srp=="SRP165679"]
keep <- lab[class %in% c("NN","PN","PP")]
X <- X[, keep$external_id]
grp <- factor(keep$class, levels=c("NN","PN","PP"))
stage <- as.integer(grp)-1L   # NN=0, PN=1, PP=2

# 1) LINEAR TREND via limma: design with numeric stage covariate
des <- model.matrix(~ stage)
fit <- lmFit(X, des); fit <- eBayes(fit, trend=TRUE)
tt <- topTable(fit, coef="stage", number=Inf, sort.by="none")
trend <- data.table(gene=rownames(X), trend_slope=tt$logFC, trend_t=tt$t,
                    trend_P=tt$P.Value, trend_FDR=tt$adj.P.Val)

# 2) Group means for monotonicity classification
mu <- sapply(c("NN","PN","PP"), function(g) rowMeans(X[, grp==g, drop=FALSE]))
colnames(mu) <- c("mu_NN","mu_PN","mu_PP")
trend <- cbind(trend, mu)
trend[, mono_up   := mu_NN < mu_PN & mu_PN < mu_PP]
trend[, mono_down := mu_NN > mu_PN & mu_PN > mu_PP]
trend[, monotonic := mono_up | mono_down]
trend[, pn_frac := round((mu_PN-mu_NN)/(mu_PP-mu_NN), 3)]

# Restrict to lesional DE genes
ppnn <- as.data.table(de$PPvsNN[["SRP165679"]])
de_genes <- ppnn[adjP<0.05 & abs(logFC)>1, gene]

g <- trend[gene %in% de_genes]

# Spearman rho of expression vs stage per gene
lab2 <- cls[srp=="SRP165679"][class %in% c("NN","PN","PP")]
Xg <- X[g$gene, lab2$external_id]
stage2 <- as.integer(factor(lab2$class, levels=c("NN","PN","PP")))-1L
rho <- apply(Xg, 1, function(v) suppressWarnings(cor(v, stage2, method="spearman")))
g[, spearman_rho := round(rho[gene], 3)]

saveRDS(g, "gradient_de_genes.rds")