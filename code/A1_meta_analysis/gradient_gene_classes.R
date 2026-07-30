# Auto-extracted generating script
# Produces: gradient_gene_classes.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds, gradient_de_genes.rds
# Source artifact version: 666dfb99-0320-46d8-9b7e-f2b9982f55e4
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({library(data.table)})

de <- readRDS("per_study_de.rds")

cls <- readRDS("/tmp/cls_all.rds")
setDT(cls)

mat <- fread(sprintf("clust_input/SRP165679.tsv"))
genes <- mat$Genes
X <- as.matrix(mat[,-1]); rownames(X) <- genes
lab <- cls[srp=="SRP165679"]
keep <- lab[class %in% c("NN","PN","PP")]
X <- X[, keep$external_id]
grp <- factor(keep$class, levels=c("NN","PN","PP"))
stage <- as.integer(grp)-1L

suppressPackageStartupMessages({library(limma)})
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

ppnn <- as.data.table(de$PPvsNN[["SRP165679"]])
de_genes <- ppnn[adjP<0.05 & abs(logFC)>1, gene]

g <- trend[gene %in% de_genes]

Xg <- X[g$gene, ]
rho <- apply(Xg, 1, function(v) suppressWarnings(cor(v, stage, method="spearman")))
g[, spearman_rho := round(rho[gene],3)]

g[, cls := fifelse(pn_frac < 0, "PN_divergent",
            fifelse(pn_frac >= 0.50, "early_PN",
            fifelse(pn_frac >= 0.15, "progressive", "late_PP")))]

fwrite(g[order(cls, -abs(trend_slope)),
         .(gene, cls, trend_slope, trend_FDR, spearman_rho, mu_NN, mu_PN, mu_PP, pn_frac)],
       "gradient_gene_classes.csv")