# Auto-extracted generating script
# Produces: gradient_genes.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds, sample_classification.rds
# Source artifact version: 447dba47-a4a6-4cba-9aad-a91048061d7e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(limma)})

de <- readRDS("per_study_de.rds")
cls <- as.data.table(readRDS("sample_classification.rds"))

mat_file <- "clust_input/SRP165679.tsv"
mat <- fread(mat_file)
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

ppnn <- as.data.table(de$PPvsNN[["SRP165679"]])
de_genes <- ppnn[adjP<0.05 & abs(logFC)>1, gene]

g <- trend[gene %in% de_genes]

Xg <- X[g$gene, ]
stage_vec <- stage
rho <- apply(Xg, 1, function(v) suppressWarnings(cor(v, stage_vec, method="spearman")))
g[, spearman_rho := round(rho[gene],3)]

ppnn2 <- as.data.table(de$PPvsNN[["SRP165679"]])[,.(gene, les_logFC=logFC, les_FDR=adjP)]
out <- merge(g, ppnn2, by="gene")
out[, monotonic_flag := ifelse(mono_up,"up_NN<PN<PP", ifelse(mono_down,"down_NN>PN>PP","non_monotonic"))]
setcolorder(out, c("gene","les_logFC","les_FDR","trend_slope","trend_t","trend_FDR",
                   "spearman_rho","mu_NN","mu_PN","mu_PP","pn_frac","monotonic_flag"))
out <- out[order(-abs(trend_slope))]
fwrite(out[,.(gene,les_logFC,les_FDR,trend_slope,trend_t,trend_FDR,spearman_rho,
              mu_NN,mu_PN,mu_PP,pn_frac,monotonic_flag)],
       "gradient_genes.csv")