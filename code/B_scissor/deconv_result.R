# Auto-extracted generating script
# Produces: deconv_result.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): bulk_cpm_for_deconv.rds
# Source artifact version: a702dd82-4cba-46af-964f-075076747a50
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(nnls)})
sig <- readRDS("scissor_repo/results/deconv_signature.rds")
bulk <- readRDS("bulk_cpm_for_deconv.rds")

S <- sig$signature                       # genes x celltypes
B <- bulk$cpm                            # genes x samples
common <- intersect(rownames(S), rownames(B))
cat("common signature genes in bulk:", length(common), "/", nrow(S), "\n")
S <- S[common,]; B <- B[common,]
# log-scale often destabilizes NNLS proportions; use linear CPM, per-column scale
# deconvolve each bulk sample: b ~ S %*% w, w>=0, normalize to fractions
props <- t(apply(B, 2, function(b){
  w <- nnls(S, b)$x
  w/sum(w)
}))
colnames(props) <- colnames(S)
props <- as.data.frame(props)
props$tier <- bulk$cls
props$sample <- bulk$samples

# mean proportion per tier
agg <- aggregate(props[,colnames(S)], list(tier=props$tier), mean)
agg <- agg[match(c("NN","PN","PP"), agg$tier),]
cat("\n=== mean estimated cell-type proportions by tier (%) ===\n")
print(round(agg[,-1]*100,2), row.names=FALSE)
cat("tiers:", agg$tier, "\n")

# test monotonic trend across NN<PN<PP for each cell type (Jonckheere-style: lm on ordinal)
props$tord <- c(NN=0,PN=1,PP=2)[props$tier]
cat("\n=== trend test (proportion ~ ordinal tier) ===\n")
trend <- data.frame()
for (ct in colnames(S)){
  fit <- lm(props[[ct]] ~ props$tord)
  co <- summary(fit)$coefficients[2,]
  trend <- rbind(trend, data.frame(celltype=ct, slope_per_tier=co[1], p=co[4]))
}
trend$padj <- p.adjust(trend$p, "BH")
trend <- trend[order(trend$p),]
print(format(trend, digits=3), row.names=FALSE)
saveRDS(list(props=props, agg=agg, trend=trend), "deconv_result.rds")
write.csv(props, "scissor_repo/results/deconv_proportions.csv", row.names=FALSE)
write.csv(trend, "scissor_repo/results/deconv_trend.csv", row.names=FALSE)
cat("\nsaved deconv_result.rds, deconv_proportions.csv, deconv_trend.csv\n")