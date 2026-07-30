# Auto-extracted generating script
# Produces: deconv_trend.csv
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): bulk_cpm_for_deconv.rds
# Source artifact version: dadc31cd-b197-4d4a-885f-b20a6ce686e3
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(nnls)})

sig_data <- readRDS("scissor_repo/results/deconv_signature.rds")
bulk <- readRDS("bulk_cpm_for_deconv.rds")

S <- sig_data$signature
B <- bulk$cpm
common <- intersect(rownames(S), rownames(B))
S <- S[common,]; B <- B[common,]

props <- t(apply(B, 2, function(b){
  w <- nnls(S, b)$x
  w/sum(w)
}))
colnames(props) <- colnames(S)
props <- as.data.frame(props)
props$tier <- bulk$cls
props$sample <- bulk$samples

props$tord <- c(NN=0,PN=1,PP=2)[props$tier]

trend <- data.frame()
for (ct in colnames(S)){
  fit <- lm(props[[ct]] ~ props$tord)
  co <- summary(fit)$coefficients[2,]
  trend <- rbind(trend, data.frame(celltype=ct, slope_per_tier=co[1], p=co[4]))
}
trend$padj <- p.adjust(trend$p, "BH")
trend <- trend[order(trend$p),]

write.csv(trend, "deconv_trend.csv", row.names=FALSE)