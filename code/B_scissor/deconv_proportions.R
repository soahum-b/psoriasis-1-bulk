# Auto-extracted generating script
# Produces: deconv_proportions.csv
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): bulk_cpm_for_deconv.rds
# Source artifact version: 3612ff33-d9e1-4c1e-a9f4-3eb2409eb32f
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

write.csv(props, "deconv_proportions.csv", row.names=FALSE)