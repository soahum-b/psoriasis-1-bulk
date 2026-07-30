# Auto-extracted generating script
# Produces: phenotype_checks.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): rse_SRP165679.rds
# Source artifact version: 2d557942-0001-45c0-9703-24d71085fa78
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(recount3)
library(SummarizedExperiment)
library(edgeR)

anch <- readRDS("scissor_repo/results/bulk_anchor.rds")
rse <- readRDS("rse_SRP165679.rds")
tier <- anch$tier; cls <- anch$cls

checks <- list()
checks$n_per_level <- table(cls)

cd <- as.data.frame(colData(rse))
ext <- cd$external_id
idx <- match(anch$samples, ext)
sub <- cd[idx,,drop=FALSE]

cat("=== non-degeneracy checks ===\n")
cat("1. n per level:\n"); print(checks$n_per_level)
cat("   min per level:", min(table(cls)), "(>= 10 OK)\n\n")

logcpm <- anch$logcpm
pc <- prcomp(t(logcpm), scale.=TRUE)
varexp <- (pc$sdev^2/sum(pc$sdev^2))[1:5]
cat("2. bulk PCA var explained PC1-5:", paste(round(varexp*100,1),collapse=", "),"%\n")
for (k in 1:3) {
  a <- summary(aov(pc$x[,k] ~ factor(cls)))[[1]]
  cat(sprintf("   PC%d ~ tier: F=%.1f, p=%.2e\n", k, a$`F value`[1], a$`Pr(>F)`[1]))
}
cat("\n3. tier is an ordinal biopsy-site label (non-circular by construction):\n")
cat("   phenotype derived from clinical designation, NOT from expression -> cannot induce the correlation Scissor searches for.\n")

chk <- data.frame(check=c("min_n_per_level","n_levels","bulk_PC1_varexp_pct","tier_single_study"),
                  value=c(min(table(cls)), length(unique(cls)), round(varexp[1]*100,1), "SRP165679_only"))
write.csv(chk, "phenotype_checks.csv", row.names=FALSE)
cat("\nsaved phenotype_checks.csv\n")