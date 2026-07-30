# Auto-extracted generating script
# Produces: permutation_null.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): scissor_result.rds, scissor_inputs.rds, scissor_run.R, scissor_glmnet_solver.R, reference_subset20k.rds
# Source artifact version: 2030cff3-14ca-482a-a00c-9d7c116853b7
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(glmnet); library(Matrix)})
source("scissor_glmnet_solver.R")
source("scissor_run.R")
inp <- readRDS("scissor_inputs.rds")
res <- readRDS("scissor_result.rds")
sub_md <- readRDS("reference_subset20k.rds")@meta.data
alpha <- res$chosen$alpha
real_frac <- res$chosen$frac
tiers <- inp$y

set.seed(202)
n_perm <- 30
null_frac <- numeric(n_perm); null_dir <- numeric(n_perm)
real_pos_tier <- mean(sub_md[res$chosen$pos,"tier"])
real_neg_tier <- mean(sub_md[res$chosen$neg,"tier"])
real_gap <- real_pos_tier - real_neg_tier

t0 <- Sys.time()
for (i in 1:n_perm){
  yp <- tiers[sample(length(tiers))]
  fit <- net_enet_path(inp$X, yp, inp$B, alpha=alpha, target_frac=0.15, seed=123)
  b <- fit$beta
  pos <- names(b)[b>0]; neg <- names(b)[b<0]
  null_frac[i] <- (length(pos)+length(neg))/ncol(inp$X)
  pt <- if(length(pos)) mean(sub_md[pos,"tier"]) else NA
  nt <- if(length(neg)) mean(sub_md[neg,"tier"]) else NA
  null_dir[i] <- pt - nt
  if (i %% 10==0) cat("perm", i, "/", n_perm, "\n")
}
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
cat(sprintf("\nREAL: selected frac=%.2f%%  pos-neg tier gap=%.3f\n", real_frac*100, real_gap))
cat(sprintf("NULL: frac mean=%.2f%% (sd %.2f)  |  tier-gap mean=%.3f (sd %.3f)\n",
            mean(null_frac)*100, sd(null_frac)*100, mean(null_dir,na.rm=T), sd(null_dir,na.rm=T)))
cat(sprintf("p(null tier-gap >= real): %.3f\n", mean(null_dir >= real_gap, na.rm=TRUE)))
saveRDS(list(real_frac=real_frac, real_gap=real_gap, null_frac=null_frac,
             null_dir=null_dir, n_perm=n_perm), "permutation_null.rds")
write.csv(data.frame(perm=1:n_perm, null_selected_frac=null_frac, null_tier_gap=null_dir),
          "scissor_repo/results/permutation_null.csv", row.names=FALSE)
cat("saved permutation_null\n")