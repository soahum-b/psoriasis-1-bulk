# Auto-extracted generating script
# Produces: reliability_test.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): scissor_result.rds, scissor_inputs.rds, scissor_glmnet_solver.R, reference_subset20k.rds
# Source artifact version: 36776557-a136-4970-adf5-243164d1151d
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(glmnet); library(Matrix)})

# scissor_glmnet_solver.R (inline from artifact)
source("scissor_glmnet_solver.R")

# scissor_reliability.R (inline)
.cv_mse_lm <- function(X, y, B, alpha, cell_num, nfold, foldid) {
  lambda2 <- (1 - alpha)
  mse <- numeric(nfold)
  for (j in seq_len(nfold)) {
    te <- which(foldid == j); tr <- setdiff(seq_len(nrow(X)), te)
    Xtr <- as(as.matrix(X[tr, , drop=FALSE]), "CsparseMatrix")
    Xa <- rbind(Xtr, sqrt(lambda2) * B)
    ya <- c(y[tr], rep(0, nrow(B)))
    gf <- glmnet(Xa, ya, alpha = 1, nlambda = 100, standardize = FALSE, intercept = TRUE)
    nz <- gf$df
    idx <- which.min(abs(nz - cell_num))
    beta <- gf$beta[, idx]; a0 <- gf$a0[idx]
    pred <- as.numeric(X[te, , drop=FALSE] %*% beta) + a0
    mse[j] <- mean((y[te] - pred)^2)
  }
  mse
}

reliability_test_glmnet <- function(inp, alpha, cell_num, n = 100, nfold = 10,
                                     seed = 1, verbose = TRUE) {
  X <- inp$X; y <- inp$y; B <- inp$B
  m <- nrow(X)
  set.seed(seed)
  foldid <- sample(cut(seq(m), breaks = nfold, labels = FALSE))

  real_mse <- .cv_mse_lm(X, y, B, alpha, cell_num, nfold, foldid)
  statistic <- mean(real_mse)
  if (verbose) cat(sprintf("real mean CV-MSE = %.4f\n", statistic))

  back <- numeric(n)
  for (i in seq_len(n)) {
    set.seed(i + 100)
    yp <- y[sample(m)]
    back[i] <- mean(.cv_mse_lm(X, yp, B, alpha, cell_num, nfold, foldid))
    if (verbose && i %% 20 == 0) cat(sprintf("  perm %d/%d  mean-null-MSE so far=%.4f\n",
                                             i, n, mean(back[1:i])))
  }
  p <- sum(back < statistic) / n
  list(statistic = statistic, p = p, real_mse = real_mse, background = back,
       alpha = alpha, cell_num = cell_num, n = n)
}

inp <- readRDS("scissor_inputs.rds")
res <- readRDS("scissor_result.rds")
cell_num <- length(res$chosen$pos) + length(res$chosen$neg)
alpha <- res$chosen$alpha
t0 <- Sys.time()
rt <- reliability_test_glmnet(inp, alpha=alpha, cell_num=cell_num, n=100, nfold=10, verbose=TRUE)
cat(sprintf("\nDONE in %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
cat(sprintf("statistic(real MSE)=%.4f  null mean=%.4f  p=%.3f\n",
            rt$statistic, mean(rt$background), rt$p))
saveRDS(rt, "reliability_test.rds")
cat("saved reliability_test.rds\n")