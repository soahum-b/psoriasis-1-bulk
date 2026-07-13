# scissor_reliability.R
# Faithful port of Scissor::reliability.test() / test_lm() (Gaussian) using the
# pure-R glmnet network solver. Cross-validated MSE of predicting the held-out
# bulk phenotype from the cell coefficients, real label vs n permuted labels.
# p = fraction of permutations whose mean CV-MSE is below the real mean CV-MSE.
# A reliable association => real MSE low, few permutations beat it => small p.

suppressMessages({library(glmnet); library(Matrix)})
source("code/scissor_glmnet_solver.R")

# One CV pass: returns per-fold test MSE. Coefficients are taken at the lambda
# on the path whose nzero is closest to `cell_num` (mirrors test_lm's index rule).
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
