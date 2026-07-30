# Auto-extracted generating script
# Produces: scissor_tuning.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 2f894659-849e-4dfa-8961-9068c0dbb96a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(Seurat); library(Matrix); library(glmnet)})

# scissor_glmnet_solver.R contents
suppressMessages(library(glmnet))
suppressMessages(library(Matrix))

build_edge_aug <- function(A) {
  A <- as(A, "CsparseMatrix")
  diag(A) <- 0
  A@x[A@x != 0] <- 1
  d <- Matrix::rowSums(A); d[d == 0] <- 1
  At <- as(triu(A), "TsparseMatrix")
  ei <- At@i + 1L; ej <- At@j + 1L
  ne <- length(ei); p <- nrow(A)
  B <- sparseMatrix(
    i = rep(seq_len(ne), 2L),
    j = c(ei, ej),
    x = c(1/sqrt(d[ei]), -1/sqrt(d[ej])),
    dims = c(ne, p)
  )
  list(B = B, n_edge = ne, degree = d)
}

net_enet_fit <- function(X, y, B, alpha, nfolds = 10, seed = 123) {
  set.seed(seed)
  N <- nrow(X); p <- ncol(X)
  lambda2 <- (1 - alpha)
  Xs <- as(as.matrix(X), "CsparseMatrix")
  Xa <- rbind(Xs, sqrt(lambda2) * B)
  ya <- c(y, rep(0, nrow(B)))
  cvf <- cv.glmnet(Xa, ya, alpha = 1, nfolds = nfolds, standardize = FALSE,
                   intercept = TRUE)
  beta <- as.numeric(coef(cvf, s = "lambda.min"))[-1]
  names(beta) <- colnames(X)
  list(beta = beta, lambda.min = cvf$lambda.min, cvf = cvf)
}

net_enet_path <- function(X, y, B, alpha, target_frac = 0.20,
                          nlambda = 100, seed = 123) {
  set.seed(seed)
  p <- ncol(X)
  lambda2 <- (1 - alpha)
  Xs <- as(as.matrix(X), "CsparseMatrix")
  Xa <- rbind(Xs, sqrt(lambda2) * B)
  ya <- c(y, rep(0, nrow(B)))
  gf <- glmnet(Xa, ya, alpha = 1, nlambda = nlambda, standardize = FALSE,
               intercept = TRUE)
  Bmat <- as.matrix(gf$beta)
  fracs <- colMeans(Bmat != 0)
  ok <- which(fracs <= target_frac & fracs > 0)
  idx <- if (length(ok)) ok[which.max(fracs[ok])] else which.min(abs(fracs - target_frac))
  beta <- Bmat[, idx]; names(beta) <- colnames(X)
  list(beta = beta, lambda = gf$lambda[idx], frac = fracs[idx],
       path = data.frame(lambda = gf$lambda, frac = fracs), idx = idx)
}

run_scissor <- function(inp, alpha_grid = c(0.05,0.1,0.2,0.3,0.4,0.5),
                        cutoff = 0.20, target_frac = 0.15, seed = 123,
                        verbose = TRUE) {
  X <- inp$X; y <- inp$y; B <- inp$B
  tuning <- data.frame(); fits <- list()
  for (a in alpha_grid) {
    fit <- net_enet_path(X, y, B, alpha = a, target_frac = target_frac, seed = seed)
    beta <- fit$beta
    pos <- names(beta)[beta > 0]; neg <- names(beta)[beta < 0]
    frac <- (length(pos) + length(neg)) / ncol(X)
    tuning <- rbind(tuning, data.frame(alpha = a, n_pos = length(pos),
                    n_neg = length(neg), selected_frac = frac,
                    lambda = fit$lambda))
    fits[[as.character(a)]] <- list(alpha = a, beta = beta, pos = pos, neg = neg,
                                    frac = frac, lambda = fit$lambda, path = fit$path)
    if (verbose) cat(sprintf("alpha=%.2f  pos=%d neg=%d  frac=%.3f%%  lambda=%.4g\n",
                             a, length(pos), length(neg), frac*100, fit$lambda))
  }
  elig <- tuning[tuning$selected_frac <= cutoff & tuning$selected_frac > 0, ]
  if (nrow(elig) == 0) elig <- tuning
  elig <- elig[order(abs(elig$selected_frac - target_frac), -elig$alpha), ]
  chosen <- fits[[as.character(elig$alpha[1])]]
  list(tuning = tuning, chosen = chosen, fits = fits)
}

inp <- readRDS("scissor_repo/results/scissor_inputs.rds")

t0 <- Sys.time()
res <- run_scissor(inp, alpha_grid=c(0.05,0.1,0.2,0.3,0.4,0.5),
                   cutoff=0.20, target_frac=0.15, seed=123, verbose=TRUE)
cat("done in", round(difftime(Sys.time(),t0,units="mins"),2),"min\n")
ch <- res$chosen
cat(sprintf("\nCHOSEN: alpha=%.2f  frac=%.2f%%  pos=%d  neg=%d  lambda=%.4g\n",
            ch$alpha, ch$frac*100, length(ch$pos), length(ch$neg), ch$lambda))
saveRDS(res, "scissor_tuning.rds")
write.csv(res$tuning, "scissor_repo/results/scissor_alpha_tuning.csv", row.names=FALSE)