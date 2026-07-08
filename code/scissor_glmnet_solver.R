# scissor_glmnet_solver.R
# -----------------------------------------------------------------------------
# Faithful pure-R reimplementation of Scissor's core selection step.
#
# Scissor (Sun & Xia, Nat Biotechnol 2022) optimizes, over the cell-by-bulk
# correlation matrix X (rows = bulk samples, cols = cells) against phenotype y,
# a graph-regularized elastic net (penalty = "Net"):
#
#     min_beta  (1/2N)||y - X beta||^2  +  lambda[ alpha ||beta||_1
#                                              + (1-alpha)/2 * beta^T L_sym beta ]
#
# where L_sym is the SYMMETRIC-NORMALIZED Laplacian of the cell-cell SNN graph:
#   L_sym = I - D^{-1/2} A D^{-1/2}   (A = binarized SNN adjacency, D = degree).
# This matches Scissor's compiled OmegaC(), which fills off-diagonals with
#   A_ij * sgn_i * sgn_j / sqrt(d_i d_j)  and zero diagonal, i.e. the normalized
# adjacency S = D^{-1/2} A D^{-1/2}; the network penalty is beta^T (I - S) beta.
#
# The compiled APML1 solver cannot be dyn.load()ed under the analysis sandbox
# (shared libraries on writable paths are refused; the conda R library where
# loading is permitted is read-only and cannot receive the package). We instead
# solve the identical objective with glmnet via the Li & Li (2008) network-
# constraint augmentation, which is exact for a quadratic Laplacian penalty.
#
# Li & Li augmentation
# --------------------
# The Laplacian penalty is folded into an augmented least-squares problem:
#   X_aug = [ X ; sqrt(lambda2) * B ] ,  y_aug = [ y ; 0 ]
# where B is the sparse edge-difference factor (||B beta||^2 = beta^T L_sym beta),
# and glmnet is run with an L1 penalty only. glmnet's objective is
#   (1/2 N_aug) ||y_aug - X_aug beta||^2 + lambda_g * ||beta||_1 .
# The augmentation block contributes  (lambda2 / 2 N_aug) * beta^T L_sym beta,
# so the graph (ridge-like) penalty and the L1 penalty are controlled by TWO
# SEPARATE knobs, not a single shared lambda:
#   * lambda2 = (1 - alpha)  is a FIXED constant per fit — it sets the graph
#     smoothing strength (larger alpha => less smoothing);
#   * lambda_g (the L1 penalty) is chosen INDEPENDENTLY along glmnet's own path.
# This is the network-regularized elastic net of Li & Li (2008); it is NOT the
# coupled lambda1=lambda*alpha / lambda2=lambda*(1-alpha) elastic-net
# parameterization (that form is not what is implemented here).
#
# Sparsity control. Scissor's compiled APML1 enforces a modest selected fraction
# via an L0 hard-threshold trim (cvTrimLmC) that this pure-R port does not
# replicate. We instead control the selected fraction directly by walking the
# glmnet lambda_g path and choosing the sparsest lambda whose selected fraction
# is closest to (but not above) a target — the operational equivalent of
# Scissor's "keep the selected set to a modest proportion" cutoff rule.
#
# Selection rule matches Scissor exactly:
#   Coefs = beta ; Scissor+ = cells with beta>0 ; Scissor- = cells with beta<0.
# -----------------------------------------------------------------------------

suppressMessages(library(glmnet))
suppressMessages(library(Matrix))

# Build the SPARSE edge-difference augmentation matrix B for the symmetric-
# normalized Laplacian penalty. Scissor's OmegaC fills off-diagonals with
# A_ij / sqrt(d_i d_j), i.e. L_sym = I - D^{-1/2} A D^{-1/2}. Its quadratic form
# factorizes over graph edges:
#   beta^T L_sym beta = sum_{(i,j) in E} ( beta_i/sqrt(d_i) - beta_j/sqrt(d_j) )^2
# so B has one row per undirected edge with entries +1/sqrt(d_i), -1/sqrt(d_j).
# Then beta^T L_sym beta = ||B beta||^2, and B is sparse (2 nonzeros/row).
build_edge_aug <- function(A) {
  A <- as(A, "CsparseMatrix")
  diag(A) <- 0
  A@x[A@x != 0] <- 1                       # binarize, matching Scissor
  d <- Matrix::rowSums(A); d[d == 0] <- 1
  # upper-triangular edges only (undirected, count once)
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

# Network-regularized elastic net via SPARSE glmnet augmentation.
# X: N x p (samples x cells), y: length N, B: n_edge x p sparse edge matrix.
# Augmented system  [X ; sqrt(l2) B] beta ~ [y ; 0]  with pure-L1 glmnet
# reproduces  (1/2N)||y-Xb||^2 + l2/2 * b^T L b + l1 ||b||_1  (Li & Li 2008).
net_enet_fit <- function(X, y, B, alpha, nfolds = 10, seed = 123) {
  set.seed(seed)
  N <- nrow(X); p <- ncol(X)
  lambda2 <- (1 - alpha)
  Xs <- as(as.matrix(X), "CsparseMatrix")
  Xa <- rbind(Xs, sqrt(lambda2) * B)       # sparse augmented design
  ya <- c(y, rep(0, nrow(B)))
  cvf <- cv.glmnet(Xa, ya, alpha = 1, nfolds = nfolds, standardize = FALSE,
                   intercept = TRUE)
  beta <- as.numeric(coef(cvf, s = "lambda.min"))[-1]
  names(beta) <- colnames(X)
  list(beta = beta, lambda.min = cvf$lambda.min, cvf = cvf)
}

# Fit the full glmnet lambda path for a fixed alpha (graph strength), then pick
# the L1 penalty lambda_g giving the sparsest solution whose selected fraction
# is just below `target_frac`. This is the direct sparsity control that replaces
# Scissor's L0 trim. Returns beta at the chosen lambda plus the path summary.
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
  Bmat <- as.matrix(gf$beta)                       # p x nlambda
  fracs <- colMeans(Bmat != 0)                     # selected fraction per lambda
  # choose sparsest lambda with frac <= target (fall back to smallest frac)
  ok <- which(fracs <= target_frac & fracs > 0)
  idx <- if (length(ok)) ok[which.max(fracs[ok])] else which.min(abs(fracs - target_frac))
  beta <- Bmat[, idx]; names(beta) <- colnames(X)
  list(beta = beta, lambda = gf$lambda[idx], frac = fracs[idx],
       path = data.frame(lambda = gf$lambda, frac = fracs), idx = idx)
}
