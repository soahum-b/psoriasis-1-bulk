# Auto-extracted generating script
# Produces: scissor_result.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): reference_subset20k.rds
# Source artifact version: 5951a328-abda-43dd-9fab-3947e289b53f
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(Seurat)
library(Matrix)
library(glmnet)
library(preprocessCore)

# scissor_glmnet_solver.R content
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

# scissor_run.R content
prepare_scissor_inputs <- function(bulk_logcpm, seurat_obj, y,
                                    hvg_only = TRUE, save_file = NULL) {
  sc_data <- GetAssayData(seurat_obj, layer = "data")
  genes_sc <- rownames(sc_data)
  if (hvg_only) genes_sc <- intersect(VariableFeatures(seurat_obj), rownames(bulk_logcpm))
  common <- intersect(rownames(bulk_logcpm), genes_sc)
  stopifnot(length(common) > 100)

  Bk <- as.matrix(bulk_logcpm[common, ])
  Sc <- as.matrix(sc_data[common, ])
  comb <- cbind(Bk, Sc)
  combn <- normalize.quantiles(comb)
  rownames(combn) <- rownames(comb); colnames(combn) <- colnames(comb)
  Eb <- combn[, 1:ncol(Bk)]
  Ec <- combn[, (ncol(Bk)+1):ncol(combn)]

  X <- cor(Eb, Ec)
  qc <- quantile(X)
  gname <- grep("_snn$", names(seurat_obj@graphs), value = TRUE)[1]
  if (is.na(gname)) gname <- names(seurat_obj@graphs)[1]
  net <- as(seurat_obj@graphs[[gname]], "CsparseMatrix")
  net <- net[colnames(Ec), colnames(Ec)]
  aug <- build_edge_aug(net)

  out <- list(X = X, y = as.numeric(y), B = aug$B, n_edge = aug$n_edge,
              common = common, qc = qc, cells = colnames(Ec))
  if (!is.null(save_file)) saveRDS(out, save_file)
  out
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

# Load data
sub <- readRDS("reference_subset20k.rds")
anch <- readRDS("scissor_repo/results/bulk_anchor.rds")

# Prepare inputs
inp <- prepare_scissor_inputs(anch$logcpm, sub, y=anch$tier, hvg_only=TRUE,
                              save_file="scissor_repo/results/scissor_inputs.rds")

# Run scissor
res <- run_scissor(inp, alpha_grid=c(0.05,0.1,0.2,0.3,0.4,0.5),
                   cutoff=0.20, target_frac=0.15, seed=123, verbose=TRUE)

ch <- res$chosen
tun <- res$tuning

# Save scissor result
md <- sub@meta.data
scissor_label <- rep("Background", nrow(md)); names(scissor_label) <- rownames(md)
scissor_label[ch$pos] <- "Scissor+"; scissor_label[ch$neg] <- "Scissor-"
coefs <- ch$beta
saveRDS(list(chosen=ch, tuning=tun, scissor_label=scissor_label, coefs=coefs,
             cells=names(scissor_label)), "scissor_result.rds")