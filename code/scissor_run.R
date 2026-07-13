# scissor_run.R
# Faithful Scissor driver (Sun & Xia 2022) using the pure-R glmnet network
# solver. Mirrors Scissor::Scissor(): quantile-normalize the combined bulk+cell
# matrix over shared genes, compute the cell-by-bulk Pearson correlation, build
# the cell-cell SNN graph, then solve the graph-regularized elastic net and
# select cells by coefficient sign, walking an alpha grid until the selected
# fraction falls under `cutoff`.

suppressMessages({library(Seurat); library(Matrix); library(preprocessCore)})
source("code/scissor_glmnet_solver.R")

# Prepare Scissor regression inputs (X, y, network) — the Scissor() preamble.
prepare_scissor_inputs <- function(bulk_logcpm, seurat_obj, y,
                                    hvg_only = TRUE, save_file = NULL) {
  sc_data <- GetAssayData(seurat_obj, layer = "data")     # log-normalized
  genes_sc <- rownames(sc_data)
  if (hvg_only) genes_sc <- intersect(VariableFeatures(seurat_obj), rownames(bulk_logcpm))
  common <- intersect(rownames(bulk_logcpm), genes_sc)
  stopifnot(length(common) > 100)

  Bk <- as.matrix(bulk_logcpm[common, ])
  Sc <- as.matrix(sc_data[common, ])
  # quantile-normalize the combined matrix (Scissor step)
  comb <- cbind(Bk, Sc)
  combn <- normalize.quantiles(comb)
  rownames(combn) <- rownames(comb); colnames(combn) <- colnames(comb)
  Eb <- combn[, 1:ncol(Bk)]
  Ec <- combn[, (ncol(Bk)+1):ncol(combn)]

  X <- cor(Eb, Ec)                          # bulk x cells Pearson correlation
  qc <- quantile(X)
  # SNN network from the Seurat graph
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

# Run Scissor selection over an alpha (graph-strength) grid. For each alpha the
# L1 penalty is walked along glmnet's path and set to hit a modest selected
# fraction (target_frac), the operational form of Scissor's cutoff rule. Among
# the alphas we keep the one whose achieved fraction is closest to target while
# <= cutoff, preferring stronger graph smoothing (larger alpha) on ties.
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
  # pick alpha: achieved frac <= cutoff, closest to target, larger alpha wins ties
  elig <- tuning[tuning$selected_frac <= cutoff & tuning$selected_frac > 0, ]
  if (nrow(elig) == 0) elig <- tuning
  elig <- elig[order(abs(elig$selected_frac - target_frac), -elig$alpha), ]
  chosen <- fits[[as.character(elig$alpha[1])]]
  list(tuning = tuning, chosen = chosen, fits = fits)
}
