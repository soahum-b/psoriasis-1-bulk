# Auto-extracted generating script
# Produces: scissor_inputs.rds
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): bulk_anchor.rds, reference_subset20k.rds
# Source artifact version: a01cf77a-b6d5-4685-8769-fe9af0f893db
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

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

suppressMessages({library(Seurat); library(Matrix); library(preprocessCore)})

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

sub <- readRDS("reference_subset20k.rds")
anch <- readRDS("bulk_anchor.rds")

t0 <- Sys.time()
inp <- prepare_scissor_inputs(anch$logcpm, sub, y=anch$tier, hvg_only=TRUE,
                              save_file="scissor_inputs.rds")
cat("inputs prepared in", round(difftime(Sys.time(),t0,units="secs"),1),"s\n")
cat("X dim (bulk x cells):", paste(dim(inp$X),collapse=" x "), "\n")
cat("common genes:", length(inp$common), "| edges:", inp$n_edge, "\n")
cat("correlation quality (five-number):\n"); print(round(inp$qc,3))