#!/usr/bin/env Rscript
# run_full_census_cluster.R
# ---------------------------------------------------------------------------
# FULL-CENSUS Scissor-on-gradient run — the cluster scale-up of the local
# backbone. The local backbone (this repo's results/) ran on a stratified
# 20,023-cell subset to keep glmnet tractable on a laptop. This script runs the
# IDENTICAL pipeline on all ~89,058 QC-passed cells. It is written for a SLURM
# node with >=128 GB RAM; the glmnet augmented design at full size is the
# memory driver (89k predictors x ~3.1M SNN edges).
#
# Submit example (SLURM):
#   sbatch code/run_full_census.sbatch
#
# Everything is parameterized at the top. Set N_CELLS=Inf for the full census,
# or a number to reproduce/enlarge the subset. Outputs go to results_full/.
#
# CHECKPOINT/RESUME: every expensive stage (Scissor inputs, alpha tuning,
# reliability test, permutation null) writes its result to results_full/ and
# is SKIPPED on a rerun if that file already exists. If the job is killed
# (wall-time, node failure, `scontrol requeue`), just resubmit the same
# script — it picks up at the first incomplete stage instead of restarting
# from scratch. To force a clean rerun of a stage, delete its output file.
# ---------------------------------------------------------------------------

suppressMessages({library(Seurat); library(Matrix); library(glmnet);
                  library(preprocessCore)})
source("code/scissor_glmnet_solver.R")
source("code/scissor_run.R")
source("code/scissor_reliability.R")

## ---- parameters ----
N_CELLS      <- Inf          # Inf = full census; else stratified subset size
ALPHA_GRID   <- c(0.05,0.1,0.2,0.3,0.4,0.5)
CUTOFF       <- 0.20         # max selected fraction
TARGET_FRAC  <- 0.15         # target selected fraction (Scissor cutoff analogue)
N_PERM_REL   <- 100          # reliability-test permutations
N_PERM_NULL  <- 100          # selection-null permutations
SEED         <- 123
OUT          <- "results_full"
dir.create(OUT, showWarnings = FALSE)

## ---- load inputs ----
message("Loading processed reference + bulk anchor ...")
so   <- readRDS("results/reference_processed.rds")   # 89,058 cells
anch <- readRDS("results/bulk_anchor.rds")

## ---- optional stratified subsample (skip when N_CELLS = Inf) ----
if (is.finite(N_CELLS) && N_CELLS < ncol(so)) {
  set.seed(42)
  md <- so@meta.data; md$strata <- paste(md$celltype, md$condition, sep="|")
  fr <- N_CELLS / ncol(so)
  pick <- unlist(lapply(split(rownames(md), md$strata), function(ids)
            sample(ids, min(length(ids), max(20, round(length(ids)*fr))))))
  so <- subset(so, cells = pick)
}
message("Cells in run: ", ncol(so))
so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)

## ---- Scissor inputs (checkpointed) ----
f_inp <- file.path(OUT, "scissor_inputs.rds")
if (file.exists(f_inp)) {
  message("[checkpoint] loading existing ", f_inp)
  inp <- readRDS(f_inp)
} else {
  message("Preparing Scissor inputs (correlation + SNN edge augmentation) ...")
  inp <- prepare_scissor_inputs(anch$logcpm, so, y = anch$tier, hvg_only = TRUE,
                                save_file = f_inp)
}

## ---- alpha tuning (checkpointed) ----
f_tune <- file.path(OUT, "scissor_tuning.rds")
if (file.exists(f_tune)) {
  message("[checkpoint] loading existing ", f_tune)
  res <- readRDS(f_tune)
} else {
  message("Running alpha tuning ...")
  res <- run_scissor(inp, alpha_grid = ALPHA_GRID, cutoff = CUTOFF,
                     target_frac = TARGET_FRAC, seed = SEED)
  saveRDS(res, f_tune)
}
ch <- res$chosen
message(sprintf("Chosen alpha=%.2f  selected=%.2f%%  pos=%d neg=%d",
                ch$alpha, ch$frac*100, length(ch$pos), length(ch$neg)))

## ---- per-cell labels (cheap; always derived from ch) ----
lab <- setNames(rep("Background", ncol(so)), colnames(so))
lab[ch$pos] <- "Scissor+"; lab[ch$neg] <- "Scissor-"
so$scissor <- factor(lab, levels = c("Scissor-","Background","Scissor+"))
saveRDS(list(chosen = ch, tuning = res$tuning, scissor_label = lab,
             coefs = ch$beta), file.path(OUT, "scissor_result.rds"))

## ---- significance: reliability test (checkpointed) ----
cell_num <- length(ch$pos) + length(ch$neg)
f_rt <- file.path(OUT, "reliability_test.rds")
if (file.exists(f_rt)) {
  message("[checkpoint] loading existing ", f_rt)
  rt <- readRDS(f_rt)
} else {
  message("Reliability test (n=", N_PERM_REL, ") ...")
  rt <- reliability_test_glmnet(inp, alpha = ch$alpha, cell_num = cell_num,
                                n = N_PERM_REL, nfold = 10, seed = 1)
  saveRDS(rt, f_rt)
}
message(sprintf("Reliability: real MSE=%.4f  null mean=%.4f  p=%.3f",
                rt$statistic, mean(rt$background), rt$p))

## ---- significance: selection permutation null (checkpointed) ----
f_null <- file.path(OUT, "permutation_null.rds")
if (file.exists(f_null)) {
  message("[checkpoint] loading existing ", f_null)
  pn <- readRDS(f_null)
} else {
  message("Selection permutation null (n=", N_PERM_NULL, ") ...")
  tiers <- inp$y; nf <- numeric(N_PERM_NULL); nd <- numeric(N_PERM_NULL)
  real_gap <- mean(so@meta.data[ch$pos,"tier"]) - mean(so@meta.data[ch$neg,"tier"])
  for (i in seq_len(N_PERM_NULL)) {
    set.seed(200 + i); yp <- tiers[sample(length(tiers))]
    f <- net_enet_path(inp$X, yp, inp$B, alpha = ch$alpha, target_frac = TARGET_FRAC)
    b <- f$beta; pos <- names(b)[b>0]; neg <- names(b)[b<0]
    nf[i] <- (length(pos)+length(neg))/ncol(inp$X)
    nd[i] <- (if(length(pos)) mean(so@meta.data[pos,"tier"]) else NA) -
             (if(length(neg)) mean(so@meta.data[neg,"tier"]) else NA)
  }
  pn <- list(real_gap = real_gap, null_frac = nf, null_dir = nd)
  saveRDS(pn, f_null)
}
message(sprintf("Selection null: real gap=%.3f  null gap mean=%.3f  p=%.3f",
                pn$real_gap, mean(pn$null_dir, na.rm=TRUE),
                mean(pn$null_dir >= pn$real_gap, na.rm=TRUE)))

## ---- characterization (checkpointed) ----
f_de <- file.path(OUT, "gradient_program_DE.csv")
if (file.exists(f_de)) {
  message("[checkpoint] ", f_de, " already exists, skipping DE + full-object save")
} else {
  message("Cell-type enrichment + gradient DE program ...")
  Idents(so) <- so$scissor
  de <- FindMarkers(so, ident.1 = "Scissor+", ident.2 = "Background",
                    logfc.threshold = 0.1, min.pct = 0.1)
  de$gene <- rownames(de)
  write.csv(de, f_de, row.names = FALSE)
  saveRDS(so, file.path(OUT, "reference_scissor_full.rds"))
}
message("FULL-CENSUS RUN COMPLETE. Outputs in ", OUT, "/")
