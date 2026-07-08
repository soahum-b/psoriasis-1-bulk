#!/usr/bin/env Rscript
# run_pipeline.R
# End-to-end LOCAL backbone driver (laptop-scale, 20k-cell stratified subset).
# Runs every step from raw data to results/ and figures/. For the full
# 89k-cell census use code/run_full_census_cluster.R instead.
#   Rscript code/run_pipeline.R              # full rebuild
#   Rscript code/run_pipeline.R --from-checkpoint   # skip Steps 1-2 if the
#                                                    # processed reference exists
# Assumes code/00_download_data.R has already populated data/.

suppressMessages({library(Seurat); library(Matrix); library(glmnet);
                  library(preprocessCore); library(org.Hs.eg.db);
                  library(AnnotationDbi); library(edgeR);
                  library(SummarizedExperiment); library(data.table);
                  library(ggplot2); library(patchwork)})
source("code/scissor_glmnet_solver.R")
source("code/scissor_run.R")
source("code/scissor_reliability.R")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
set.seed(123)
from_ckpt <- "--from-checkpoint" %in% commandArgs(TRUE)

## ------------------------------------------------------------------ Steps 1-2
if (from_ckpt && file.exists("results/reference_processed.rds")) {
  message("== Steps 1-2 skipped (--from-checkpoint): loading processed reference ==")
  so <- readRDS("results/reference_processed.rds")
} else {
  message("== Step 1: assemble scRNA reference from data/raw/*.csv.gz ==")
  man <- read.csv("data/sample_manifest.csv", stringsAsFactors = FALSE)
  mats <- lapply(seq_len(nrow(man)), function(i) {
    f <- file.path("data/raw", man$file[i])
    m <- as.matrix(fread(f), rownames = 1)
    m <- as(m, "CsparseMatrix")
    colnames(m) <- paste0(man$gsm[i], "_", colnames(m))   # unique cell ids
    m
  })
  genes <- rownames(mats[[1]])
  stopifnot(all(vapply(mats, function(m) identical(rownames(m), genes), logical(1))))
  counts <- do.call(cbind, mats)
  meta <- data.frame(row.names = colnames(counts))
  cellsplit <- sub("_.*$", "", colnames(counts))
  mi <- match(cellsplit, man$gsm)
  meta$donor <- man$sample[mi]; meta$condition <- man$condition[mi]
  meta$tier <- c(NN = 0, PN = 1, PP = 2)[meta$condition]
  so <- CreateSeuratObject(counts, meta.data = meta, project = "Ma2023_psoriasis")
  saveRDS(so, "results/reference_raw.rds")
  message("  merged: ", nrow(so), " genes x ", ncol(so), " cells")

  message("== Step 2: Ensembl->symbol, QC, normalize, cluster, annotate ==")
  sym <- mapIds(org.Hs.eg.db, keys = sub("\\.\\d+$", "", rownames(so)),
                column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  keep <- !is.na(sym)
  ck <- GetAssayData(so, layer = "counts")[keep, ]; syms <- factor(sym[keep])
  G <- sparseMatrix(i = as.integer(syms), j = seq_along(syms), x = 1,
                    dims = c(nlevels(syms), length(syms)),
                    dimnames = list(levels(syms), NULL))
  cs <- as(G %*% ck, "CsparseMatrix")
  so <- CreateSeuratObject(cs, meta.data = so@meta.data)
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- subset(so, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                            nCount_RNA >= 500 & percent.mt < 20)
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
  so <- ScaleData(so, verbose = FALSE)
  so <- RunPCA(so, npcs = 30, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)
  so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
  so <- RunUMAP(so, dims = 1:20, verbose = FALSE)
  markers <- list(
    Keratinocyte=c("KRT14","KRT5","KRT10","KRT1","KRT15","KRTDAP","SBSN"),
    Fibroblast=c("COL1A1","COL1A2","COL3A1","PDGFRA","LUM","DCN"),
    Myeloid=c("LYZ","CD68","CD14","AIF1","ITGAX","FCGR3A"),
    DC=c("CD1C","CLEC9A","LAMP3","CD207"),
    Tcell=c("CD3D","CD3E","CD2","TRAC","IL7R"), NK=c("NKG7","GNLY","KLRD1"),
    Bcell=c("MS4A1","CD79A","IGHG1","MZB1"),
    Endothelial=c("PECAM1","VWF","CLDN5","CCL21"),
    Melanocyte=c("MLANA","PMEL","TYRP1","DCT"), Mast=c("TPSAB1","TPSB2","CPA3"),
    Muscle=c("ACTA2","TAGLN","MYH11","DES"))
  markers <- lapply(markers, intersect, rownames(so))
  avg <- AverageExpression(so, features = unique(unlist(markers)),
                           group.by = "seurat_clusters", layer = "data")$RNA
  zz <- t(scale(t(as.matrix(avg))))
  score <- sapply(markers, function(g) colMeans(zz[intersect(g,rownames(zz)),,drop=FALSE], na.rm=TRUE))
  assign <- colnames(score)[apply(score, 1, which.max)]
  names(assign) <- sub("^g", "", rownames(score))
  so$celltype <- factor(unname(assign[as.character(so$seurat_clusters)]))
  saveRDS(so, "results/reference_processed.rds")
  message("  processed: ", ncol(so), " cells, ", nlevels(so$celltype), " lineages")
}

## -------------------------------------------------------------- Step 3 (anchor)
message("== Step 3: build bulk phenotype anchor (SRP165679) ==")
if (!file.exists("results/bulk_anchor.rds")) {
  rse <- readRDS("data/rse_SRP165679.rds")
  # sample_classification.rds ships in data/ or is fetched with the recount pull
  sc <- readRDS("data/sample_classification.rds")
  sc <- sc[sc$srp == "SRP165679" & sc$class %in% c("NN","PN","PP"), ]
  ext <- colData(rse)$external_id; kp <- ext %in% sc$external_id
  rse <- rse[, kp]; cls <- sc$class[match(ext[kp], sc$external_id)]
  d <- DGEList(assay(rse,"counts"), group = cls); d <- calcNormFactors(d)
  d <- d[filterByExpr(d, group = cls), , keep.lib.sizes = FALSE]
  lc <- cpm(d, log = TRUE, prior.count = 1); rownames(lc) <- rowData(rse)$gene_name
  lc <- lc[order(rowMeans(lc), decreasing = TRUE), ]
  lc <- lc[!duplicated(rownames(lc)) & !is.na(rownames(lc)) & rownames(lc) != "", ]
  saveRDS(list(logcpm = lc, tier = c(NN=0,PN=1,PP=2)[cls], cls = cls,
               samples = ext[kp], srp = "SRP165679"), "results/bulk_anchor.rds")
}
anch <- readRDS("results/bulk_anchor.rds")

## ------------------------------------------------------- Step 4 (20k subset)
message("== Step 4: stratified 20k subsample + SNN ==")
set.seed(42)
md <- so@meta.data; md$strata <- paste(md$celltype, md$condition, sep = "|")
fr <- 20000 / ncol(so)
pick <- unlist(lapply(split(rownames(md), md$strata), function(ids)
          sample(ids, min(length(ids), max(20, round(length(ids)*fr))))))
sub <- subset(so, cells = pick)
sub <- FindNeighbors(sub, dims = 1:20, verbose = FALSE)
saveRDS(sub, "results/reference_subset20k.rds")

## ----------------------------------------------- Step 5 (Scissor selection)
message("== Step 5: Scissor inputs + alpha tuning ==")
inp <- prepare_scissor_inputs(anch$logcpm, sub, y = anch$tier, hvg_only = TRUE,
                              save_file = "results/scissor_inputs.rds")
res <- run_scissor(inp, alpha_grid = c(0.05,0.1,0.2,0.3,0.4,0.5),
                   cutoff = 0.20, target_frac = 0.15, seed = 123)
ch <- res$chosen
lab <- setNames(rep("Background", ncol(sub)), colnames(sub))
lab[ch$pos] <- "Scissor+"; lab[ch$neg] <- "Scissor-"
saveRDS(list(chosen = ch, tuning = res$tuning, scissor_label = lab,
             coefs = ch$beta, cells = names(lab)), "results/scissor_result.rds")
write.csv(res$tuning, "results/scissor_alpha_tuning.csv", row.names = FALSE)
message(sprintf("  chosen alpha=%.2f  selected=%.2f%%", ch$alpha, ch$frac*100))

## --------------------------------------------- Step 6 (reliability + null)
message("== Step 6: reliability test + permutation null ==")
cell_num <- length(ch$pos) + length(ch$neg)
rt <- reliability_test_glmnet(inp, ch$alpha, cell_num, n = 100, nfold = 10, seed = 1)
saveRDS(rt, "results/reliability_test.rds")
tiers <- inp$y; np <- 30; nd <- numeric(np)
real_gap <- mean(sub@meta.data[ch$pos,"tier"]) - mean(sub@meta.data[ch$neg,"tier"])
for (i in seq_len(np)) { set.seed(200+i); yp <- tiers[sample(length(tiers))]
  f <- net_enet_path(inp$X, yp, inp$B, alpha = ch$alpha, target_frac = 0.15)
  b <- f$beta; pos <- names(b)[b>0]; neg <- names(b)[b<0]
  nd[i] <- (if(length(pos)) mean(sub@meta.data[pos,"tier"]) else NA) -
           (if(length(neg)) mean(sub@meta.data[neg,"tier"]) else NA) }
saveRDS(list(real_gap = real_gap, null_dir = nd), "results/permutation_null.rds")
message(sprintf("  reliability p=%.3f | null gap p=%.3f",
                rt$p, mean(nd >= real_gap, na.rm = TRUE)))

## ---------------------------------------------- Step 7 (characterization)
message("== Step 7: cell-type enrichment + gradient DE program ==")
sub$scissor <- factor(lab, levels = c("Scissor-","Background","Scissor+"))
Idents(sub) <- sub$scissor
de <- FindMarkers(sub, ident.1 = "Scissor+", ident.2 = "Background",
                  logfc.threshold = 0.1, min.pct = 0.1)
de$gene <- rownames(de); write.csv(de, "results/gradient_program_DE.csv", row.names = FALSE)
saveRDS(sub, "results/reference_subset20k_umap.rds")
message("PIPELINE COMPLETE. Results in results/, checkpoints saved per step.")
