#!/usr/bin/env Rscript
# paired_blood_skin.R — within-patient blood<->skin signature analysis
#
# ERP110814 (whole blood) and ERP110816 (skin) are the SAME 10 patients
# (individual IDs AFOX02..AFOY11) at weeks 0/1/12 on etanercept. This makes a
# within-patient design possible with no cross-study batch confound -- the
# reason the earlier ERP110814-vs-GTEx case/control contrast was abandoned.
#
# Usage:  Rscript paired_blood_skin.R [outdir]
# Env:    scissor-r or r441_env (needs recount3, SummarizedExperiment, edgeR)

suppressPackageStartupMessages({
  library(recount3); library(SummarizedExperiment); library(edgeR)
})
outdir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

get_rse <- function(p) create_rse_manual(
  project = p, project_home = "data_sources/sra",
  organism = "human", annotation = "gencode_v26", type = "gene")

fld <- function(x, k) sub(paste0(".*", k, ";;([^|]*).*"), "\\1",
                          ifelse(grepl(k, x), x, NA))

meta <- function(rse) {
  m <- as.data.frame(colData(rse))
  data.frame(run   = m$external_id,
             indiv = fld(m$sra.sample_attributes, "individual"),
             time  = sub(".*Experimental Factor: time;;([^|]*).*", "\\1",
                         m$sra.experiment_attributes),
             site  = fld(m$sra.sample_attributes, "sampling site"),
             stringsAsFactors = FALSE)
}

# log-CPM after filterByExpr + TMM, rownames = gene symbols
norm_cpm <- function(rse) {
  cts <- assay(rse, "raw_counts"); rownames(cts) <- rowData(rse)$gene_name
  d <- DGEList(cts); d <- d[filterByExpr(d), , keep.lib.sizes = FALSE]
  cpm(calcNormFactors(d, method = "TMM"), log = TRUE, prior.count = 1)
}

# a priori signatures (see notes/references.md for provenance)
sigs <- list(
  blood_platelet = c("ITGA2B","ITGB3","PF4","PPBP","GP1BA","GP9","TBXAS1",
                     "PTGS1","SELP","TREML1","CLEC1B","MYL9"),   # Garshick 2020
  blood_ifn      = c("IFI6","IFI44","IFI44L","ISG15","MX1","MX2","OAS1","OAS2",
                     "OAS3","IFIT1","IFIT3","RSAD2","STAT1"),
  skin_endo      = c("PECAM1","CDH5","VWF","CLDN5","EGFL7","ACKR1","PLVAP",
                     "AQP1","RAMP2","EMCN","KDR","FLT1"),
  skin_il1_infl  = c("IL1R1","CASP1","CASP4","GSDMD","PYCARD","AIM2","IL6",
                     "NFKB1","ICAM1","SELE"),
  skin_il36      = c("IL36A","IL36B","IL36G","IL36RN"))

score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  colMeans(t(scale(t(mat[g, , drop = FALSE]))), na.rm = TRUE)
}

rse_b <- get_rse("ERP110814"); rse_s <- get_rse("ERP110816")
mb <- meta(rse_b); ms <- meta(rse_s)
stopifnot(setequal(na.omit(mb$indiv), na.omit(ms$indiv)))   # same patients
lb <- norm_cpm(rse_b); ls_ <- norm_cpm(rse_s)

sc_b <- data.frame(mb, platelet = score(lb, sigs$blood_platelet),
                       ifn      = score(lb, sigs$blood_ifn))
sc_s <- data.frame(ms, endo = score(ls_, sigs$skin_endo),
                       il1  = score(ls_, sigs$skin_il1_infl),
                       il36 = score(ls_, sigs$skin_il36))

# lesional skin, averaged per patient x timepoint, joined to blood
agg <- aggregate(cbind(endo, il1, il36) ~ indiv + time + site, sc_s, mean)
les <- agg[agg$site == "lesion", c("indiv","time","endo","il1","il36")]
names(les)[3:5] <- paste0("skin_", names(les)[3:5])
m <- merge(sc_b[, c("indiv","time","platelet","ifn")], les,
           by = c("indiv","time"))
m$tnum <- as.numeric(as.character(m$time))
vars <- c("platelet","ifn","skin_endo","skin_il1","skin_il36")

# patient-centring removes between-patient variation -> within-patient inference
cen <- m; for (v in vars) cen[[v]] <- cen[[v]] - ave(cen[[v]], cen$indiv)

tr <- do.call(rbind, lapply(vars, function(v) {
  f <- summary(lm(cen[[v]] ~ cen$tnum))$coefficients
  data.frame(signature = v, slope_per_week = f[2,1], p = f[2,4])
}))
tr$p_BH <- p.adjust(tr$p, "BH")

pr <- expand.grid(blood = c("platelet","ifn"),
                  skin  = c("skin_endo","skin_il1","skin_il36"),
                  stringsAsFactors = FALSE)
wc <- do.call(rbind, lapply(seq_len(nrow(pr)), function(i) {
  ct <- suppressWarnings(cor.test(cen[[pr$blood[i]]], cen[[pr$skin[i]]],
                                 method = "spearman"))
  data.frame(pr[i,], rho_within = ct$estimate, p = ct$p.value)
}))
wc$p_BH <- p.adjust(wc$p, "BH")

write.csv(m,  file.path(outdir, "paired_blood_skin_scores.csv"), row.names = FALSE)
write.csv(tr, file.path(outdir, "paired_trajectories.csv"),      row.names = FALSE)
write.csv(wc, file.path(outdir, "paired_within_patient_corr.csv"), row.names = FALSE)
cat("obs:", nrow(m), " patients:", length(unique(m$indiv)), "\n"); print(tr); print(wc)
