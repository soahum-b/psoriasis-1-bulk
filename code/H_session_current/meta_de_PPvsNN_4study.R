# Auto-extracted generating script
# Produces: meta_de_PPvsNN_4study.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds, SRP065812_rse_clean.rds
# Source artifact version: 1ec80908-9a53-4351-8900-a793b995dabd
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import subprocess, shutil, os

# Write R script
r_code = r"""
suppressMessages({library(edgeR); library(limma); library(SummarizedExperiment)})
rse812 <- readRDS("SRP065812_rse_clean.rds")

de_one <- function(rse, cls, cA, cB) {
  keep_s <- which(cls %in% c(cA, cB))
  grp <- factor(cls[keep_s], levels = c(cB, cA))
  counts <- assay(rse, "counts")[, keep_s]
  gn <- rowData(rse)$gene_name
  dge <- DGEList(counts = counts, genes = data.frame(gene_name = gn))
  keep <- filterByExpr(dge, group = grp); dge <- dge[keep,, keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, "TMM")
  des <- model.matrix(~ grp)
  v <- voom(dge, des); fit <- eBayes(lmFit(v, des))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = dge$genes$gene_name, logFC = tt$logFC,
             SE = tt$logFC / tt$t, t = tt$t, P = tt$P.Value,
             adjP = p.adjust(tt$P.Value, "BH"), AveExpr = tt$AveExpr,
             n1 = sum(grp==cA), n2 = sum(grp==cB))
}
cls812 <- colData(rse812)$grp
de812 <- de_one(rse812, cls812, "PP", "NN")

psd <- readRDS("per_study_de.rds")
de_list4 <- c(psd$PPvsNN, list(SRP065812 = de812))

meta_dl <- function(de_list) {
  genes <- sort(unique(unlist(lapply(de_list, function(d) d$gene))))
  Y <- sapply(de_list, function(d) d$logFC[match(genes, d$gene)])
  V <- sapply(de_list, function(d) (d$SE[match(genes, d$gene)])^2)
  rownames(Y) <- rownames(V) <- genes
  out <- t(apply(cbind(Y, V), 1, function(row) {
    k <- ncol(Y); yi <- row[1:k]; vi <- row[(k+1):(2*k)]
    ok <- is.finite(yi) & is.finite(vi) & vi > 0
    if (sum(ok) < 2) return(c(logFC=NA,SE=NA,P=NA,I2=NA,k=sum(ok)))
    yi <- yi[ok]; vi <- vi[ok]; wi <- 1/vi
    ybar <- sum(wi*yi)/sum(wi)
    Q <- sum(wi*(yi-ybar)^2); df <- length(yi)-1
    C <- sum(wi) - sum(wi^2)/sum(wi)
    tau2 <- max(0, (Q-df)/C)
    ws <- 1/(vi+tau2); mu <- sum(ws*yi)/sum(ws); se <- sqrt(1/sum(ws))
    I2 <- max(0,(Q-df)/Q)*100
    c(logFC=mu, SE=se, P=2*pnorm(-abs(mu/se)), I2=I2, k=length(yi))
  }))
  d <- data.frame(gene=genes, out); d$FDR <- p.adjust(d$P, "BH"); d
}
meta4 <- meta_dl(de_list4)
write.csv(meta4, "meta_de_PPvsNN_4study.csv", row.names=FALSE)
cat("saved 4-study results\n")
"""

shutil.copy(host.artifact_path("8b4e71d3-763b-4c68-bcb0-0e21d1c0b496"), "per_study_de.rds")
shutil.copy(host.artifact_path("d767443f-052b-4328-8d7f-c9c264d11b0f"), "SRP065812_rse_clean.rds")

with open("run_meta4.R", "w") as f:
    f.write(r_code)

result = subprocess.run(["Rscript", "run_meta4.R"], capture_output=True, text=True)
print(result.stdout)
print(result.stderr)