# Auto-extracted generating script
# Produces: meta_de_PPvsNN.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 668b37c3-d3fd-4cbc-bcb1-10ffdc22e8d3
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(edgeR)
  library(limma)
  library(SummarizedExperiment)
  library(data.table)
})

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/")

# Download gene-level RSE for each study
get_rse <- function(srp) {
  rse <- create_rse(subset(available_projects(organism="human"), project==srp & project_type=="data_sources"),
                    type="gene")
  assay(rse, "counts") <- transform_counts(rse)
  rse
}

studies <- c("SRP035988","SRP165679","SRP076982","SRP126422","SRP016583")
rses <- list()
for (s in studies) {
  message("downloading ", s, " ...")
  rses[[s]] <- get_rse(s)
  message("  ", s, ": ", nrow(rses[[s]]), " genes x ", ncol(rses[[s]]), " samples")
}

# Load metadata for each study
mds <- list()
for (s in studies) {
  mds[[s]] <- read_metadata(file_retrieve(locate_url(s, "data_sources/sra", "metadata", organism="human")))
}

# Classify samples per study
classify_samples <- function(srp, md) {
  a  <- md$sra.sample_attributes
  st <- md$sra.sample_title
  ext <- md$external_id
  cls <- rep("DROP", nrow(md))
  if (srp=="SRP035988") {
    cls[grepl("lesional psoriatic", a)] <- "PP"
    cls[grepl("normal skin", a)]        <- "NN"
  } else if (srp=="SRP165679") {
    dis  <- ifelse(grepl("^AD",st),"AD", ifelse(grepl("^PSO",st),"PSO","HC"))
    skin <- ifelse(grepl("non-lesional",st),"PN", ifelse(grepl("lesional",st),"PP","NN"))
    cls[dis=="PSO" & skin=="PP"] <- "PP"
    cls[dis=="PSO" & skin=="PN"] <- "PN"
    cls[dis=="HC"]               <- "NN"
  } else if (srp=="SRP076982") {
    cls[grepl("tissue type;;psoriasis", a)]   <- "PP"
    cls[grepl("tissue type;;uninvolved", a)]  <- "PN"
  } else if (srp=="SRP126422") {
    skinbx <- grepl("skin biopsy", a, ignore.case=TRUE)
    cls[skinbx & grepl("\\(PP\\)", a)] <- "PP"
    cls[skinbx & grepl("\\(NN\\)", a)] <- "NN"
    cls[skinbx & grepl("\\(PN\\)", a)] <- "PN"
  } else if (srp=="SRP016583") {
    cls[grepl("group;;lesional", a, ignore.case=TRUE)]     <- "PP"
    cls[grepl("group;;non-lesional", a, ignore.case=TRUE)] <- "PN"
  }
  data.frame(external_id=ext, srp=srp, class=cls, stringsAsFactors=FALSE)
}

cls_all <- rbindlist(lapply(studies, function(s) classify_samples(s, mds[[s]])))
setDT(cls_all)

# Run limma-voom DE for a given contrast pair within a study
run_de <- function(rse, srp, cls_dt, g1, g2) {
  cmap <- cls_dt[cls_dt$srp==srp,]
  cvec <- cmap$class[match(colnames(rse), cmap$external_id)]
  keep_s <- which(cvec %in% c(g1,g2))
  if (length(keep_s) < 6) return(NULL)
  grp <- factor(cvec[keep_s], levels=c(g2,g1))
  if (length(unique(grp)) < 2 || min(table(grp)) < 3) return(NULL)
  counts <- assay(rse, "counts")[, keep_s]
  gn <- rowData(rse)$gene_name
  dge <- DGEList(counts=counts, genes=data.frame(gene_name=gn))
  keep <- filterByExpr(dge, group=grp)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge, method="TMM")
  sym <- dge$genes$gene_name
  design <- model.matrix(~grp)
  v <- voom(dge, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
  tt$gene_name <- sym
  tt$SE <- sqrt(fit$s2.post) * fit$stdev.unscaled[,2]
  dt <- as.data.table(tt)
  dt <- dt[order(-AveExpr)][!duplicated(gene_name)]
  data.frame(gene=dt$gene_name, logFC=dt$logFC, SE=dt$SE, t=dt$t,
             P=dt$P.Value, adjP=dt$adj.P.Val, AveExpr=dt$AveExpr,
             n1=sum(grp==g1), n2=sum(grp==g2), stringsAsFactors=FALSE)
}

contrasts <- list(PPvsNN=c("PP","NN"), PNvsNN=c("PN","NN"), PNvsPP=c("PN","PP"))
de <- list()
for (cn in names(contrasts)) {
  g <- contrasts[[cn]]
  de[[cn]] <- list()
  for (s in names(rses)) {
    r <- run_de(rses[[s]], s, cls_all, g[1], g[2])
    if (!is.null(r)) { de[[cn]][[s]] <- r }
  }
  cat(sprintf("%s: studies contributing = %s\n", cn, paste(names(de[[cn]]), collapse=", ")))
}

# Build gene x study matrices of logFC (yi) and variance (vi=SE^2) for a contrast
build_mats <- function(delist) {
  studies <- names(delist)
  genes <- sort(unique(unlist(lapply(delist, function(x) x$gene))))
  Y <- matrix(NA, length(genes), length(studies), dimnames=list(genes, studies))
  V <- Y
  for (s in studies) {
    x <- delist[[s]]
    idx <- match(x$gene, genes)
    Y[idx, s] <- x$logFC
    V[idx, s] <- x$SE^2
  }
  list(Y=Y, V=V, studies=studies, genes=genes)
}

# Vectorized DerSimonian-Laird random-effects across rows (genes)
dl_meta <- function(Y, V) {
  W <- 1/V
  W[is.na(Y)] <- NA
  k <- rowSums(!is.na(Y))
  sumW  <- rowSums(W, na.rm=TRUE)
  sumWY <- rowSums(W*Y, na.rm=TRUE)
  mu_fe <- sumWY / sumW
  Q <- rowSums(W*(Y-mu_fe)^2, na.rm=TRUE)
  sumW2 <- rowSums(W^2, na.rm=TRUE)
  C <- sumW - sumW2/sumW
  tau2 <- pmax(0, (Q-(k-1))/C)
  tau2[k<2] <- NA
  Wr <- 1/(V + matrix(tau2, nrow(Y), ncol(Y)))
  Wr[is.na(Y)] <- NA
  sumWr  <- rowSums(Wr, na.rm=TRUE)
  mu_re  <- rowSums(Wr*Y, na.rm=TRUE)/sumWr
  se_re  <- sqrt(1/sumWr)
  z <- mu_re/se_re
  p <- 2*pnorm(-abs(z))
  I2 <- ifelse(Q>0 & k>1, pmax(0, 100*(Q-(k-1))/Q), 0)
  data.table(gene=rownames(Y), k=k, logFC=mu_re, SE=se_re, z=z, P=p,
             tau2=tau2, I2=I2, Q=Q)
}

meta_res <- list()
for (cn in names(de)) {
  m <- build_mats(de[[cn]])
  res <- dl_meta(m$Y, m$V)
  res <- res[k>=2]
  res[, FDR := p.adjust(P, "BH")]
  res <- res[order(P)]
  meta_res[[cn]] <- res
  cat(sprintf("%s: %d genes in >=2 studies; %d FDR<0.05\n", cn, nrow(res), sum(res$FDR<0.05)))
}

# Write meta_de_PPvsNN.csv
r <- copy(meta_res$PPvsNN)
r[, direction := ifelse(logFC>0,"up","down")]
setcolorder(r, c("gene","k","logFC","SE","z","P","FDR","I2","tau2","Q","direction"))
r[, `:=`(logFC=round(logFC,3), SE=round(SE,3), z=round(z,2),
         I2=round(I2,1), tau2=round(tau2,4), Q=round(Q,2))]
fwrite(r, "meta_de_PPvsNN.csv")