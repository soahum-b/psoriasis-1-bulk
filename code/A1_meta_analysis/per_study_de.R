# Auto-extracted generating script
# Produces: per_study_de.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): dge_filt_norm.rds
# Source artifact version: 8b4e71d3-763b-4c68-bcb0-0e21d1c0b496
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

studies <- c("SRP035988","SRP165679","SRP076982","SRP126422","SRP016583")

get_rse <- function(srp) {
  rse <- create_rse(subset(available_projects(organism="human"), project==srp & project_type=="data_sources"),
                    type="gene")
  assay(rse, "counts") <- transform_counts(rse)
  rse
}
rses <- list()
for (s in studies) {
  message("downloading ", s, " ...")
  rses[[s]] <- get_rse(s)
  message("  ", s, ": ", nrow(rses[[s]]), " genes x ", ncol(rses[[s]]), " samples")
}

mds <- lapply(studies, function(srp) {
  read_metadata(file_retrieve(locate_url(srp, "data_sources/sra", "metadata", organism="human")))
})
names(mds) <- studies

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

run_de <- function(rse, srp, cls_dt, g1, g2) {
  cmap <- cls_dt[cls_dt$srp==srp, ]
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
saveRDS(de, "per_study_de.rds")
cat("\nsaved per_study_de.rds\n")