# Auto-extracted generating script
# Produces: camera_perstudy_k4_full.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): study_inventory.csv
# Source artifact version: 2b8808b1-4117-4389-b89c-9e3e7714bacf
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(edgeR)
  library(limma)
  library(msigdbr)
  library(data.table)
})

studies <- c("SRP035988", "SRP165679", "SRP126422", "SRP065812")

classify <- function(srp, md) {
  a <- md$sra.sample_attributes; st <- md$sra.sample_title; ext <- md$external_id
  cls <- rep("DROP", nrow(md))
  if (srp == "SRP035988") {
    cls[grepl("lesional psoriatic", a)] <- "PP"
    cls[grepl("normal skin", a)] <- "NN"
  } else if (srp == "SRP165679") {
    dis <- ifelse(grepl("^PSO", st), "PSO", ifelse(grepl("^AD", st), "AD", "HC"))
    skin <- ifelse(grepl("non-lesional", st), "PN", ifelse(grepl("lesional", st), "PP", "NN"))
    cls[dis == "PSO" & skin == "PP"] <- "PP"
    cls[dis == "HC"] <- "NN"
  } else if (srp == "SRP126422") {
    sb <- grepl("skin biopsy", a, ignore.case = TRUE)
    cls[sb & grepl("\\(PP\\)", a)] <- "PP"
    cls[sb & grepl("\\(NN\\)", a)] <- "NN"
  } else if (srp == "SRP065812") {
    cls[grepl("before adalimumab", a, ignore.case = TRUE) & grepl("lesional", a, ignore.case = TRUE)] <- "PP"
    cls[grepl("healthy control", a, ignore.case = TRUE) | grepl("normal skin", a, ignore.case = TRUE)] <- "NN"
    cls[grepl("after|post", a, ignore.case = TRUE)] <- "DROP"
  }
  data.table(external_id = ext, srp = srp, class = cls)
}

get_rse <- function(srp) {
  rse <- create_rse(subset(available_projects(organism = "human"),
                           project == srp & project_type == "data_sources"), type = "gene")
  assay(rse, "counts") <- transform_counts(rse)
  rse
}

load_md <- function(srp) read_metadata(file_retrieve(locate_url(srp, "data_sources/sra", "metadata", organism = "human")))

H <- msigdbr(species = "Homo sapiens", collection = "H")
Hlist <- split(H$gene_symbol, H$gs_name)

cam_one <- function(srp) {
  md <- load_md(srp)
  cl <- classify(srp, md)
  rse <- get_rse(srp)
  cvec <- cl$class[match(colnames(rse), cl$external_id)]
  keep_s <- which(cvec %in% c("NN", "PP"))
  grp <- factor(cvec[keep_s], levels = c("NN", "PP"))
  if (length(unique(grp)) < 2 || min(table(grp)) < 3) return(NULL)
  dge <- DGEList(counts = assay(rse, "counts")[, keep_s],
                 genes = data.frame(gene_name = rowData(rse)$gene_name))
  keep <- filterByExpr(dge, group = grp)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, "TMM")
  design <- model.matrix(~grp)
  v <- voom(dge, design)
  idx <- ids2indices(Hlist, dge$genes$gene_name)
  cam <- camera(v, idx, design, contrast = 2)
  dt <- as.data.table(cam, keep.rownames = "pathway")
  dt[, `:=`(study = srp, n = length(grp), nNN = sum(grp == "NN"), nPP = sum(grp == "PP"))]
  dt
}

set.seed(1)
cam_all <- rbindlist(lapply(studies[1:3], cam_one), fill = TRUE)

# SRP065812 with corrected classifier (already handled in classify())
md <- load_md("SRP065812")
cl <- classify("SRP065812", md)
rse <- get_rse("SRP065812")
cvec <- cl$class[match(colnames(rse), cl$external_id)]
keep_s <- which(cvec %in% c("NN", "PP"))
grp <- factor(cvec[keep_s], levels = c("NN", "PP"))
dge <- DGEList(counts = assay(rse, "counts")[, keep_s],
               genes = data.frame(gene_name = rowData(rse)$gene_name))
keep <- filterByExpr(dge, group = grp)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, "TMM")
design <- model.matrix(~grp)
v <- voom(dge, design)
idx <- ids2indices(Hlist, dge$genes$gene_name)
cam <- camera(v, idx, design, contrast = 2)
dt <- as.data.table(cam, keep.rownames = "pathway")
dt[, `:=`(study = "SRP065812", n = length(grp), nNN = sum(grp == "NN"), nPP = sum(grp == "PP"))]

cam_all <- cam_all[study != "SRP065812"]
cam_all <- rbindlist(list(cam_all, dt), fill = TRUE)
saveRDS(cam_all, "camera_perstudy_k4_full.rds")