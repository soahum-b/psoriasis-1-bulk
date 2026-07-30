# Auto-extracted generating script
# Produces: meta_pathway_camera_PPvsNN.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): study_inventory.csv
# Source artifact version: ec629b64-7d48-4981-825b-69ede10987a4
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(edgeR)
  library(limma)
  library(msigdbr)
  library(data.table)
  library(SummarizedExperiment)
  library(recount3)
})

# Load study inventory
inv <- fread("study_inventory.csv")

# Sample classification
cls_all_raw <- list()

# Fetch metadata for each study
load_md <- function(srp) {
  read_metadata(file_retrieve(locate_url(srp, "data_sources/sra", "metadata", organism = "human")))
}

studies <- c("SRP035988", "SRP165679", "SRP076982", "SRP126422", "SRP016583")
mds <- lapply(studies, load_md)
names(mds) <- studies

classify_samples <- function(srp, md) {
  a   <- md$sra.sample_attributes
  st  <- md$sra.sample_title
  ext <- md$external_id
  cls <- rep("DROP", nrow(md))
  if (srp == "SRP035988") {
    cls[grepl("lesional psoriatic", a)] <- "PP"
    cls[grepl("normal skin", a)]        <- "NN"
  } else if (srp == "SRP165679") {
    dis  <- ifelse(grepl("^AD",  st), "AD",  ifelse(grepl("^PSO", st), "PSO", "HC"))
    skin <- ifelse(grepl("non-lesional", st), "PN", ifelse(grepl("lesional", st), "PP", "NN"))
    cls[dis == "PSO" & skin == "PP"] <- "PP"
    cls[dis == "PSO" & skin == "PN"] <- "PN"
    cls[dis == "HC"]                 <- "NN"
  } else if (srp == "SRP076982") {
    cls[grepl("tissue type;;psoriasis",  a)] <- "PP"
    cls[grepl("tissue type;;uninvolved", a)] <- "PN"
  } else if (srp == "SRP126422") {
    skinbx <- grepl("skin biopsy", a, ignore.case = TRUE)
    cls[skinbx & grepl("\\(PP\\)", a)] <- "PP"
    cls[skinbx & grepl("\\(NN\\)", a)] <- "NN"
    cls[skinbx & grepl("\\(PN\\)", a)] <- "PN"
  } else if (srp == "SRP016583") {
    cls[grepl("group;;lesional",     a, ignore.case = TRUE)] <- "PP"
    cls[grepl("group;;non-lesional", a, ignore.case = TRUE)] <- "PN"
  }
  data.frame(external_id = ext, srp = srp, class = cls, stringsAsFactors = FALSE)
}

cls_all <- rbindlist(lapply(studies, function(s) classify_samples(s, mds[[s]])))
setDT(cls_all)

# Download gene-level RSE for each study
get_rse <- function(srp) {
  rse <- create_rse(
    subset(available_projects(organism = "human"), project == srp & project_type == "data_sources"),
    type = "gene"
  )
  assay(rse, "counts") <- transform_counts(rse)
  rse
}

rses <- list()
for (s in studies) {
  rses[[s]] <- get_rse(s)
}

# Hallmark gene sets
H     <- msigdbr(species = "Homo sapiens", category = "H")
Hlist <- split(H$gene_symbol, H$gs_name)

# Per-study CAMERA for PPvsNN
camera_study <- function(rse, srp, g1 = "PP", g2 = "NN") {
  cmap <- cls_all[cls_all$srp == srp]
  cvec <- cmap$class[match(colnames(rse), cmap$external_id)]
  keep_s <- which(cvec %in% c(g1, g2))
  grp <- factor(cvec[keep_s], levels = c(g2, g1))
  if (length(unique(grp)) < 2 || min(table(grp)) < 3) return(NULL)
  counts <- assay(rse, "counts")[, keep_s]
  gn <- rowData(rse)$gene_name
  dge <- DGEList(counts = counts, genes = data.frame(gene_name = gn))
  keep <- filterByExpr(dge, group = grp)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, "TMM")
  design <- model.matrix(~grp)
  v      <- voom(dge, design)
  idx    <- ids2indices(Hlist, dge$genes$gene_name)
  cam    <- camera(v, idx, design, contrast = 2)
  dt     <- as.data.table(cam, keep.rownames = "pathway")
  dt[, study := srp]
  dt
}

studies_ppnn <- c("SRP035988", "SRP165679", "SRP126422")
cam_all <- rbindlist(lapply(studies_ppnn, function(s) camera_study(rses[[s]], s)), fill = TRUE)

# Pool by Stouffer weighted-Z
cam_all[, Pfloor := pmax(PValue, 1e-300)]
cam_all[, zraw   := qnorm(1 - pmin(Pfloor, 0.9999999) / 2)]
cam_all[, z      := ifelse(Direction == "Up", 1, -1) * pmin(zraw, 38)]
cam_all[!is.finite(z), z := ifelse(Direction == "Up", 1, -1) * 38]
cam_all[, w := sqrt(inv$total[match(study, inv$srp)])]

pooled <- cam_all[, .(
  k          = .N,
  Zc         = sum(w * z) / sqrt(sum(w^2)),
  nUp        = sum(Direction == "Up"),
  meanNGenes = round(mean(NGenes))
), by = pathway]

pooled[, p_comb := 2 * pnorm(-abs(pmin(Zc, 38)))]
pooled[, FDR    := p.adjust(p_comb, "BH")]
pooled <- pooled[order(-Zc)]

fwrite(pooled, "meta_pathway_camera_PPvsNN.csv")