# Auto-extracted generating script
# Produces: camera_perstudy_PPvsNN.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds, dge_filt_norm.rds
# Source artifact version: b1161da7-572a-42da-a4e7-e1e5d6d9b7e6
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR); library(limma); library(msigdbr); library(data.table); library(SummarizedExperiment); library(recount3)})

cls_all_raw <- readRDS("per_study_de.rds")
de <- cls_all_raw
# cls_all_raw is per_study_de.rds - but we need cls_all for sample classification
# We need to reconstruct cls_all from the metadata

# Load the RSEs for the 3 studies used in camera_study PPvsNN
studies_ppnn <- c("SRP035988","SRP165679","SRP126422")

get_rse <- function(srp) {
  rse <- create_rse(subset(available_projects(organism="human"), project==srp & project_type=="data_sources"),
                    type="gene")
  assay(rse, "counts") <- transform_counts(rse)
  rse
}

rses <- list()
for (s in studies_ppnn) {
  message("downloading ", s, " ...")
  rses[[s]] <- get_rse(s)
}

mds <- lapply(studies_ppnn, function(srp)
  read_metadata(file_retrieve(locate_url(srp, "data_sources/sra", "metadata", organism="human")))
)
names(mds) <- studies_ppnn

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
  } else if (srp=="SRP126422") {
    skinbx <- grepl("skin biopsy", a, ignore.case=TRUE)
    cls[skinbx & grepl("\\(PP\\)", a)] <- "PP"
    cls[skinbx & grepl("\\(NN\\)", a)] <- "NN"
    cls[skinbx & grepl("\\(PN\\)", a)] <- "PN"
  }
  data.frame(external_id=ext, srp=srp, class=cls, stringsAsFactors=FALSE)
}

cls_all <- rbindlist(lapply(studies_ppnn, function(s) classify_samples(s, mds[[s]])))
setDT(cls_all)

# Hallmark sets (symbols)
H <- msigdbr(species="Homo sapiens", category="H")
Hlist <- split(H$gene_symbol, H$gs_name)

# per-study CAMERA for PPvsNN
camera_study <- function(rse, srp, g1="PP", g2="NN") {
  cmap <- cls_all[cls_all$srp==srp]
  cvec <- cmap$class[match(colnames(rse), cmap$external_id)]
  keep_s <- which(cvec %in% c(g1,g2))
  grp <- factor(cvec[keep_s], levels=c(g2,g1))
  if (length(unique(grp))<2 || min(table(grp))<3) return(NULL)
  counts <- assay(rse,"counts")[,keep_s]
  gn <- rowData(rse)$gene_name
  dge <- DGEList(counts=counts, genes=data.frame(gene_name=gn))
  keep <- filterByExpr(dge, group=grp); dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge,"TMM")
  design <- model.matrix(~grp); v <- voom(dge, design)
  idx <- ids2indices(Hlist, dge$genes$gene_name)
  cam <- camera(v, idx, design, contrast=2)
  dt <- as.data.table(cam, keep.rownames="pathway")
  dt[, study:=srp]; dt
}

cam_all <- rbindlist(lapply(studies_ppnn, function(s) camera_study(rses[[s]], s)), fill=TRUE)
saveRDS(cam_all, "camera_perstudy_PPvsNN.rds")