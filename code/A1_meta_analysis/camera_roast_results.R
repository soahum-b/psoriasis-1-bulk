# Auto-extracted generating script
# Produces: camera_roast_results.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_results_all.rds, collectri_regulon.rds, de_results_full.rds, dge_filt_norm.rds
# Source artifact version: 0144430b-37ba-456c-9103-4469bb684952
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR); library(limma); library(msigdbr); library(fgsea); library(data.table)})

NNblue <- "#4C72B0"; PPred <- "#C44E52"

dge <- readRDS("dge_filt_norm.rds")
de_tab <- readRDS("de_results_full.rds")
col <- readRDS("collectri_regulon.rds")

grp <- factor(dge$samples$group, levels=c("NN","PP"))
design <- model.matrix(~grp)
colnames(design) <- c("Intercept","PPvsNN")

v <- voom(dge, design, plot=FALSE)

h_df <- msigdbr(species="Homo sapiens", category="H")
H <- split(h_df$gene_symbol, h_df$gs_name)
r_df <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
RE <- split(r_df$gene_symbol, r_df$gs_name)

stat3_targets <- unique(col$target[col$source=="STAT3"])
stat3_measured <- intersect(stat3_targets, rownames(v$E))

idx_H  <- ids2indices(H,  rownames(v$E))
idx_RE <- ids2indices(RE, rownames(v$E))
idx_stat3 <- list(STAT3_regulon = which(rownames(v$E) %in% stat3_measured))

cam_H <- camera(v, idx_H, design, contrast="PPvsNN")
cam_stat3 <- camera(v, idx_stat3, design, contrast="PPvsNN")

set.seed(1)
roast_stat3 <- roast(v, idx_stat3$STAT3_regulon, design, contrast="PPvsNN", nrot=9999)

saveRDS(list(cam_H=cam_H, cam_stat3=cam_stat3, roast_stat3=roast_stat3,
             stat3_measured=stat3_measured, idx_H=idx_H, idx_RE=idx_RE),
        "camera_roast_results.rds")