# Auto-extracted generating script
# Produces: method_comparison.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_ranks.rds, gsea_results_all.rds, collectri_regulon.rds
# Source artifact version: 6cb23db2-4db8-43a5-82fb-64a055bf9e7e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR); library(limma); library(msigdbr); library(fgsea); library(data.table)})

dge <- readRDS("gsea_ranks.rds")

grp <- factor(dge$samples$group, levels=c("NN","PP"))
design <- model.matrix(~grp)
colnames(design) <- c("Intercept","PPvsNN")

v <- voom(dge, design, plot=FALSE)

col <- readRDS("collectri_regulon.rds")

stat3_targets <- unique(col$target[col$source=="STAT3"])
stat3_measured <- intersect(stat3_targets, rownames(v$E))

h_df <- msigdbr(species="Homo sapiens", category="H")
H <- split(h_df$gene_symbol, h_df$gs_name)
idx_H <- ids2indices(H, rownames(v$E))

idx_stat3 <- list(STAT3_regulon = which(rownames(v$E) %in% stat3_measured))

cam_H <- camera(v, idx_H, design, contrast="PPvsNN")
cam_stat3 <- camera(v, idx_stat3, design, contrast="PPvsNN")

set.seed(1)
roast_stat3 <- roast(v, idx_stat3$STAT3_regulon, design, contrast="PPvsNN", nrot=9999)

ranks <- readRDS("gsea_ranks.rds")

gsea_all <- readRDS("gsea_results_all.rds")

set.seed(1)
fg_stat3 <- fgsea(pathways=list(STAT3_regulon=stat3_measured), stats=ranks,
                  minSize=5, maxSize=1000)

gseaH <- as.data.frame(gsea_all$H)[,c("pathway","NES","padj")]
camH  <- data.frame(pathway=rownames(cam_H), cam_dir=cam_H$Direction,
                    cam_FDR=cam_H$FDR, cam_NGenes=cam_H$NGenes)
mrg <- merge(gseaH, camH, by="pathway")
mrg$gsea_sig <- mrg$padj < 0.05
mrg$cam_sig  <- mrg$cam_FDR < 0.05
mrg$agree <- (mrg$NES>0 & mrg$cam_dir=="Up" & mrg$gsea_sig & mrg$cam_sig) |
             (mrg$NES<0 & mrg$cam_dir=="Down" & mrg$gsea_sig & mrg$cam_sig)

saveRDS(list(mrg=mrg, fg_stat3=fg_stat3), "method_comparison.rds")