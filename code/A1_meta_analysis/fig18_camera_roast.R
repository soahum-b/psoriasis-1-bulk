# Auto-extracted generating script
# Produces: fig18_camera_roast.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_ranks.rds, gsea_results_all.rds, collectri_regulon.rds, de_results_full.rds, dge_filt_norm.rds
# Source artifact version: 49b8307e-418d-437d-b0f4-bc0bc0f75d3b
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR); library(limma); library(msigdbr); library(fgsea); library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})

NNblue <- "#4C72B0"; PPred <- "#C44E52"

dge <- readRDS("dge_filt_norm.rds")
de_tab <- readRDS("de_results_full.rds")
col <- readRDS("collectri_regulon.rds")
ranks <- readRDS("gsea_ranks.rds")
gsea_all <- readRDS("gsea_results_all.rds")

grp <- factor(dge$samples$group, levels=c("NN","PP"))
design <- model.matrix(~grp)
colnames(design) <- c("Intercept","PPvsNN")

v <- voom(dge, design, plot=FALSE)

h_df <- msigdbr(species="Homo sapiens", category="H")
H <- split(h_df$gene_symbol, h_df$gs_name)
idx_H <- ids2indices(H, rownames(v$E))

stat3_targets <- unique(col$target[col$source=="STAT3"])
stat3_measured <- intersect(stat3_targets, rownames(v$E))
idx_stat3 <- list(STAT3_regulon = which(rownames(v$E) %in% stat3_measured))

cam_H <- camera(v, idx_H, design, contrast="PPvsNN")
cam_stat3 <- camera(v, idx_stat3, design, contrast="PPvsNN")

set.seed(1)
roast_stat3 <- roast(v, idx_stat3$STAT3_regulon, design, contrast="PPvsNN", nrot=9999)

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

grey <- "#B0B0B0"

mrg$neglog_cam <- -log10(mrg$cam_FDR)
mrg$cat <- ifelse(mrg$agree, "Both sig, same direction",
           ifelse(mrg$gsea_sig | mrg$cam_sig, "One method sig", "Neither sig"))
lab <- mrg[mrg$agree & (mrg$neglog_cam > 8 | abs(mrg$NES) > 2.3),]

pA <- ggplot(mrg, aes(NES, neglog_cam)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color=grey) +
  geom_vline(xintercept=0, color=grey) +
  geom_point(aes(color=cat, size=cam_NGenes), alpha=0.8) +
  geom_text_repel(data=lab, aes(label=gsub("HALLMARK_","",pathway)),
                  size=2.5, max.overlaps=20, min.segment.length=0) +
  scale_color_manual(values=c("Both sig, same direction"=PPred,
                              "One method sig"="#DD8452","Neither sig"=grey)) +
  scale_size_continuous(range=c(1.5,6), name="Set size") +
  labs(x="GSEA NES (fgsea, competitive)", y="CAMERA  -log10(FDR)",
       title="A. Hallmark: GSEA vs CAMERA concordance",
       subtitle="27/50 sets significant in both, same direction; CAMERA corrects inter-gene correlation",
       color=NULL) +
  theme_bw(base_size=10) + theme(legend.position="right",
       plot.subtitle=element_text(size=7.5))

db <- data.frame(
  method=factor(c("GSEA (competitive)","CAMERA (competitive,\ncorrelation-corrected)","ROAST (self-contained)"),
                levels=c("GSEA (competitive)","CAMERA (competitive,\ncorrelation-corrected)","ROAST (self-contained)")),
  p=c(fg_stat3$pval[1], cam_stat3$PValue[1], roast_stat3$p.value["Up","P.Value"]),
  type=c("competitive","competitive","self-contained"))
db$negp <- -log10(db$p)

pB <- ggplot(db, aes(method, negp, fill=type)) +
  geom_col(width=0.6) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color=grey) +
  geom_text(aes(label=sprintf("p=%.1e", p)), vjust=-0.4, size=3) +
  scale_fill_manual(values=c("competitive"=NNblue,"self-contained"=PPred), name="Test type") +
  labs(x=NULL, y="-log10(p)", title="B. STAT3 regulon (371 targets): three tests agree",
       subtitle="All UP in lesional; CAMERA p >> GSEA p because it discounts co-regulation") +
  ylim(0, max(db$negp)*1.15) +
  theme_bw(base_size=10) + theme(legend.position="right",
       axis.text.x=element_text(size=7.5), plot.subtitle=element_text(size=7.5))

fig18 <- pA / pB + plot_layout(heights=c(1.25,1))
ggsave("fig18_camera_roast.png", fig18, width=9, height=9, dpi=150)