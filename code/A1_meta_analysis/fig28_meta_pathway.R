# Auto-extracted generating script
# Produces: fig28_meta_pathway.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 4f5dfa41-c806-4789-a8c8-18f9b1f52bb6
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(edgeR)
library(limma)
library(msigdbr)
library(recount3)
library(SummarizedExperiment)
library(data.table)
library(ggplot2)
library(ggrepel)

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/")

cls_all_raw <- readRDS("/tmp/cls_all.rds")
setDT(cls_all_raw)

rses <- readRDS("/tmp/rses5.rds")

H <- msigdbr(species="Homo sapiens", category="H")
Hlist <- split(H$gene_symbol, H$gs_name)

camera_study <- function(rse, srp, g1="PP", g2="NN") {
  cmap <- cls_all_raw[cls_all_raw$srp==srp]
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

studies_ppnn <- c("SRP035988","SRP165679","SRP126422")
cam_all <- rbindlist(lapply(studies_ppnn, function(s) camera_study(rses[[s]], s)), fill=TRUE)

inv <- fread("study_inventory.csv")

cam_all[, Pfloor := pmax(PValue, 1e-300)]
cam_all[, zraw := qnorm(1 - pmin(Pfloor,0.9999999)/2)]
cam_all[, z := ifelse(Direction=="Up",1,-1) * pmin(zraw, 38)]
cam_all[!is.finite(z), z := ifelse(Direction=="Up",1,-1)*38]
cam_all[, w := sqrt(inv$total[match(study, inv$srp)])]
pooled <- cam_all[, .(k=.N, Zc=sum(w*z)/sqrt(sum(w^2)),
                      nUp=sum(Direction=="Up"), meanNGenes=round(mean(NGenes))), by=pathway]
pooled[, p_comb := 2*pnorm(-abs(pmin(Zc,38)))][, FDR := p.adjust(p_comb,"BH")]
pooled <- pooled[order(-Zc)]

NNblue<-"#4C72B0"; PPred<-"#C44E52"; green<-"#55A868"; grey<-"#B0B0B0"

top <- pooled[order(-Zc)][1:20]
top[, pw := sub("HALLMARK_","",pathway)]
top[, pw := gsub("_"," ", pw)]
top[, is_stat3 := grepl("IL6 JAK STAT3", pw)]
top[, Zc_disp := pmin(Zc, 38)]
top <- top[order(Zc_disp)]
top[, pw := factor(pw, levels=pw)]

p <- ggplot(top, aes(Zc_disp, pw, fill=is_stat3)) +
  geom_col(width=0.72) +
  geom_text(aes(label=sprintf("k=%d", k)), hjust=-0.15, size=2.6, color="#555555") +
  scale_fill_manual(values=c("FALSE"=NNblue,"TRUE"=PPred), guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x="pooled Stouffer Z (capped at 38; all pathways up in lesional)", y=NULL,
       title="Figure 28. Pooled Hallmark pathway enrichment (PP vs NN, k=3 studies)",
       subtitle="CAMERA per study, Stouffer-combined. IL6-JAK-STAT3 highlighted (red). Interferon + proliferation dominate.") +
  theme_bw(base_size=9) +
  theme(plot.subtitle=element_text(size=6.8), axis.text.y=element_text(size=7.2))
ggsave("fig28_meta_pathway.png", p, width=9.5, height=5, dpi=150)