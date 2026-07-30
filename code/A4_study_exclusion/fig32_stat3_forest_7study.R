# Auto-extracted generating script
# Produces: fig32_stat3_forest_7study.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 59d0881b-b88c-4cea-a384-b97e0c6416b2
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(metafor)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

theme_set(theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=12), plot.title.position="plot"))
BLUE<-"#4C72B0"; RED<-"#C44E52"; GREY<-"#8C8C8C"; GOLD<-"#DD8452"

# --- Load recount3 data ---
ap <- available_projects(organism="human")
get_rse <- function(srp){
  pj <- subset(ap, project==srp & project_type=="data_sources")
  rse <- create_rse(pj); assay(rse,"counts") <- transform_counts(rse); rse
}
rse_s3 <- get_rse("SRP154474")

# --- S3 classification ---
at3 <- tolower(as.character(colData(rse_s3)$sra.sample_attributes))
cls_s3 <- rep(NA_character_, length(at3))
cls_s3[grepl("healthy control", at3)]              <- "NN"
cls_s3[grepl("conventional psoriatic skin", at3)]  <- "PP"
cls_s3[grepl("scalp psoriatic skin", at3)]         <- "SCALP"
cls_s3[grepl("palmoplantar", at3)]                 <- "PALMO"

# --- de_one function ---
de_one <- function(rse, cls, cA, cB){
  keep_s <- which(cls %in% c(cA,cB))
  if(length(unique(cls[keep_s]))<2) return(NULL)
  grp <- factor(cls[keep_s], levels=c(cB,cA))
  if(min(table(grp))<2) return(NULL)
  counts <- assay(rse,"counts")[,keep_s]
  gn <- rowData(rse)$gene_name
  dge <- DGEList(counts=counts, genes=data.frame(gene_name=gn))
  keep <- filterByExpr(dge, group=grp); dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge,"TMM")
  des <- model.matrix(~grp); v <- voom(dge,des); fit <- eBayes(lmFit(v,des))
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
  data.frame(gene=dge$genes$gene_name, logFC=tt$logFC, SE=tt$logFC/tt$t,
             t=tt$t, P=tt$P.Value, adjP=p.adjust(tt$P.Value,"BH"),
             AveExpr=tt$AveExpr, n1=sum(grp==cA), n2=sum(grp==cB))
}

de_s3_PPvsNN <- de_one(rse_s3, cls_s3, "PP","NN")

# --- Load saved 5-study per-study DE and append S3 ---
psd_old <- readRDS("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/per_study_de.rds")

psd_new <- psd_old
psd_new$PPvsNN$SRP154474 <- de_s3_PPvsNN

# --- STAT3 forest data ---
stat3_rows <- function(psd){
  do.call(rbind, lapply(names(psd$PPvsNN), function(s){
    d<-psd$PPvsNN[[s]]; r<-d[d$gene=="STAT3",][1,]
    data.frame(study=s, yi=r$logFC, sei=r$SE, n1=r$n1, n2=r$n2)
  }))
}
s3_old <- stat3_rows(psd_old)
s3_new <- stat3_rows(psd_new)

r_old  <- rma(yi=yi, sei=sei, data=s3_old, method="DL")
r_new  <- rma(yi=yi, sei=sei, data=s3_new, method="DL")
rk_old <- rma(yi=yi, sei=sei, data=s3_old, method="REML", test="knha")
rk_new <- rma(yi=yi, sei=sei, data=s3_new, method="REML", test="knha")

# --- Diamond helper ---
dia <- function(res,y,labtxt){ data.frame(
  x=c(res$ci.lb,res$b,res$ci.ub,res$b), y=c(y,y+0.28,y,y-0.28), lab=labtxt) }

# --- FIG: STAT3 4-study forest (PP vs NN) ---
lab <- c(SRP035988="SRP035988 (Li, 95v83)", SRP165679="SRP165679 (Tsoi, 28v38)",
         SRP126422="SRP126422 (biopsy, 4v4)", SRP154474="SRP154474 (conv. plaque, 8v9) NEW")
fdf <- s3_new
fdf$study_lab <- lab[fdf$study]
fdf$ci.lb <- fdf$yi-1.96*fdf$sei; fdf$ci.ub <- fdf$yi+1.96*fdf$sei
fdf$w <- 1/(fdf$sei^2); fdf$w <- fdf$w/max(fdf$w)
fdf$is_new <- grepl("NEW", fdf$study_lab)
fdf <- fdf[order(fdf$yi),]; fdf$y <- seq_len(nrow(fdf))
yD1 <- 0; yD2 <- -1.1
d_new_dl <- dia(r_new, yD1, "RE pooled (DL, 4-study)")
d_new_hk <- dia(rk_new, yD2, "RE pooled (HKSJ, 4-study)")

figA <- ggplot() +
  geom_vline(xintercept=0, linetype=2, color=GREY) +
  geom_errorbarh(data=fdf, aes(y=y, xmin=ci.lb, xmax=ci.ub, color=is_new), height=0, linewidth=0.7) +
  geom_point(data=fdf, aes(y=y, x=yi, size=w, color=is_new)) +
  geom_polygon(data=d_new_dl, aes(x=x,y=y), fill=RED, color="black", linewidth=0.3) +
  geom_polygon(data=d_new_hk, aes(x=x,y=y), fill=GOLD, color="black", linewidth=0.3) +
  geom_text(data=fdf, aes(y=y, x=yi, label=sprintf("%+.2f",yi)), vjust=-1.2, size=3) +
  annotate("text", x=r_new$b, y=yD1-0.55, size=3,
           label=sprintf("DL %+.2f [%.2f, %.2f], p=%.1e", r_new$b,r_new$ci.lb,r_new$ci.ub,r_new$pval)) +
  annotate("text", x=rk_new$b, y=yD2-0.55, size=3,
           label=sprintf("HKSJ %+.2f [%.2f, %.2f], p=%.3f", rk_new$b,rk_new$ci.lb,rk_new$ci.ub,rk_new$pval)) +
  scale_color_manual(values=c(`FALSE`="#333333",`TRUE`=RED), guide="none") +
  scale_size(range=c(2.5,7), guide="none") +
  scale_y_continuous(breaks=fdf$y, labels=fdf$study_lab,
        sec.axis=dup_axis(breaks=c(yD1,yD2), labels=c("RE pooled (DL)","RE pooled (HKSJ)"))) +
  labs(x="log2 fold-change, lesional (PP) vs healthy (NN)", y=NULL,
       title="STAT3 up-regulation holds and sharpens with a 4th cohort",
       subtitle="Adding SRP154474 moves the conservative HKSJ interval off zero (was 0.059 at 3 studies)") +
  coord_cartesian(xlim=c(-0.4,2.6), ylim=c(yD2-0.9, max(fdf$y)+0.7)) +
  theme(axis.text.y=element_text(size=9))
ggsave("fig32_stat3_forest_7study.png", figA, width=9, height=5, dpi=150)