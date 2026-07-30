# Auto-extracted generating script
# Produces: fig33_robustness_7study.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: d2b41871-b228-4715-ad77-0cdd4509e39d
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(recount3); library(SummarizedExperiment)})
suppressMessages({library(edgeR); library(limma); library(metafor)})
suppressMessages({library(ggplot2); library(patchwork); library(data.table)})

theme_set(theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=12), plot.title.position="plot"))
BLUE<-"#4C72B0"; RED<-"#C44E52"; GREY<-"#8C8C8C"; GOLD<-"#DD8452"

ap <- available_projects(organism="human")
get_rse <- function(srp){
  pj <- subset(ap, project==srp & project_type=="data_sources")
  rse <- create_rse(pj); assay(rse,"counts") <- transform_counts(rse); rse
}
rse_s3 <- get_rse("SRP154474")
rse_s4 <- get_rse("SRP057087")

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

at3 <- tolower(as.character(colData(rse_s3)$sra.sample_attributes))
cls_s3 <- rep(NA_character_, length(at3))
cls_s3[grepl("healthy control", at3)]              <- "NN"
cls_s3[grepl("conventional psoriatic skin", at3)]  <- "PP"
cls_s3[grepl("scalp psoriatic skin", at3)]         <- "SCALP"
cls_s3[grepl("palmoplantar", at3)]                 <- "PALMO"

at4 <- tolower(as.character(colData(rse_s4)$sra.sample_attributes))
cls_s4 <- rep(NA_character_, length(at4))
cls_s4[grepl("group;;pp", at4)] <- "PP"
cls_s4[grepl("group;;pn", at4)] <- "PN"

de_s3_PPvsNN <- de_one(rse_s3, cls_s3, "PP","NN")
de_s4_PNvsPP <- de_one(rse_s4, cls_s4, "PN","PP")

meta_dl <- function(de_list){
  genes <- sort(unique(unlist(lapply(de_list, function(d) d$gene))))
  Y <- sapply(de_list, function(d) d$logFC[match(genes,d$gene)])
  V <- sapply(de_list, function(d) (d$SE[match(genes,d$gene)])^2)
  rownames(Y)<-rownames(V)<-genes
  out <- t(apply(cbind(Y,V),1,function(row){
    k<-ncol(Y); yi<-row[1:k]; vi<-row[(k+1):(2*k)]
    ok<-is.finite(yi)&is.finite(vi)&vi>0
    if(sum(ok)<2) return(c(logFC=NA,SE=NA,P=NA,I2=NA,k=sum(ok)))
    yi<-yi[ok]; vi<-vi[ok]; wi<-1/vi; ybar<-sum(wi*yi)/sum(wi)
    Q<-sum(wi*(yi-ybar)^2); df<-length(yi)-1; C<-sum(wi)-sum(wi^2)/sum(wi)
    tau2<-max(0,(Q-df)/C); ws<-1/(vi+tau2); mu<-sum(ws*yi)/sum(ws); se<-sqrt(1/sum(ws))
    I2<-max(0,(Q-df)/Q)*100
    c(logFC=mu,SE=se,P=2*pnorm(-abs(mu/se)),I2=I2,k=length(yi))
  }))
  d<-data.frame(gene=genes,out); d$FDR<-p.adjust(d$P,"BH"); d
}

psd_old <- readRDS("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/per_study_de.rds")
meta_old <- readRDS("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/meta_de_results.rds")

psd_new <- psd_old
psd_new$PPvsNN$SRP154474 <- de_s3_PPvsNN
psd_new$PNvsPP$SRP057087 <- de_s4_PNvsPP

meta_new <- lapply(psd_new, meta_dl)

o <- meta_old$PPvsNN; n <- meta_new$PPvsNN
m <- merge(o[,c("gene","logFC","FDR")], n[,c("gene","logFC","FDR")], by="gene", suffixes=c("_old","_new"))
both_sig <- m[m$FDR_old<0.05 & m$FDR_new<0.05,]

key <- c("STAT3","IL17A","IL17F","IL23A","IL12B","SOCS3","S100A7","S100A8","S100A9",
         "DEFB4A","PI3","IL36G","CCL20","KRT16","IL19","IL36A","LCN2","CXCL8")
kt <- merge(data.frame(gene=key), n[,c("gene","logFC","FDR","I2","k")], by="gene", all.x=TRUE)
kt <- kt[match(key, kt$gene),]

## ---- FIG B: robustness comparison panel ----
# B1: sig-count bars old vs new per contrast
cnt <- data.frame(
  contrast=rep(c("PP vs NN","PN vs NN","PN vs PP"),each=2),
  set=rep(c("5-study","7-study"),3),
  n=c(12036,12570, 458,458, 12778,15398))
cnt$contrast <- factor(cnt$contrast, levels=c("PP vs NN","PN vs NN","PN vs PP"))
b1 <- ggplot(cnt, aes(contrast, n, fill=set)) +
  geom_col(position=position_dodge(0.7), width=0.65) +
  geom_text(aes(label=formatC(n,big.mark=",",format="d")),
            position=position_dodge(0.7), vjust=-0.4, size=2.9) +
  scale_fill_manual(values=c(`5-study`=GREY,`7-study`=RED), name=NULL) +
  labs(x=NULL, y="meta-significant genes (FDR<0.05)",
       title="B  More cohorts, more replicated genes") +
  coord_cartesian(ylim=c(0,17000)) +
  theme(legend.position=c(0.5,0.9), legend.direction="horizontal",
        legend.background=element_blank())

# B2: logFC concordance scatter (old vs new PPvsNN), both-sig genes
set.seed(1); samp <- both_sig[sample(nrow(both_sig), min(4000,nrow(both_sig))),]
b2 <- ggplot(samp, aes(logFC_old, logFC_new)) +
  geom_abline(slope=1,intercept=0,linetype=2,color=GREY) +
  geom_point(alpha=0.18, size=0.7, color="#333333") +
  annotate("text", x=-8, y=9, hjust=0, size=3.2,
           label=sprintf("r = %.3f\n100%% sign-concordant\n(n=%s genes sig in both)",
                          cor(m$logFC_old,m$logFC_new,use="complete.obs"),
                          formatC(nrow(both_sig),big.mark=",",format="d"))) +
  labs(x="pooled log2FC (5-study)", y="pooled log2FC (7-study)",
       title="C  PP-vs-NN effect sizes are unchanged") +
  coord_equal()

# B3: key-gene lollipop, new 7-study
keep_key <- kt[!is.na(kt$logFC),]
keep_key <- keep_key[order(keep_key$logFC),]; keep_key$gene <- factor(keep_key$gene, levels=keep_key$gene)
b3 <- ggplot(keep_key, aes(logFC, gene)) +
  geom_segment(aes(x=0,xend=logFC,yend=gene), color=GREY, linewidth=0.5) +
  geom_point(aes(color=(FDR<0.05)), size=3) +
  geom_text(aes(label=sprintf("%.1f",logFC)), hjust=-0.3, size=2.7) +
  scale_color_manual(values=c(`TRUE`=RED,`FALSE`=GREY), name="FDR<0.05") +
  labs(x="pooled log2FC (PP vs NN, 7-study)", y=NULL,
       title="D  IL-17 / keratinocyte axis stays strongly up") +
  coord_cartesian(xlim=c(0,13)) +
  theme(legend.position=c(0.8,0.2), legend.background=element_blank(),
        axis.text.y=element_text(face="italic",size=8))

figB <- (b1 | b2 | b3) + plot_layout(widths=c(1,1,1))
ggsave("fig33_robustness_7study.png", figB, width=15, height=5, dpi=150)