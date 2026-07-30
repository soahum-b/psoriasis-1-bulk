# Auto-extracted generating script
# Produces: fig23_perstudy_volcano.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: ddaf47ef-2f99-49e7-9e02-d4e45ebad803
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(edgeR)
  library(limma)
  library(SummarizedExperiment)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/")

cls_all <- readRDS("/tmp/cls_all.rds")
setDT(cls_all)
studies <- c("SRP035988","SRP165679","SRP076982","SRP126422","SRP016583")

rses <- readRDS("/tmp/rses5.rds")

run_de <- function(rse, srp, g1, g2) {
  cmap <- cls_all[cls_all$srp==srp]
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

de_PPvsNN <- list()
for (s in names(rses)) {
  r <- run_de(rses[[s]], s, "PP", "NN")
  if (!is.null(r)) { de_PPvsNN[[s]] <- r }
}

NNblue <- "#4C72B0"; PPred <- "#C44E52"; grey <- "#B0B0B0"

dfl <- rbindlist(lapply(names(de_PPvsNN), function(s){
  x <- as.data.table(de_PPvsNN[[s]]); x[, study:=s]; x
}))
dfl[, sig := adjP<0.05 & abs(logFC)>1]
dfl[, neglogP := -log10(P)]
dfl[neglogP>50, neglogP:=50]
labs <- dfl[, .(lab=sprintf("%s (%d PP vs %d NN)", study[1], n1[1], n2[1])), by=study]
dfl[, panel := labs$lab[match(study, labs$study)]]

hi <- c("STAT3","SOCS3","S100A7","DEFB4A","IL17A")
dfl[, showlab := ifelse(gene %in% hi & sig, gene, "")]

p <- ggplot(dfl, aes(logFC, neglogP)) +
  geom_point(aes(color=sig), size=0.5, alpha=0.4) +
  geom_point(data=dfl[gene %in% hi], aes(logFC,neglogP), color="black", size=1.4) +
  geom_text_repel(aes(label=showlab), size=2.4, max.overlaps=30,
                           min.segment.length=0, segment.size=0.3) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color=grey, linewidth=0.3) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color=grey, linewidth=0.3) +
  scale_color_manual(values=c("FALSE"=grey,"TRUE"=PPred), guide="none") +
  facet_wrap(~panel, nrow=1, scales="free_x") +
  labs(x="log2 fold-change (lesional vs healthy)", y="-log10 P (capped at 50)",
       title="Figure 23. Per-study lesional-vs-healthy (PP vs NN) signature",
       subtitle="STAT3 and core psoriasis genes (black) up in every study; red = FDR<0.05 & |logFC|>1") +
  theme_bw(base_size=9) +
  theme(plot.subtitle=element_text(size=7.5), strip.text=element_text(size=7))
ggsave("fig23_perstudy_volcano.png", p, width=11, height=4, dpi=150)