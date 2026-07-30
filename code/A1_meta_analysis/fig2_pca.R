# Auto-extracted generating script
# Produces: fig2_pca.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 7d1e4543-64b0-4a29-aa39-368b74d849fe
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({ library(recount3); library(edgeR); library(dplyr); library(ggplot2); library(tidyr) })

human_projects <- available_projects()
proj_info <- subset(human_projects, project == "SRP035988" & project_type == "data_sources")
rse_gene <- create_rse(proj_info)
assay(rse_gene, "counts") <- transform_counts(rse_gene)

colData(rse_gene)$group <- NA
colData(rse_gene)$group[grepl("normal_skin",   colData(rse_gene)$sra.sample_attributes, ignore.case=TRUE)] <- "NN"
colData(rse_gene)$group[grepl("Psoriasis_skin",colData(rse_gene)$sra.sample_attributes, ignore.case=TRUE)] <- "PP"
grp <- factor(colData(rse_gene)$group, levels=c("NN","PP"))

myCounts <- assay(rse_gene, "counts")
gs <- rowData(rse_gene)$gene_name
ok <- !is.na(gs) & gs != ""
myCounts <- rowsum(myCounts[ok,], group = gs[ok])

dge <- DGEList(counts = round(myCounts), group = grp)
keep_fbe <- filterByExpr(dge, group = grp)
dge.filt <- dge[keep_fbe, , keep.lib.sizes = FALSE]
dge.filt.norm <- calcNormFactors(dge.filt, method = "TMM")

log2cpm.after  <- cpm(dge.filt.norm, log=TRUE, normalized.lib.sizes=TRUE)

pal <- c(NN="#4C72B0", PP="#C44E52")

pca <- prcomp(t(log2cpm.after), scale.=FALSE)
pcv <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
pcadf <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=grp)
p2 <- ggplot(pcadf, aes(PC1, PC2, color=group)) +
  geom_point(size=3, alpha=0.85) +
  scale_color_manual(values=pal, labels=c(NN="Normal (NN)", PP="Lesional (PP)")) +
  labs(title="Lesional and normal skin separate on PC1",
       subtitle="PCA of TMM-normalized log2-CPM (24,528 genes, 178 samples)",
       x=paste0("PC1 (",pcv[1],"%)"), y=paste0("PC2 (",pcv[2],"%)"), color=NULL) +
  theme_bw(base_size=12) +
  theme(legend.position="top", plot.title=element_text(face="bold"))
ggsave("fig2_pca.png", p2, width=7, height=6, dpi=300)