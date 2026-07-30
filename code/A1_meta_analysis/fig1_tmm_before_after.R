# Auto-extracted generating script
# Produces: fig1_tmm_before_after.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 8301867f-21af-402e-a8fa-188586f7573b
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

log2cpm.before <- cpm(dge.filt,      log=TRUE, normalized.lib.sizes=FALSE)
log2cpm.after  <- cpm(dge.filt.norm, log=TRUE, normalized.lib.sizes=TRUE)

set.seed(1)
mk_long <- function(m, lab){
  df <- as.data.frame(m)
  df$gene <- rownames(df)
  tidyr::pivot_longer(df, -gene, names_to="sample", values_to="log2cpm") |>
    dplyr::mutate(stage=lab)
}
long.before <- mk_long(log2cpm.before, "Before TMM (raw library size)")
long.after  <- mk_long(log2cpm.after,  "After TMM")
plotdf <- dplyr::bind_rows(long.before, long.after)
grp_map <- setNames(as.character(grp), colnames(dge.filt.norm))
plotdf$group <- grp_map[plotdf$sample]
plotdf$stage <- factor(plotdf$stage, levels=c("Before TMM (raw library size)","After TMM"))

samp_order <- colnames(dge.filt.norm)[order(grp)]
plotdf$sample <- factor(plotdf$sample, levels=samp_order)

pal <- c(NN="#4C72B0", PP="#C44E52")
p1 <- ggplot(plotdf, aes(x=sample, y=log2cpm, fill=group)) +
  geom_boxplot(outlier.shape=NA, linewidth=0.2, width=0.8) +
  facet_wrap(~stage, ncol=1) +
  scale_fill_manual(values=pal, labels=c(NN="Normal (NN)", PP="Lesional (PP)")) +
  labs(title="TMM normalization aligns sample expression distributions",
       subtitle="Each box = one sample's log2-CPM across 24,528 genes; medians should line up after TMM",
       x="Samples (n=178, ordered by group)", y=expression(log[2]~CPM), fill=NULL) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid.major.x=element_blank(), legend.position="top",
        plot.title=element_text(face="bold"))
ggsave("fig1_tmm_before_after.png", p1, width=10, height=7, dpi=300)