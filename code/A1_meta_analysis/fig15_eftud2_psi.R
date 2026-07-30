# Auto-extracted generating script
# Produces: fig15_eftud2_psi.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): dge_filt_norm.rds
# Source artifact version: c1389080-7d31-45c9-938c-10a81b9ad8e8
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(edgeR)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(httr)
  library(jsonlite)
})

# --- fetch junction-level RSE ---
ap <- available_projects()
proj_info <- subset(ap, project == "SRP035988" & project_type == "data_sources")
rse_jxn <- create_rse(proj_info, type = "jxn")

# --- re-derive NN/PP groups ---
attr_raw <- colData(rse_jxn)$sra.sample_attributes
grp <- ifelse(grepl("normal", attr_raw, ignore.case=TRUE), "NN",
       ifelse(grepl("Psoriasis", attr_raw, ignore.case=TRUE), "PP", NA))

# --- locate STAT3-locus junctions ---
gr <- rowRanges(rse_jxn)
stat3_win <- GRanges("chr17", IRanges(42313000, 42389000))
hits <- which(overlapsAny(gr, stat3_win) & as.character(strand(gr)) == "-")

# --- alpha and beta distinguishing junctions ---
i_alpha <- which(start(gr)==42316902 & end(gr)==42317181)
i_beta  <- which(start(gr)==42316852 & end(gr)==42317181)

ca <- as.numeric(assay(rse_jxn)[i_alpha, ])
cb <- as.numeric(assay(rse_jxn)[i_beta, ])

dt <- data.table(sample=colnames(rse_jxn), grp=grp, a=ca, b=cb)
dt[, depth := a+b][, psi_beta := b/depth]

# --- load filtered+normalized gene-level DGEList checkpoint ---
dge <- readRDS("dge_filt_norm.rds")
logcpm <- cpm(dge, log=TRUE)

# --- build merged per-sample table ---
m <- copy(dt)
setkey(m, sample)
m[, EFTUD2 := logcpm["EFTUD2", sample]]
m[, STAT3  := logcpm["STAT3",  sample]]

m$grp <- factor(m$grp, levels=c("NN","PP"))

col_nn <- "#4C72B0"; col_pp <- "#C44E52"

rho_pp <- cor.test(m[grp=="PP"]$psi_beta, m[grp=="PP"]$EFTUD2, method="spearman")$estimate
rho_nn <- cor.test(m[grp=="NN"]$psi_beta, m[grp=="NN"]$EFTUD2, method="spearman")$estimate

# Panel A: scatter PSI_beta vs EFTUD2, colored by group, per-group trend lines
pA <- ggplot(m, aes(EFTUD2, psi_beta*100, color=grp)) +
  geom_point(size=1.5, alpha=0.7) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.8, alpha=0.15) +
  scale_color_manual(values=c(NN=col_nn, PP=col_pp), name=NULL) +
  labs(title="Within-group: higher EFTUD2 tracks LOWER STAT3-beta",
       x="EFTUD2 expression (log2 CPM)", y="PSI beta (%)") +
  annotate("text", x=min(m$EFTUD2), y=22,
           label=sprintf("NN rho=%.2f  PP rho=%.2f", rho_nn, rho_pp),
           hjust=0, size=3.3) +
  theme_bw(base_size=12) + theme(plot.title=element_text(size=10.5), legend.position="top")

# Panel B: EFTUD2 is itself up in PP (the confound that flips the overall corr)
pB <- ggplot(m, aes(grp, EFTUD2, fill=grp)) +
  geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.55) +
  geom_jitter(aes(color=grp), width=0.12, size=1.3, alpha=0.65) +
  scale_fill_manual(values=c(NN=col_nn, PP=col_pp), guide="none") +
  scale_color_manual(values=c(NN=col_nn, PP=col_pp), guide="none") +
  labs(title="EFTUD2 up in lesional skin (p = 7e-25)", x=NULL, y="EFTUD2 (log2 CPM)") +
  theme_bw(base_size=12) + theme(plot.title=element_text(size=10.5))

fig <- pA | pB
ggsave("fig15_eftud2_psi.png", fig, width=10, height=4.6, dpi=200)