# Auto-extracted generating script
# Produces: fig14_stat3_isoform.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds
# Source artifact version: ba9bd6d6-df30-447c-93a3-901b08291c02
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(data.table)
  library(httr)
  library(jsonlite)
  library(ggplot2)
  library(patchwork)
})

# --- fetch junction-level RSE ---
ap <- available_projects()
proj_info <- subset(ap, project == "SRP035988" & project_type == "data_sources")
rse_jxn <- create_rse(proj_info, type = "jxn")

# --- re-derive NN/PP groups from sample_attributes ---
attr_raw <- colData(rse_jxn)$sra.sample_attributes
grp <- ifelse(grepl("normal", attr_raw, ignore.case=TRUE), "NN",
       ifelse(grepl("Psoriasis", attr_raw, ignore.case=TRUE), "PP", NA))

# --- STAT3 locus, GRCh38 minus strand ---
gr <- rowRanges(rse_jxn)
stat3_win <- GRanges("chr17", IRanges(42313000, 42389000))
hits <- which(overlapsAny(gr, stat3_win) & as.character(strand(gr)) == "-")

# --- locate the two distinguishing junctions ---
i_alpha <- which(start(gr)==42316902 & end(gr)==42317181)
i_beta  <- which(start(gr)==42316852 & end(gr)==42317181)

ca <- as.numeric(assay(rse_jxn)[i_alpha, ])
cb <- as.numeric(assay(rse_jxn)[i_beta, ])

dt <- data.table(sample=colnames(rse_jxn), grp=grp, a=ca, b=cb)
dt[, depth := a+b][, psi_beta := b/depth]

col_nn <- "#4C72B0"; col_pp <- "#C44E52"; cols <- c(NN=col_nn, PP=col_pp)
dt$grp <- factor(dt$grp, levels=c("NN","PP"))
pv <- wilcox.test(psi_beta ~ grp, data=dt)$p.value

# Panel A: PSI_beta by group
pA <- ggplot(dt, aes(grp, psi_beta*100, fill=grp)) +
  geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.55) +
  geom_jitter(aes(color=grp), width=0.12, size=1.4, alpha=0.7) +
  scale_fill_manual(values=cols, guide="none") + scale_color_manual(values=cols, guide="none") +
  annotate("text", x=1.5, y=max(dt$psi_beta*100)*1.02,
           label=sprintf("Wilcoxon p = %.3f", pv), size=3.4) +
  labs(title="STAT3-beta inclusion rises in lesional skin",
       x=NULL, y="PSI beta = beta / (alpha + beta)  (%)") +
  theme_bw(base_size=12) + theme(plot.title=element_text(size=11))

# Panel B: absolute alpha and beta junction counts (log10), both groups
long <- rbind(data.table(grp=dt$grp, iso="STAT3 alpha", cnt=dt$a),
              data.table(grp=dt$grp, iso="STAT3 beta",  cnt=dt$b))
long$iso <- factor(long$iso, levels=c("STAT3 alpha","STAT3 beta"))
pB <- ggplot(long, aes(iso, cnt+1, fill=grp)) +
  geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.55, position=position_dodge(0.7)) +
  scale_y_log10() + scale_fill_manual(values=cols, name=NULL) +
  labs(title="Both isoforms up in lesional skin; alpha dominant",
       x=NULL, y="junction reads + 1 (log10)") +
  theme_bw(base_size=12) + theme(plot.title=element_text(size=11), legend.position="top")

fig <- pA | pB
ggsave("fig14_stat3_isoform.png", fig, width=10, height=4.6, dpi=200)