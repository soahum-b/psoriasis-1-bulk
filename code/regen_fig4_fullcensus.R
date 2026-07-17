#!/usr/bin/env Rscript
# regen_fig4_fullcensus.R — regenerate WHITEPAPER §4 figures from the 89k full-census object.
# Replaces the 20k-backbone figures. Outputs to figures_full/.
suppressMessages({library(Seurat); library(ggplot2); library(patchwork); library(Matrix)})
OUT <- "figures_full"; dir.create(OUT, showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))
pal <- c("Scissor-"="#2c7fb8","Background"="#bdbdbd","Scissor+"="#d7301f")
tierpal <- c("NN"="#4575b4","PN"="#fdae61","PP"="#d73027")
msg <- function(...) cat(sprintf(...), "\n")

## ---- load object + result tables ----
so   <- readRDS("results_full/reference_scissor_full.rds")
tun  <- readRDS("results_full/scissor_tuning.rds")
rel  <- readRDS("results_full/reliability_test.rds")
pn   <- readRDS("results_full/permutation_null.rds")
de   <- read.csv("results_full/gradient_program_DE_BH.csv", stringsAsFactors = FALSE)
enr  <- read.csv("results_full/celltype_enrichment_BH.csv", stringsAsFactors = FALSE)
msg("cells=%d  scissor levels: %s", ncol(so), paste(levels(so$scissor), collapse="/"))
if (!"umap" %in% Reductions(so)) { msg("no umap -> computing"); so <- RunUMAP(so, dims=1:20, verbose=FALSE) }

## ---- Fig A: alpha tuning ----
tt <- tun$tuning
pA <- ggplot(tt, aes(alpha, selected_frac*100)) +
  geom_hline(yintercept=20, linetype=2, color="grey50") +
  geom_line(color="#377eb8") + geom_point(size=3, color="#377eb8") +
  geom_point(data=tt[which.min(abs(tt$selected_frac-0.15)),], size=5, shape=21, fill="#d7301f") +
  labs(x="alpha (graph-smoothing strength)", y="% cells selected",
       title="Full-census alpha tuning", subtitle="chosen alpha=0.20, 14.67% selected (<20% cutoff)")
ggsave(file.path(OUT,"fig_alpha_tuning.png"), pA, width=6, height=4.5, dpi=150)

## ---- Fig B: Scissor UMAP ----
df <- data.frame(so@reductions$umap@cell.embeddings, scissor=so$scissor)
names(df)[1:2] <- c("UMAP1","UMAP2")
df <- df[order(df$scissor),]  # draw selected on top
pB <- ggplot(df, aes(UMAP1, UMAP2, color=scissor)) +
  geom_point(size=0.15, alpha=0.6) + scale_color_manual(values=pal) +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1))) +
  labs(title="Scissor selection on the full-census reference (89,058 cells)", color=NULL) +
  theme(legend.position="bottom")
ggsave(file.path(OUT,"fig_scissor_umap.png"), pB, width=7, height=7, dpi=150)

## ---- Fig C: composition (celltype enrichment + tier fractions) ----
enr$celltype <- factor(enr$celltype, levels=enr$celltype[order(enr$OR)])
pC1 <- ggplot(enr, aes(celltype, OR, fill=OR>1)) + geom_col() + coord_flip() +
  geom_hline(yintercept=1, linetype=2) +
  scale_fill_manual(values=c("TRUE"="#d7301f","FALSE"="#2c7fb8"), guide="none") +
  labs(x=NULL, y="odds ratio in Scissor+", title="Cell-type enrichment (full census)",
       subtitle="Endothelial OR 5.24, top lineage")
comp <- as.data.frame(prop.table(table(so$scissor, so$condition), 1))
names(comp) <- c("scissor","tier","frac")
pC2 <- ggplot(comp, aes(scissor, frac, fill=tier)) + geom_col() +
  scale_fill_manual(values=tierpal) + labs(x=NULL, y="fraction of cells", fill="tier",
       title="Tier composition by Scissor class", subtitle="Scissor+ shifts toward PN/PP")
ggsave(file.path(OUT,"fig_scissor_composition.png"), pC1|pC2, width=11, height=4.6, dpi=150)

## ---- Fig D: significance controls ----
d1 <- data.frame(mse=rel$background)
pD1 <- ggplot(d1, aes(mse)) + geom_histogram(bins=30, fill="grey70", color="white") +
  geom_vline(xintercept=rel$statistic, color="#d7301f", linewidth=1.2) +
  labs(x="mean CV-MSE", y="permutations",
       title="Reliability test", subtitle=sprintf("real %.3f vs null; p=%.3f", rel$statistic, rel$p))
d2 <- data.frame(gap=pn$null_dir)
pD2 <- ggplot(d2, aes(gap)) + geom_histogram(bins=25, fill="grey70", color="white") +
  geom_vline(xintercept=pn$real_gap, color="#d7301f", linewidth=1.2) +
  labs(x="pos - neg mean tier gap", y="permutations",
       title="Selection permutation null",
       subtitle=sprintf("real %.3f vs null mean %.3f", pn$real_gap, mean(pn$null_dir, na.rm=TRUE)))
ggsave(file.path(OUT,"fig_significance.png"), pD1|pD2, width=11, height=4.6, dpi=150)

## ---- Fig E: gradient program volcano (BH-FDR) ----
de$sig <- de$fdr_BH < 0.05
de$neglog <- -log10(pmax(de$fdr_BH, 1e-300))
lab <- de[de$gene %in% c("CCL14","ACKR1","RAMP3","PLVAP","APLNR","CYTL1","SPNS2","PECAM1","VWF","STAT3"),]
pE <- ggplot(de, aes(avg_log2FC, neglog, color=sig)) +
  geom_point(size=0.5, alpha=0.5) +
  scale_color_manual(values=c("TRUE"="#d7301f","FALSE"="grey70"), guide="none") +
  ggrepel::geom_text_repel(data=lab, aes(label=gene), color="black", size=3, max.overlaps=20) +
  labs(x="avg log2FC (Scissor+ vs background)", y="-log10 BH-FDR",
       title="Gradient-tracking gene program (full census, BH-FDR)",
       subtitle=sprintf("%d/%d genes at BH q<0.05; vascular-led", sum(de$sig), nrow(de)))
ggsave(file.path(OUT,"fig_gradient_program.png"), pE, width=7, height=5.5, dpi=150)

## ---- Fig F: STAT3 across Scissor classes (now n.s.) ----
st3 <- FetchData(so, vars="STAT3", layer="data")[,1]
dfF <- data.frame(expr=st3, scissor=so$scissor)
s <- de[de$gene=="STAT3",]
sub <- if(nrow(s)) sprintf("log2FC %.2f, raw p=%.2f, BH q=%.2f -> n.s. at full census", s$avg_log2FC, s$p_val, s$fdr_BH) else "STAT3"
pF <- ggplot(dfF, aes(scissor, expr, fill=scissor)) +
  geom_violin(scale="width", trim=TRUE) +
  geom_boxplot(width=0.12, outlier.shape=NA, fill="white") +
  scale_fill_manual(values=pal, guide="none") +
  labs(x=NULL, y="STAT3 log-normalized expression",
       title="STAT3 across Scissor classes (full census)", subtitle=sub)
ggsave(file.path(OUT,"fig_stat3.png"), pF, width=6, height=5, dpi=150)

msg("DONE. figures in %s/", OUT)
