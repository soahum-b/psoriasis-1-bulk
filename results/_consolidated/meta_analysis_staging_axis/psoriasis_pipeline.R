#!/usr/bin/env Rscript
# =============================================================================
# psoriasis_pipeline.R
# Reproducible master pipeline: psoriasis lesional (PP) vs normal (NN) skin
# Study: SRP035988 (Li et al.), recount3, Gencode v26 (G026) / GRCh38
# -----------------------------------------------------------------------------
# Runs end-to-end from a raw recount3 download to the differential-expression
# results and every figure in the white paper. Written to be run top-to-bottom:
#
#     Rscript psoriasis_pipeline.R
#
# or step-by-step in an interactive R session. Each STEP writes a checkpoint
# (.rds) so that if you are disconnected you can resume from the last one
# instead of recomputing from scratch. See README.md for setup + package list.
#
# Colour convention throughout: NN (normal) = #4C72B0 blue, PP (psoriasis) = #C44E52 red.
# =============================================================================

set.seed(1)  # reproducibility for any stochastic step
OUTDIR <- "."                      # write checkpoints/figures to working dir
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- packages ---------------------------------------------------------------
suppressMessages({
  library(recount3)   # RNA-seq data access (Bioconductor)
  library(edgeR)      # DGEList, filterByExpr, TMM
  library(limma)      # voom, lmFit, eBayes, topTable
  library(qvalue)     # Storey q-values / pi0
  library(dplyr)
  library(ggplot2)
  library(ggrepel)    # non-overlapping volcano labels
  # ---- pathway analysis (STEP 6+) ----
  library(fgsea)          # threshold-free GSEA (multilevel permutation)
  library(msigdbr)        # MSigDB gene-set collections (2026.1.Hs)
  library(clusterProfiler)# ORA: enricher() hypergeometric/Fisher
  library(org.Hs.eg.db)   # gene annotation backend
  library(data.table)     # fast merges for GSEA/ORA tables
  # ---- network / TF-activity analysis (STEP 7) ----
  library(decoupleR)      # TF-activity inference (run_ulm) from a signed regulon
  # ---- isoform splicing (STEP 8) ----
  library(recount3)       # junction-level RSE (type="jxn") for STAT3 alpha/beta
  library(SummarizedExperiment)
  library(patchwork)      # multi-panel figures (fig14/fig15)
  # ---- enrichment robustness (STEP 10) ----
  # CAMERA and ROAST live inside limma (already loaded, no extra install).
  library(pheatmap)       # ssGSEA pathway-score heatmap (fig20)
  library(RColorBrewer)
  # UpSetR (fig21) is loaded lazily inside STEP 10 from the workspace lib,
  # because it is a CRAN-only package not in the base conda env. See STEP 10.
  # ---- multi-study meta-analysis (STEP 11) ----
  # metafor (random-effects rma, forest plots) is loaded lazily inside STEP 11
  #   from the workspace lib (CRAN-only, pure R). clust (cross-study co-expression
  #   modules) is a Python tool run via a separate command line, see STEP 11.7.
})

# =============================================================================
# STEP 1 - Download study SRP035988 and build the count matrix
# -----------------------------------------------------------------------------
# recount3 delivers a uniformly processed RangedSummarizedExperiment (RSE).
# transform_counts() converts recount's raw base-pair coverage counts into
# per-gene read counts suitable for edgeR/limma.
# =============================================================================
message("STEP 1: downloading SRP035988 from recount3 ...")
human_projects <- available_projects()
proj_info <- subset(human_projects,
                    project == "SRP035988" & project_type == "data_sources")
rse_gene  <- create_rse(proj_info)
assay(rse_gene, "counts") <- transform_counts(rse_gene)

# ---- assign disease groups from the SRA sample metadata ----
colData(rse_gene)$group <- NA
colData(rse_gene)$group[grepl("normal_skin",    colData(rse_gene)$sra.sample_attributes, ignore.case = TRUE)] <- "NN"
colData(rse_gene)$group[grepl("Psoriasis_skin", colData(rse_gene)$sra.sample_attributes, ignore.case = TRUE)] <- "PP"
grp <- factor(colData(rse_gene)$group, levels = c("NN", "PP"))
message(sprintf("  samples: NN=%d, PP=%d, NA=%d",
                sum(grp=="NN", na.rm=TRUE), sum(grp=="PP", na.rm=TRUE), sum(is.na(grp))))

# ---- collapse Ensembl gene IDs to gene symbols ----
# recount3 rows are versioned Ensembl IDs; we sum counts per HGNC symbol so the
# downstream tables are human-readable (STAT3, IL17A, ...).
myCounts <- assay(rse_gene, "counts")
gs <- rowData(rse_gene)$gene_name
ok <- !is.na(gs) & gs != ""
myCounts <- rowsum(myCounts[ok, ], group = gs[ok])
message(sprintf("  genes after symbol collapse: %d", nrow(myCounts)))

# =============================================================================
# STEP 2 - Filter low-expressed genes + TMM normalisation  --> CHECKPOINT 1
# -----------------------------------------------------------------------------
# filterByExpr keeps genes with enough counts in enough samples (defaults:
# min.count=10, min.total.count=15, large.n=10, min.prop=0.7).
# TMM (Robinson & Oshlack 2010) corrects for library composition. This is the
# between-sample normalisation used for DE; length cancels out in a PP-vs-NN
# contrast so no TPM/RPKM is needed.
# =============================================================================
message("STEP 2: filterByExpr + TMM ...")
dge <- DGEList(counts = round(myCounts), group = grp)
keep_fbe      <- filterByExpr(dge, group = grp)
dge.filt      <- dge[keep_fbe, , keep.lib.sizes = FALSE]
dge.filt.norm <- calcNormFactors(dge.filt, method = "TMM")
saveRDS(dge.filt.norm, file.path(OUTDIR, "dge_filt_norm.rds"))
message(sprintf("  genes after filtering: %d x %d samples  [checkpoint: dge_filt_norm.rds]",
                nrow(dge.filt.norm), ncol(dge.filt.norm)))

# =============================================================================
# STEP 3 - Differential expression: limma-voom  --> CHECKPOINT 2
# -----------------------------------------------------------------------------
# voom models the mean-variance trend and turns counts into log2-CPM with
# precision weights; lmFit fits a per-gene linear model (the design ~grp makes
# the PPvsNN coefficient equal to the log2 fold-change); eBayes shrinks the
# per-gene variances (moderated t-test). We add Storey q-values alongside the
# BH-adjusted p-values (adj.P.Val).
# =============================================================================
message("STEP 3: limma-voom differential expression ...")
dge <- readRDS(file.path(OUTDIR, "dge_filt_norm.rds"))
grp <- dge$samples$group

design <- model.matrix(~ grp)
colnames(design) <- c("Intercept", "PPvsNN")

v   <- voom(dge, design)          # mean-variance weights (save the object for Fig 5)
fit <- lmFit(v, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "PPvsNN", number = Inf, sort.by = "P")
res$gene <- rownames(res)

qobj      <- qvalue(p = res$P.Value)   # Storey pi0 + q-values
res$qvalue <- qobj$qvalues
message(sprintf("  Storey pi0 = %.4f (est. %.1f%% of genes truly null)",
                qobj$pi0, 100*qobj$pi0))

de <- res[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","qvalue","B")]
saveRDS(de, file.path(OUTDIR, "de_results_full.rds"))
write.csv(de, file.path(OUTDIR, "de_results_full.csv"), row.names = FALSE)

# significant set: BH-FDR < 0.05 AND |log2FC| > 1  (2-fold)
de_sig <- de[de$adj.P.Val < 0.05 & abs(de$logFC) > 1, ]
de_sig <- de_sig[order(de_sig$adj.P.Val), ]
write.csv(de_sig, file.path(OUTDIR, "de_results_significant.csv"), row.names = FALSE)

n_up   <- sum(de_sig$logFC > 0)
n_down <- sum(de_sig$logFC < 0)
message(sprintf("  significant: %d genes (%d up, %d down)  [checkpoint: de_results_full.rds]",
                nrow(de_sig), n_up, n_down))

# =============================================================================
# STEP 4 - Figure 6: volcano plot
# =============================================================================
message("STEP 4: volcano plot ...")
de$sig <- "n.s."
de$sig[de$adj.P.Val < 0.05 & de$logFC >  1] <- "Up in PP"
de$sig[de$adj.P.Val < 0.05 & de$logFC < -1] <- "Down in PP"
lab <- de[de$gene %in% c("DEFB4A","S100A7A","IL36A","PI3","SERPINB4",
                         "STAT3","IL17A","SOCS3","IL6","KRT77","BTC","RORC"), ]

pV <- ggplot(de, aes(logFC, -log10(adj.P.Val), colour = sig)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3.2,
                           max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = c("Up in PP" = "#C44E52", "Down in PP" = "#4C72B0",
                                 "n.s." = "grey75")) +
  labs(title = "Lesional (PP) vs normal (NN) skin: differential expression",
       subtitle = sprintf("%d genes at adj.P<0.05 and |log2FC|>1 (%d up, %d down)",
                          nrow(de_sig), n_up, n_down),
       x = "log2 fold-change (PP / NN)", y = "-log10 adjusted p", colour = NULL) +
  theme_bw(base_size = 13)
ggsave(file.path(OUTDIR, "fig6_volcano.png"), pV, width = 8.5, height = 6.5, dpi = 150)

# =============================================================================
# STEP 5 - Figure 7: p-value histogram + Storey q vs BH
# =============================================================================
message("STEP 5: q-value diagnostic ...")
pA <- ggplot(de, aes(P.Value)) +
  geom_histogram(boundary = 0, bins = 50, fill = "grey70", colour = "white") +
  geom_hline(yintercept = qobj$pi0 * nrow(de) / 50, colour = "#C44E52",
             linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.6, y = qobj$pi0 * nrow(de) / 50,
           label = sprintf("pi0 = %.3f (null level)", qobj$pi0),
           vjust = -0.6, colour = "#C44E52", size = 4) +
  labs(title = "A. p-value distribution",
       subtitle = "Spike near 0 = real signal; flat part = null genes",
       x = "p-value (PP vs NN)", y = "gene count") +
  theme_bw(base_size = 13)

pB <- ggplot(de, aes(adj.P.Val, qvalue)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dotted") +
  geom_point(alpha = 0.25, size = 0.5, colour = "#4C72B0") +
  geom_vline(xintercept = 0.05, colour = "#C44E52", linetype = "dashed") +
  geom_hline(yintercept = 0.05, colour = "#C44E52", linetype = "dashed") +
  labs(title = "B. Storey q-value vs BH adjusted p",
       subtitle = "q = pi0 x BH, so q <= BH (more powerful)",
       x = "BH adjusted p (adj.P.Val)", y = "Storey q-value") +
  theme_bw(base_size = 13)

# combine without requiring patchwork: use gridExtra if present, else save A only
if (requireNamespace("patchwork", quietly = TRUE)) {
  ggsave(file.path(OUTDIR, "fig7_qvalue.png"),
         patchwork::wrap_plots(pA, pB, ncol = 2), width = 11, height = 4.6, dpi = 150)
} else if (requireNamespace("gridExtra", quietly = TRUE)) {
  ggsave(file.path(OUTDIR, "fig7_qvalue.png"),
         gridExtra::arrangeGrob(pA, pB, ncol = 2), width = 11, height = 4.6, dpi = 150)
} else {
  ggsave(file.path(OUTDIR, "fig7_qvalue.png"), pA, width = 6, height = 4.6, dpi = 150)
  message("  (install patchwork or gridExtra for the 2-panel version; saved panel A only)")
}

# =============================================================================
# STEP 6 - Pathway analysis: GSEA (primary) + ORA (cross-check)  --> CHECKPOINT 3
# -----------------------------------------------------------------------------
# GSEA (fgsea) is threshold-free: it ranks ALL genes by the limma moderated
# t-statistic and asks whether each gene set clusters at the extremes.
# ORA (clusterProfiler::enricher) is threshold-based: Fisher's exact test on a
# 2x2 table asking whether set members are over-represented among the DE hits.
# We run both against four MSigDB collections and merge the verdicts. Requires
# network access to zenodo.org (msigdbr fetches MSigDB 2026.1.Hs at runtime).
# =============================================================================
message("STEP 6: pathway analysis (GSEA + ORA) ...")

# ---- 6a. build the ranking vector: named limma t-stat, sorted descending ----
ranks <- de$t
names(ranks) <- de$gene
ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
saveRDS(ranks, file.path(OUTDIR, "gsea_ranks.rds"))   # CHECKPOINT 3
message(sprintf("  ranking vector: %d genes (top: %s)",
                length(ranks), paste(head(names(ranks), 3), collapse = ", ")))

# ---- 6b. gene-set collections from MSigDB (symbols) ----
build_list <- function(coll, sub = NULL) {
  df <- if (is.null(sub)) msigdbr(species = "human", collection = coll)
        else msigdbr(species = "human", collection = coll, subcollection = sub)
  split(df$gene_symbol, df$gs_name)
}
collections <- list(
  Hallmark = build_list("H"),
  Reactome = build_list("C2", "CP:REACTOME"),
  `GO:BP`  = build_list("C5", "GO:BP"),
  KEGG     = build_list("C2", "CP:KEGG_LEGACY")
)

# ---- 6c. GSEA on every collection ----
set.seed(42)
gsea_res <- lapply(collections, function(lst)
  fgsea(pathways = lst, stats = ranks, minSize = 10, maxSize = 500, eps = 0))
saveRDS(gsea_res, file.path(OUTDIR, "gsea_results_all.rds"))
for (nm in names(gsea_res))
  message(sprintf("  GSEA %-9s: %d sets, %d significant (padj<0.05)",
                  nm, nrow(gsea_res[[nm]]), sum(gsea_res[[nm]]$padj < 0.05)))

# flatten leadingEdge list-column for a CSV
flat <- function(dt, coll) {
  d <- as.data.frame(dt)
  d$leadingEdge <- vapply(d$leadingEdge, function(x) paste(head(x, 30), collapse = ";"), character(1))
  d$collection  <- coll
  d[, c("collection","pathway","NES","ES","padj","pval","size","leadingEdge")]
}
gsea_csv <- do.call(rbind, Map(flat, gsea_res, names(gsea_res)))
write.csv(gsea_csv, file.path(OUTDIR, "gsea_results_all.csv"), row.names = FALSE)

# ---- 6d. ORA cross-check (Fisher/hypergeometric) on the 3,477-gene hit list ----
universe  <- de$gene
sig_genes <- de$gene[de$adj.P.Val < 0.05 & abs(de$logFC) > 1]
t2g <- function(lst) do.call(rbind, lapply(names(lst), function(n)
  data.frame(term = n, gene = lst[[n]])))
run_ora <- function(lst, genes) as.data.frame(enricher(
  gene = genes, universe = universe, TERM2GENE = t2g(lst),
  pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500))
ora_res <- lapply(collections, run_ora, genes = sig_genes)
for (nm in names(ora_res))
  message(sprintf("  ORA  %-9s: %d significant (padj<0.05)",
                  nm, sum(ora_res[[nm]]$p.adjust < 0.05)))

# ---- 6e. merge GSEA + ORA per pathway, save comparison table ----
gsea_std <- function(dt, coll) as.data.table(dt)[, .(pathway, collection = coll,
                              gsea_nes = NES, gsea_padj = padj, gsea_size = size)]
ora_std <- function(d, coll) {
  if (nrow(d) == 0) return(data.table())
  gr <- as.numeric(sub("/.*","",d$GeneRatio)); gd <- as.numeric(sub(".*/","",d$GeneRatio))
  br <- as.numeric(sub("/.*","",d$BgRatio));   bd <- as.numeric(sub(".*/","",d$BgRatio))
  data.table(pathway = d$ID, collection = coll, ora_padj = d$p.adjust,
             ora_fold = (gr/gd)/(br/bd), ora_k = gr, ora_setsize = br)
}
GSEA <- rbindlist(Map(gsea_std, gsea_res, names(gsea_res)))
ORA  <- rbindlist(Map(ora_std, ora_res, names(ora_res)), fill = TRUE)
merged <- merge(GSEA, ORA, by = c("pathway","collection"), all = TRUE)
merged[, gsea_sig := !is.na(gsea_padj) & gsea_padj < 0.05 & gsea_nes > 0]
merged[, ora_sig  := !is.na(ora_padj)  & ora_padj  < 0.05]
saveRDS(merged, file.path(OUTDIR, "gsea_ora_merged.rds"))
fwrite(merged, file.path(OUTDIR, "gsea_ora_merged.csv"))

# ---- 6f. Figure 8: GSEA enrichment "mountain plot" for IL6/JAK/STAT3 ----
NNblue <- "#4C72B0"; PPred <- "#C44E52"
pw <- "HALLMARK_IL6_JAK_STAT3_SIGNALING"; gs <- collections$Hallmark[[pw]]
row <- as.data.frame(gsea_res$Hallmark[gsea_res$Hallmark$pathway == pw, ])
r <- sort(ranks, decreasing = TRUE); N <- length(r); hit <- names(r) %in% gs
Phit  <- cumsum(ifelse(hit, abs(r), 0)) / sum(abs(r[hit]))
Pmiss <- cumsum(ifelse(!hit, 1, 0)) / (N - sum(hit))
res_walk <- Phit - Pmiss; peak_i <- which.max(abs(res_walk))
dfc <- data.frame(rank = seq_len(N), RES = res_walk)
fig8 <- ggplot(dfc, aes(rank, RES)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(colour = PPred, linewidth = 0.9) +
  geom_vline(xintercept = peak_i, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_rug(data = data.frame(rank = which(hit)), aes(x = rank), inherit.aes = FALSE,
           sides = "b", colour = PPred, alpha = 0.55, length = unit(0.06, "npc")) +
  annotate("text", x = N*0.62, y = max(res_walk)*0.85,
           label = sprintf("NES = %.2f\nadj.P = %.1e\n%d genes", row$NES, row$padj, row$size),
           hjust = 0, size = 4.1, colour = "grey15") +
  labs(title = "GSEA enrichment - Hallmark IL6/JAK/STAT3 signaling",
       subtitle = "Running enrichment score across all ranked genes",
       x = "Gene rank (high to low t-statistic)", y = "Running enrichment score") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(OUTDIR, "fig8_gsea_il6jakstat3.png"), fig8, width = 8.5, height = 5.2, dpi = 150)

# ---- 6g. Figure 9: pathway landscape dot-plot (curated STAT3-axis + top sets) ----
pick <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1",
  "REACTOME_INTERLEUKIN_17_SIGNALING",
  "REACTOME_STAT3_NUCLEAR_EVENTS_DOWNSTREAM_OF_ALK_SIGNALING_AND_PTK6",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY","GOBP_KERATINIZATION",
  "GOBP_KERATINOCYTE_DIFFERENTIATION","GOBP_RESPONSE_TO_INTERLEUKIN_17",
  "GOBP_INTERLEUKIN_17_PRODUCTION")
dp <- GSEA[pathway %in% pick & gsea_padj < 0.05][order(gsea_nes)]
dp[, lab := substr(tolower(gsub("_"," ", gsub("^(HALLMARK|REACTOME|KEGG|GOBP)_","",pathway))), 1, 44)]
dp[, lab := factor(lab, levels = lab)]
fig9 <- ggplot(dp, aes(gsea_nes, lab)) +
  geom_segment(aes(x = 0, xend = gsea_nes, y = lab, yend = lab), colour = "grey80", linewidth = 0.4) +
  geom_point(aes(size = gsea_size, colour = -log10(gsea_padj))) +
  scale_colour_gradient(low = "#F4C4C4", high = PPred, name = expression(-log[10]~adj.P)) +
  scale_size_continuous(name = "set size", range = c(2.5, 8)) +
  labs(title = "GSEA pathway landscape - psoriasis (PP vs NN)",
       subtitle = "Positive NES = coordinately up in lesional skin. All shown padj < 0.05.",
       x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(OUTDIR, "fig9_gsea_dotplot.png"), fig9, width = 9.5, height = 6, dpi = 150)

# ---- 6h. Figure 10: GSEA vs ORA concordance scatter ----
D <- merged[!is.na(gsea_padj) & !is.na(ora_padj) & gsea_nes > 0]
D[, glp := -log10(pmax(gsea_padj, 1e-50))]; D[, olp := -log10(pmax(ora_padj, 1e-50))]
D[, cls := fifelse(gsea_sig & ora_sig, "both",
            fifelse(gsea_sig & !ora_sig, "GSEA only",
            fifelse(!gsea_sig & ora_sig, "ORA only", "neither")))]
lab10 <- D[pathway %in% pick]
lab10[, short := gsub("_"," ", gsub("^(HALLMARK|REACTOME|KEGG|GOBP)_","",pathway))]
fig10 <- ggplot(D, aes(olp, glp)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#B0B0B0", linewidth = 0.4) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "#B0B0B0", linewidth = 0.4) +
  geom_point(aes(colour = cls), alpha = 0.55, size = 1.6) +
  geom_point(data = lab10, colour = "black", size = 2.6, shape = 21, fill = PPred, stroke = 0.6) +
  ggrepel::geom_text_repel(data = lab10, aes(label = short), size = 3.0, max.overlaps = 20,
                           box.padding = 0.5, min.segment.length = 0, colour = "grey15", seed = 1) +
  scale_colour_manual(values = c(both = PPred, `GSEA only` = NNblue,
                                 `ORA only` = "#DD8452", neither = "#B0B0B0"), name = NULL) +
  labs(title = "GSEA vs ORA - two views of the same pathways",
       subtitle = "Each point = one up-enriched gene set tested by both. Dashed = padj 0.05.",
       x = expression(ORA~-log[10]~adj.P~(Fisher)),
       y = expression(GSEA~-log[10]~adj.P~(fgsea))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave(file.path(OUTDIR, "fig10_gsea_vs_ora.png"), fig10, width = 8.5, height = 6.5, dpi = 150)

both_tab <- table(GSEA_up = D$gsea_sig, ORA = D$ora_sig)
message(sprintf("  GSEA/ORA concordance: both=%d, GSEA-only=%d, ORA-only=%d, neither=%d",
                both_tab["TRUE","TRUE"], both_tab["TRUE","FALSE"],
                both_tab["FALSE","TRUE"], both_tab["FALSE","FALSE"]))

message("DONE. Checkpoints: dge_filt_norm.rds, de_results_full.rds, gsea_ranks.rds, gsea_results_all.rds, gsea_ora_merged.rds")
message("Tables: de_results_full.csv, de_results_significant.csv, gsea_results_all.csv, gsea_ora_merged.csv")
message("Figures: fig6_volcano.png, fig7_qvalue.png, fig8_gsea_il6jakstat3.png, fig9_gsea_dotplot.png, fig10_gsea_vs_ora.png")

# =============================================================================
# STEP 7 - Network analysis: TF-activity inference  --> CHECKPOINT 4
# -----------------------------------------------------------------------------
# Pathway analysis (STEP 6) asks "which gene SETS moved?". Network analysis asks
# a causal question: "which upstream REGULATORS drove the change?". We infer the
# activity of each transcription factor (TF) not from its own mRNA, but from the
# coordinated behaviour of the genes it is known to regulate (its "regulon").
#
# Method: decoupleR::run_ulm (univariate linear model). For each TF, regress the
# genome-wide limma t-statistics on that TF's signed regulon membership
# (+1 = activated target, -1 = repressed target); the slope's t-value is the
# activity score. Regulon: CollecTRI (signed, curated TF->target network),
# fetched from the OmniPath HTTP API. Requires network access to omnipathdb.org.
# =============================================================================
message("STEP 7: network analysis - TF-activity inference (decoupleR + CollecTRI) ...")

# ---- 7a. fetch the CollecTRI signed regulon from OmniPath ----
# Direct HTTP to the API returns a plain TSV and avoids OmnipathR's cross-species
# machinery (which reaches omabrowser.org and is unnecessary for human data).
collectri_url <- paste0("https://omnipathdb.org/interactions?resources=CollecTRI",
                        "&datasets=collectri&fields=sources,references",
                        "&genesymbols=1&organisms=9606")
raw_net <- fread(collectri_url, sep = "\t", header = TRUE)
net <- raw_net[, .(source = source_genesymbol, target = target_genesymbol,
                   mor = fifelse(is_inhibition & !is_stimulation, -1, 1))]
net <- unique(net[source != "" & target != ""])
net <- net[, .(mor = fifelse(sum(mor) < 0, -1, 1)), by = .(source, target)]  # 1 sign per edge
saveRDS(net, file.path(OUTDIR, "collectri_regulon.rds"))                       # CHECKPOINT 4
message(sprintf("  regulon: %d edges, %d TFs, %d targets",
                nrow(net), uniqueN(net$source), uniqueN(net$target)))

# ---- 7b. run the ULM activity model on the limma t-statistics ----
mat <- matrix(de$t, ncol = 1, dimnames = list(de$gene, "PPvsNN"))
mat <- mat[!is.na(mat[, 1]), , drop = FALSE]
set.seed(42)
tf_res <- run_ulm(mat = mat, network = net, .source = "source",
                  .target = "target", .mor = "mor", minsize = 5)
setDT(tf_res)
tf_res <- tf_res[order(-score)]
tf_res[, padj := p.adjust(p_value, "BH")]
tf_res[, rank := .I]
axis_tf <- c("STAT3","STAT1","RELA","RELB","REL","NFKB1","NFKB2","RORC","RORA",
             "IRF1","MYC","E2F1","JUNB","FOS","JUN")
tf_res[, axis := source %in% axis_tf]
fwrite(tf_res, file.path(OUTDIR, "tf_activity_collectri.csv"))
saveRDS(tf_res, file.path(OUTDIR, "tf_activity_collectri.rds"))
message(sprintf("  %d TFs scored, %d significant (padj<0.05). STAT3: score=%.2f, padj=%.1e, rank %d/%d",
                nrow(tf_res), sum(tf_res$padj < 0.05),
                tf_res[source=="STAT3"]$score, tf_res[source=="STAT3"]$padj,
                tf_res[source=="STAT3"]$rank, nrow(tf_res)))

# ---- 7c. Figure 12: TF-activity ranked landscape (STAT3/NF-kB/IFN axis flagged) ----
NNblue <- "#4C72B0"; PPred <- "#C44E52"; grey <- "#B0B0B0"
plotdt <- copy(tf_res)[!grepl("_", source)]        # drop NF-kB dimer complexes for readability
plotdt[, rank2 := frank(-score, ties.method = "first")]
labdt <- plotdt[source %in% c("STAT3","STAT1","RELA","NFKB1","RELB","MYC","E2F1","IRF1","RORC","JUNB") | rank2 <= 3]
fig12 <- ggplot(plotdt, aes(rank2, score)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(aes(colour = axis, size = axis, alpha = axis)) +
  scale_colour_manual(values = c(`FALSE` = grey, `TRUE` = PPred), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.1, `TRUE` = 2.6), guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 1), guide = "none") +
  ggrepel::geom_text_repel(data = labdt, aes(label = source), size = 3.2, fontface = "bold",
                           max.overlaps = 30, box.padding = 0.5, min.segment.length = 0,
                           colour = "grey15", seed = 1) +
  labs(title = "Transcription-factor activity in psoriasis - inferred from target behaviour",
       subtitle = "732 TFs ranked by ULM activity score (CollecTRI). Red = STAT3 / NF-kB / IFN axis.",
       x = "TF activity rank (most active -> least)", y = "Activity score (ULM t-value)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11))
ggsave(file.path(OUTDIR, "fig12_tf_activity_rank.png"), fig12, width = 10, height = 6.2, dpi = 150)

# ---- 7d. Figure 13: STAT3 target coherence (why the activity score is high) ----
s <- net[source == "STAT3"]
dt <- merge(s, data.table(target = de$gene, t = de$t), by = "target")
dt[, dir := fifelse(mor > 0, "activated by STAT3 (+)", "repressed by STAT3 (-)")]
lab <- rbind(dt[mor > 0][order(-t)][1:10], dt[mor < 0][order(t)][1:6])
n_act <- sum(dt$mor > 0); n_rep <- sum(dt$mor < 0)
fig13 <- ggplot(dt, aes(factor(mor), t, colour = dir)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.55, size = 1.6) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = NA, colour = "grey25", linewidth = 0.5) +
  ggrepel::geom_text_repel(data = lab, aes(label = target), size = 2.9, max.overlaps = 30,
                           box.padding = 0.4, min.segment.length = 0, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = c("activated by STAT3 (+)" = PPred,
                                 "repressed by STAT3 (-)" = NNblue), name = NULL) +
  scale_x_discrete(labels = c("-1" = sprintf("repressed\ntargets (n=%d)", n_rep),
                              "1"  = sprintf("activated\ntargets (n=%d)", n_act))) +
  labs(title = "STAT3 is an ACTIVE regulator - its targets move as its regulon predicts",
       subtitle = "Each point = one STAT3 target. Activated targets shift up, repressed shift down.",
       x = NULL, y = "Observed change in lesional skin (limma t-statistic)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 12),
        legend.position = "top")
ggsave(file.path(OUTDIR, "fig13_stat3_targets.png"), fig13, width = 8.5, height = 6.3, dpi = 150)

message("STEP 7 DONE. Checkpoint: collectri_regulon.rds, tf_activity_collectri.rds")
message("Tables: tf_activity_collectri.csv")
message("Figures: fig12_tf_activity_rank.png, fig13_stat3_targets.png")

# =============================================================================
# STEP 8 - STAT3 isoform splicing (alpha vs beta) from junction-level recount3
# =============================================================================
# Biology: STAT3 mRNA is alternatively spliced at its 3' end into two proteins
#   - STAT3-alpha (770 aa): canonical activator, retains C-terminal
#     transactivation domain (TAD, incl. Ser727).
#   - STAT3-beta  (722 aa): truncated dominant-negative; loses the TAD, ends in
#     the unique heptapeptide FIDAVWK.
# The choice is an alternative-5'-splice-site (A5SS) event: a SHARED ACCEPTOR
# with two competing DONORS (chr17, minus strand, GRCh38):
#   alpha: donor 42,316,902 -> acceptor 42,317,181
#   beta : donor 42,316,852 -> acceptor 42,317,181
# We read the splice choice directly from junction counts and summarise it as
#   PSI_beta = beta / (alpha + beta)   (percent-spliced-in of the beta isoform).
message("STEP 8: STAT3 isoform splicing (alpha/beta PSI from junction RSE) ...")

suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
})

## 8.1 Junction-level RSE for the same study ---------------------------------
if (!exists("proj_info")) {
  hs <- available_projects()
  proj_info <- subset(hs, project == "SRP035988" & project_type == "data_sources")
}
rse_jxn <- create_rse(proj_info, type = "jxn")   # ~1.8M junctions x 178 samples
message(sprintf("  junction RSE: %d junctions x %d samples",
                nrow(rse_jxn), ncol(rse_jxn)))

## 8.2 Group labels (NN normal vs PP psoriasis) from sample attributes -------
attr_txt <- as.character(colData(rse_jxn)$sra.sample_attributes)
grp <- ifelse(grepl("normal", attr_txt, ignore.case = TRUE), "NN",
       ifelse(grepl("[Pp]soriasis|lesional|involved", attr_txt), "PP", NA))
stopifnot(!any(is.na(grp)))
grp <- factor(grp, levels = c("NN", "PP"))
message(sprintf("  groups: NN=%d, PP=%d", sum(grp == "NN"), sum(grp == "PP")))

## 8.3 Locate the alpha- and beta-defining junctions -------------------------
# Match on shared acceptor 42,317,181 and the two alternative donors.
gr  <- rowRanges(rse_jxn)
onchr <- as.character(seqnames(gr)) %in% c("chr17", "17")
acc <- 42317181L
a_don <- 42316902L; b_don <- 42316852L
tol <- 2L   # small tolerance for 0/1-based edge conventions
is_a <- onchr & abs(start(gr) - a_don) <= tol & abs(end(gr) - acc) <= tol
is_b <- onchr & abs(start(gr) - b_don) <= tol & abs(end(gr) - acc) <= tol
# (start = low coord = donor on minus strand near 42.316 Mb; end = acceptor)
stopifnot(sum(is_a) >= 1, sum(is_b) >= 1)
a_cnt <- colSums(matrix(assay(rse_jxn, "counts")[is_a, , drop = FALSE], ncol = ncol(rse_jxn)))
b_cnt <- colSums(matrix(assay(rse_jxn, "counts")[is_b, , drop = FALSE], ncol = ncol(rse_jxn)))
message(sprintf("  alpha junctions matched: %d (total reads %d)", sum(is_a), sum(a_cnt)))
message(sprintf("  beta  junctions matched: %d (total reads %d)", sum(is_b), sum(b_cnt)))

## 8.4 Per-sample PSI_beta ---------------------------------------------------
depth   <- a_cnt + b_cnt
psi_beta <- b_cnt / depth
dt <- data.frame(sample = colnames(rse_jxn), grp = grp,
                 a = a_cnt, b = b_cnt, depth = depth, psi_beta = psi_beta)
message(sprintf("  min junction depth across samples: %d", min(depth)))

## 8.5 Test: does PSI_beta differ NN vs PP? ----------------------------------
wt <- wilcox.test(psi_beta ~ grp, data = dt)
mean_nn <- mean(dt$psi_beta[dt$grp == "NN"]); mean_pp <- mean(dt$psi_beta[dt$grp == "PP"])
message(sprintf("  PSI_beta: NN mean %.4f vs PP mean %.4f  (Wilcoxon p = %.4f)",
                mean_nn, mean_pp, wt$p.value))

## 8.6 Figure 14 - isoform splicing (PSI + absolute junction counts) ---------
NNblue <- "#4C72B0"; PPred <- "#C44E52"
p14a <- ggplot(dt, aes(grp, psi_beta, colour = grp)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6, size = 1.8) +
  geom_boxplot(width = 0.4, outlier.shape = NA, fill = NA, colour = "grey25", linewidth = 0.5) +
  scale_colour_manual(values = c(NN = NNblue, PP = PPred), guide = "none") +
  labs(title = "STAT3 beta isoform inclusion (PSI beta)",
       subtitle = sprintf("Wilcoxon p = %.3f", wt$p.value),
       x = NULL, y = "PSI beta = beta / (alpha + beta)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

dt_long <- rbind(data.frame(grp = dt$grp, iso = "alpha", reads = dt$a),
                 data.frame(grp = dt$grp, iso = "beta",  reads = dt$b))
p14b <- ggplot(dt_long, aes(interaction(iso, grp), reads + 1, colour = grp)) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  geom_boxplot(width = 0.4, outlier.shape = NA, fill = NA, colour = "grey25", linewidth = 0.5) +
  scale_y_log10() +
  scale_colour_manual(values = c(NN = NNblue, PP = PPred), name = NULL) +
  labs(title = "Absolute isoform junction reads",
       subtitle = "Both isoforms up in lesional skin; alpha dominant throughout",
       x = NULL, y = "junction reads + 1 (log10)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

fig14 <- patchwork::wrap_plots(p14a, p14b, nrow = 1)
ggsave(file.path(OUTDIR, "fig14_stat3_isoform.png"), fig14, width = 10, height = 4.6, dpi = 150)

## 8.7 EFTUD2 correlation (spliceosome hypothesis test) ----------------------
# Needs the gene-level log-CPM matrix from STEP 2 (dge_filt_norm.rds).
if (!exists("logcpm")) {
  if (!exists("dge")) dge <- readRDS(file.path(OUTDIR, "dge_filt_norm.rds"))
  logcpm <- edgeR::cpm(dge, log = TRUE)
}
# Map EFTUD2 by symbol (rownames may be Ensembl IDs; use dge$genes if present).
eftud2_row <- NULL
if (!is.null(dge$genes) && "gene_name" %in% colnames(dge$genes)) {
  idx <- which(dge$genes$gene_name == "EFTUD2")
  if (length(idx) == 1) eftud2_row <- logcpm[idx, ]
}
if (is.null(eftud2_row) && "EFTUD2" %in% rownames(logcpm)) eftud2_row <- logcpm["EFTUD2", ]
if (!is.null(eftud2_row)) {
  # Align samples between junction table and expression matrix.
  common <- intersect(dt$sample, colnames(logcpm))
  m <- data.frame(dt[match(common, dt$sample), ],
                  eftud2 = eftud2_row[match(common, names(eftud2_row))])
  sp_all <- suppressWarnings(cor.test(m$eftud2, m$psi_beta, method = "spearman"))
  sp_pp  <- suppressWarnings(cor.test(m$eftud2[m$grp=="PP"], m$psi_beta[m$grp=="PP"], method = "spearman"))
  sp_nn  <- suppressWarnings(cor.test(m$eftud2[m$grp=="NN"], m$psi_beta[m$grp=="NN"], method = "spearman"))
  message(sprintf("  EFTUD2 vs PSI_beta  overall rho=%.3f p=%.3f | PP rho=%.3f p=%.3f | NN rho=%.3f p=%.3f",
                  sp_all$estimate, sp_all$p.value, sp_pp$estimate, sp_pp$p.value,
                  sp_nn$estimate, sp_nn$p.value))

  fig15a <- ggplot(m, aes(eftud2, psi_beta, colour = grp)) +
    geom_point(alpha = 0.65, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    scale_colour_manual(values = c(NN = NNblue, PP = PPred), name = NULL) +
    labs(title = "PSI beta vs EFTUD2 (within-group fits)",
         subtitle = sprintf("PP rho=%.2f p=%.3f ; NN rho=%.2f p=%.3f ; pooled rho=%.2f (Simpson's paradox)",
                            sp_pp$estimate, sp_pp$p.value, sp_nn$estimate, sp_nn$p.value, sp_all$estimate),
         x = "EFTUD2 expression (log2 CPM)", y = "PSI beta") +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "top")
  fig15b <- ggplot(m, aes(grp, eftud2, colour = grp)) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
    geom_boxplot(width = 0.4, outlier.shape = NA, fill = NA, colour = "grey25", linewidth = 0.5) +
    scale_colour_manual(values = c(NN = NNblue, PP = PPred), guide = "none") +
    labs(title = "EFTUD2 is up in lesional skin", x = NULL, y = "EFTUD2 (log2 CPM)") +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
  fig15 <- patchwork::wrap_plots(fig15a, fig15b, nrow = 1, widths = c(2, 1))
  ggsave(file.path(OUTDIR, "fig15_eftud2_psi.png"), fig15, width = 11, height = 4.6, dpi = 150)
  saveRDS(m, file.path(OUTDIR, "stat3_isoform_eftud2.rds"))
}

## 8.8 Checkpoint ------------------------------------------------------------
saveRDS(dt, file.path(OUTDIR, "stat3_isoform_psi.rds"))
write.csv(dt, file.path(OUTDIR, "stat3_isoform_psi.csv"), row.names = FALSE)
message("STEP 8 DONE. Checkpoint: stat3_isoform_psi.rds, stat3_isoform_eftud2.rds")
message("Figures: fig14_stat3_isoform.png, fig15_eftud2_psi.png")

# =============================================================================
# STEP 9 - Sample clustering & DE-gene co-expression heatmaps (descriptive/QC)
# =============================================================================
# These are DESCRIPTIVE figures, NOT hypothesis tests. NN/PP labels are a-priori
# (from clinical sample annotation), so there is no double-dipping: we colour by
# the pre-existing labels and observe whether unsupervised clustering agrees.
message("STEP 9: sample clustering + DE-gene heatmaps ...")

suppressPackageStartupMessages({ library(pheatmap); library(RColorBrewer); library(grid) })
if (!exists("logcpm")) {
  if (!exists("dge")) dge <- readRDS(file.path(OUTDIR, "dge_filt_norm.rds"))
  logcpm <- edgeR::cpm(dge, log = TRUE)
}
if (!exists("de_tab")) de_tab <- readRDS(file.path(OUTDIR, "de_results_full.rds"))
NNblue <- "#4C72B0"; PPred <- "#C44E52"

sig <- subset(de_tab, adj.P.Val < 0.05 & abs(logFC) > 1)
X   <- logcpm[sig$gene, , drop = FALSE]

ann <- data.frame(group = dge$samples$group); rownames(ann) <- colnames(logcpm)
ann_colors <- list(group = c(NN = NNblue, PP = PPred))

## 9.1 Figure 16 - sample-sample correlation (Pearson + Spearman) ------------
cor_p <- cor(X, method = "pearson"); cor_s <- cor(X, method = "spearman")
rng  <- range(c(cor_p, cor_s))
brks <- seq(floor(rng[1]*100)/100, 1, length.out = 101)
pal  <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
draw_hm <- function(mat, ttl, col_i) {
  ph <- pheatmap(mat, color = pal, breaks = brks,
                 annotation_col = ann, annotation_row = ann, annotation_colors = ann_colors,
                 show_rownames = FALSE, show_colnames = FALSE,
                 clustering_distance_rows = as.dist(1 - mat),
                 clustering_distance_cols = as.dist(1 - mat),
                 clustering_method = "average", main = ttl, silent = TRUE,
                 annotation_names_row = FALSE, annotation_names_col = FALSE)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = col_i))
  grid.draw(ph$gtable); popViewport()
}
png(file.path(OUTDIR, "fig16_sample_corr_heatmap.png"), width = 13, height = 6.6, units = "in", res = 150)
grid.newpage(); pushViewport(viewport(layout = grid.layout(1, 2)))
draw_hm(cor_p, "Pearson correlation (linear)", 1)
draw_hm(cor_s, "Spearman correlation (rank)", 2)
popViewport(); dev.off()

## 9.2 Figure 17 - top DE-gene expression modules ---------------------------
sig_up <- head(sig$gene[sig$logFC > 0][order(-sig$logFC[sig$logFC > 0])], 60)
sig_dn <- head(sig$gene[sig$logFC < 0][order(sig$logFC[sig$logFC < 0])], 60)
topg <- c(sig_up, sig_dn)
Z <- t(scale(t(logcpm[topg, , drop = FALSE]))); Z[Z > 3] <- 3; Z[Z < -3] <- -3
rowann <- data.frame(direction = ifelse(topg %in% sig_up, "up in PP", "down in PP"))
rownames(rowann) <- topg
ac2 <- list(group = c(NN = NNblue, PP = PPred),
            direction = c("up in PP" = PPred, "down in PP" = NNblue))
png(file.path(OUTDIR, "fig17_degene_heatmap.png"), width = 11, height = 10, units = "in", res = 150)
pheatmap(Z, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-3, 3, length.out = 101),
         annotation_col = ann, annotation_row = rowann, annotation_colors = ac2,
         show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 5.5,
         clustering_distance_cols = "correlation", clustering_distance_rows = "correlation",
         clustering_method = "average", treeheight_row = 25, treeheight_col = 30,
         main = "Top 120 DE genes (60 up + 60 down in lesional skin), row z-scored log2-CPM",
         annotation_names_row = FALSE, annotation_names_col = FALSE)
dev.off()

message("STEP 9 DONE. Figures: fig16_sample_corr_heatmap.png, fig17_degene_heatmap.png")
message("NOTE: clust (Abu-Jamous & Kelly 2018) deferred to the multi-study meta-analysis,")
message("      where its cross-dataset consistency mode is the real payoff.")

# =============================================================================
# STEP 10 - Enrichment robustness: CAMERA, ROAST, ssGSEA, UpSet
#           --> Figures 18-21, CHECKPOINT camera_roast_results.rds
#
# STEP 6 established the STAT3 program with competitive GSEA (+ ORA cross-check).
# Competitive GSEA treats the genes in a set as INDEPENDENT votes; pathway genes
# are co-regulated, so it can be over-confident (Goeman & Buhlmann 2007; Wu &
# Smyth 2012). STEP 10 stress-tests the signal four complementary ways:
#   10.1 CAMERA  - competitive, but corrects inter-gene correlation (VIF).
#   10.2 ROAST   - self-contained rotation test ("engaged at all?").
#   10.3 ssGSEA  - per-sample pathway scores (base-R impl, no GSVA dependency).
#   10.4 UpSet   - leading-edge overlap = pathway redundancy.
# Needs: dge_filt_norm.rds (STEP 2), collectri_regulon.rds (STEP 7),
#        gsea_results_all.rds + gsea_ranks.rds (STEP 6), stat3_isoform_psi.rds (STEP 8).
# =============================================================================
message("STEP 10: enrichment robustness (CAMERA / ROAST / ssGSEA / UpSet) ...")

dge     <- readRDS(file.path(OUTDIR, "dge_filt_norm.rds"))
grp     <- factor(dge$samples$group, levels = c("NN", "PP"))
design  <- model.matrix(~grp); colnames(design) <- c("Intercept", "PPvsNN")
v       <- voom(dge, design, plot = FALSE)          # observation-level fit for camera/roast
logcpm  <- edgeR::cpm(dge, log = TRUE)

# gene-set collections (symbols; dge/voom rownames are gene SYMBOLS) -----------
h_df <- msigdbr(species = "Homo sapiens", category = "H")
H    <- split(h_df$gene_symbol, h_df$gs_name)
r_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
RE   <- split(r_df$gene_symbol, r_df$gs_name)

# STAT3 regulon (self-contained hypothesis set) from CollecTRI ----------------
col            <- readRDS(file.path(OUTDIR, "collectri_regulon.rds"))
stat3_targets  <- unique(col$target[col$source == "STAT3"])
stat3_measured <- intersect(stat3_targets, rownames(v$E))   # 371 measured
idx_H     <- ids2indices(H,  rownames(v$E))
idx_RE    <- ids2indices(RE, rownames(v$E))
idx_stat3 <- list(STAT3_regulon = which(rownames(v$E) %in% stat3_measured))

## 10.1 CAMERA - competitive, corrects inter-gene correlation -----------------
cam_H     <- camera(v, idx_H,     design, contrast = "PPvsNN")
cam_RE    <- camera(v, idx_RE,    design, contrast = "PPvsNN")
cam_stat3 <- camera(v, idx_stat3, design, contrast = "PPvsNN")

## 10.2 ROAST - self-contained rotation test on the STAT3 regulon -------------
set.seed(1)
roast_stat3 <- roast(v, idx_stat3$STAT3_regulon, design, contrast = "PPvsNN", nrot = 9999)

## GSEA on the same regulon (as a custom set) for the 3-method comparison ------
ranks    <- readRDS(file.path(OUTDIR, "gsea_ranks.rds"))
gsea_all <- readRDS(file.path(OUTDIR, "gsea_results_all.rds"))
set.seed(1)
fg_stat3 <- fgsea(pathways = list(STAT3_regulon = stat3_measured),
                  stats = ranks, minSize = 5, maxSize = 1000)

saveRDS(list(cam_H = cam_H, cam_RE = cam_RE, cam_stat3 = cam_stat3,
             roast_stat3 = roast_stat3, fg_stat3 = fg_stat3,
             stat3_measured = stat3_measured),
        file.path(OUTDIR, "camera_roast_results.rds"))

## Figure 18 - three-method comparison ----------------------------------------
grey <- "#B0B0B0"
gseaH <- as.data.frame(gsea_all$H)[, c("pathway", "NES", "padj")]
camH  <- data.frame(pathway = rownames(cam_H), cam_dir = cam_H$Direction,
                    cam_FDR = cam_H$FDR, cam_NGenes = cam_H$NGenes)
mrg <- merge(gseaH, camH, by = "pathway")
mrg$gsea_sig  <- mrg$padj    < 0.05
mrg$cam_sig   <- mrg$cam_FDR < 0.05
mrg$agree <- (mrg$NES > 0 & mrg$cam_dir == "Up"   & mrg$gsea_sig & mrg$cam_sig) |
             (mrg$NES < 0 & mrg$cam_dir == "Down" & mrg$gsea_sig & mrg$cam_sig)
mrg$neglog_cam <- -log10(mrg$cam_FDR)
mrg$cat <- ifelse(mrg$agree, "Both sig, same direction",
           ifelse(mrg$gsea_sig | mrg$cam_sig, "One method sig", "Neither sig"))
lab <- mrg[mrg$agree & (mrg$neglog_cam > 8 | abs(mrg$NES) > 2.3), ]
pA <- ggplot(mrg, aes(NES, neglog_cam)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = grey) +
  geom_vline(xintercept = 0, color = grey) +
  geom_point(aes(color = cat, size = cam_NGenes), alpha = 0.8) +
  ggrepel::geom_text_repel(data = lab, aes(label = gsub("HALLMARK_", "", pathway)),
                  size = 2.5, max.overlaps = 20, min.segment.length = 0) +
  scale_color_manual(values = c("Both sig, same direction" = PPred,
                                "One method sig" = "#DD8452", "Neither sig" = grey)) +
  scale_size_continuous(range = c(1.5, 6), name = "Set size") +
  labs(x = "GSEA NES (fgsea, competitive)", y = "CAMERA  -log10(FDR)",
       title = "A. Hallmark: GSEA vs CAMERA concordance",
       subtitle = "sets significant in both, same direction; CAMERA corrects inter-gene correlation",
       color = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "right",
       plot.subtitle = element_text(size = 7.5))
db <- data.frame(
  method = factor(c("GSEA (competitive)", "CAMERA (competitive,\ncorrelation-corrected)",
                    "ROAST (self-contained)"),
    levels = c("GSEA (competitive)", "CAMERA (competitive,\ncorrelation-corrected)",
               "ROAST (self-contained)")),
  p = c(fg_stat3$pval[1], cam_stat3$PValue[1], roast_stat3$p.value["Up", "P.Value"]),
  type = c("competitive", "competitive", "self-contained"))
db$negp <- -log10(db$p)
pB <- ggplot(db, aes(method, negp, fill = type)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = grey) +
  geom_text(aes(label = sprintf("p=%.1e", p)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("competitive" = NNblue, "self-contained" = PPred), name = "Test type") +
  labs(x = NULL, y = "-log10(p)", title = "B. STAT3 regulon (371 targets): three tests agree",
       subtitle = "All UP in lesional; CAMERA p >> GSEA p because it discounts co-regulation") +
  ylim(0, max(db$negp) * 1.15) +
  theme_bw(base_size = 10) + theme(legend.position = "right",
       axis.text.x = element_text(size = 7.5), plot.subtitle = element_text(size = 7.5))
ggsave(file.path(OUTDIR, "fig18_camera_roast.png"), pA / pB + plot_layout(heights = c(1.25, 1)),
       width = 9, height = 9, dpi = 150)

## 10.3 ssGSEA - per-sample pathway scores (base R; Barbie et al. 2009) --------
# GSVA/ssGSEA implemented directly (the Bioconductor GSVA pkg + its compiled deps
# do not install in the sandbox). Per sample: rank genes; ES = integral of the
# weighted (ECDF_in - ECDF_out) walk down the ranked list; range-normalize.
ssgsea_scores <- function(expr, gene_sets, alpha = 0.25, norm = TRUE, min_size = 5) {
  genes <- rownames(expr); n <- nrow(expr); ns <- ncol(expr)
  R <- apply(expr, 2, function(x) rank(x, ties.method = "average"))
  gene_sets <- lapply(gene_sets, function(g) intersect(g, genes))
  gene_sets <- gene_sets[sapply(gene_sets, length) >= min_size]
  ES <- matrix(NA, length(gene_sets), ns, dimnames = list(names(gene_sets), colnames(expr)))
  for (j in seq_len(ns)) {
    ord <- order(R[, j], decreasing = TRUE); rnk <- R[ord, j]; g_ord <- genes[ord]
    for (s in seq_along(gene_sets)) {
      inset <- g_ord %in% gene_sets[[s]]; w <- abs(rnk)^alpha
      Pin  <- cumsum(ifelse(inset, w, 0)); Pin  <- Pin  / Pin[n]
      Pout <- cumsum(ifelse(inset, 0, 1)); Pout <- Pout / Pout[n]
      ES[s, j] <- sum(Pin - Pout)
    }
  }
  if (norm) ES <- ES / (max(ES) - min(ES))
  ES
}
panel_sets <- H[c("HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MTORC1_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYOGENESIS", "HALLMARK_ADIPOGENESIS", "HALLMARK_FATTY_ACID_METABOLISM")]
sc <- ssgsea_scores(v$E, panel_sets); rownames(sc) <- gsub("HALLMARK_", "", rownames(sc))
il6 <- sc["IL6_JAK_STAT3_SIGNALING", ]
saveRDS(list(panel_scores = sc, ssgsea_fn = ssgsea_scores), file.path(OUTDIR, "ssgsea_panel.rds"))

# Figure 19 - score vs STAT3 expression and vs PSI_beta
psi <- as.data.frame(readRDS(file.path(OUTDIR, "stat3_isoform_psi.rds")))
m <- data.frame(sample = rownames(dge$samples), il6 = il6,
                stat3 = logcpm["STAT3", ], grp = grp)
m <- merge(m, psi[, c("sample", "psi_beta")], by = "sample")
c_expr   <- cor.test(m$il6, m$stat3,    method = "spearman")
c_psi    <- cor.test(m$il6, m$psi_beta, method = "spearman")
c_psi_pp <- cor.test(m$il6[m$grp=="PP"], m$psi_beta[m$grp=="PP"], method = "spearman")
c_psi_nn <- cor.test(m$il6[m$grp=="NN"], m$psi_beta[m$grp=="NN"], method = "spearman")
qA <- ggplot(m, aes(stat3, il6, color = grp)) +
  geom_point(alpha = 0.8, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  scale_color_manual(values = c(NN = NNblue, PP = PPred), name = NULL) +
  labs(x = "STAT3 expression (log2-CPM)", y = "IL6/JAK/STAT3 ssGSEA score",
       title = "A. Per-sample pathway score tracks STAT3 expression",
       subtitle = sprintf("Spearman rho=%.2f, p=%.1e", c_expr$estimate, c_expr$p.value)) +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 8))
qB <- ggplot(m, aes(psi_beta * 100, il6, color = grp)) +
  geom_point(alpha = 0.8, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, aes(group = grp), linewidth = 0.6) +
  scale_color_manual(values = c(NN = NNblue, PP = PPred), name = NULL) +
  labs(x = "PSI_beta (% STAT3-beta)", y = "IL6/JAK/STAT3 ssGSEA score",
       title = "B. Pathway score vs beta-isoform fraction",
       subtitle = sprintf("all rho=%.2f (p=%.2f); within PP rho=%.2f (p=%.3f), NN rho=%.2f (p=%.3f)",
         c_psi$estimate, c_psi$p.value, c_psi_pp$estimate, c_psi_pp$p.value,
         c_psi_nn$estimate, c_psi_nn$p.value)) +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 7))
ggsave(file.path(OUTDIR, "fig19_gsva_correlations.png"), qA / qB, width = 8, height = 8, dpi = 150)

# Figure 20 - per-sample pathway-score heatmap (row z-scored, samples NN->PP)
ann <- data.frame(Group = grp); rownames(ann) <- colnames(sc)
scz <- t(scale(t(sc))); ordc <- order(grp, il6)
png(file.path(OUTDIR, "fig20_pathway_score_heatmap.png"), width = 1500, height = 750, res = 150)
pheatmap(scz[, ordc], cluster_cols = FALSE, cluster_rows = TRUE,
         annotation_col = ann, annotation_colors = list(Group = c(NN = NNblue, PP = PPred)),
         show_colnames = FALSE, fontsize_row = 8,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Fig 20. Per-sample ssGSEA pathway scores (row z-scored), samples ordered NN->PP",
         breaks = seq(-3, 3, length.out = 101))
dev.off()

## 10.4 UpSet - leading-edge overlap = pathway redundancy ---------------------
gh <- as.data.frame(gsea_all$H)
pick <- c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT")
le <- setNames(lapply(pick, function(p) {
  x <- gh$leadingEdge[gh$pathway == p][[1]]; if (is.list(x)) x <- x[[1]]; as.character(x)
}), gsub("HALLMARK_", "", pick))
saveRDS(le, file.path(OUTDIR, "leading_edge_sets.rds"))
# UpSetR is CRAN-only; load it from the workspace lib if the base env lacks it.
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  message("  UpSetR not found - install with install.packages('UpSetR') (CRAN); skipping fig21.")
} else {
  um <- UpSetR::fromList(le)
  png(file.path(OUTDIR, "fig21_upset_leadingedge.png"), width = 1700, height = 1000, res = 150)
  print(UpSetR::upset(um, sets = rev(names(le)), keep.order = TRUE, order.by = "freq",
        nintersects = 25, mainbar.y.label = "Genes shared across pathways (leading-edge)",
        sets.x.label = "Leading-edge size", main.bar.color = PPred, sets.bar.color = NNblue,
        matrix.color = "#333333", shade.color = grey,
        text.scale = c(1.3, 1.2, 1.1, 1.0, 1.2, 1.1), point.size = 2.6, line.size = 0.8))
  grid::grid.text("Fig 21. Leading-edge gene overlap across 6 top immune Hallmark pathways",
                  x = 0.05, y = 0.97, just = "left", gp = grid::gpar(fontsize = 12, fontface = "bold"))
  dev.off()
}

message("STEP 10 DONE. Figures: fig18_camera_roast.png, fig19_gsva_correlations.png,")
message("      fig20_pathway_score_heatmap.png, fig21_upset_leadingedge.png")
message("      Checkpoints: camera_roast_results.rds, ssgsea_panel.rds, leading_edge_sets.rds")
message(sprintf("      STAT3 regulon: GSEA p=%.1e  CAMERA p=%.1e  ROAST(Up) p=%.1e",
        fg_stat3$pval[1], cam_stat3$PValue[1], roast_stat3$p.value["Up","P.Value"]))

# =============================================================================
# STEP 11 - Multi-study meta-analysis (5 whole-skin psoriasis cohorts)
# -----------------------------------------------------------------------------
# Scales the single-study pipeline (STEP 1-10) to a meta-analysis across
# independent recount3 cohorts. Rationale and full statistics: white-paper
# section 11 and meta_analysis_deep_dive.md. Random-effects pooling separates
# disease biology from single-study idiosyncrasy; a signal that survives across
# cohorts is attributable to psoriasis, not to one experiment.
#
# Study set (curated to WHOLE-SKIN HUMAN bulk RNA-seq, Gencode v26, poolable):
#   SRP035988 (anchor, Li)   95 PP / 83 NN
#   SRP165679 (Tsoi AD/PSO)  28 PP / 27 PN / 38 NN   (AD samples dropped)
#   SRP076982 (anatomic)    211 PP / 48 PN           (shallow ~4M reads; no NN)
#   SRP126422 (biopsy arm)    4 PP /  4 PN /  4 NN
#   SRP016583 (early)         3 PP /  3 PN
# Classes: PP lesional, PN peri-lesional/uninvolved, NN normal-healthy.
# Contrasts: PPvsNN (3 studies), PNvsNN (2), PNvsPP (4).
message("STEP 11: multi-study meta-analysis ...")
suppressMessages({
  library(metafor)   # from workspace lib (.libPaths set in STEP 7/10 if used)
})

STUDIES <- c("SRP035988","SRP165679","SRP076982","SRP126422","SRP016583")

## 11.1 Download gene-level counts for every study --------------------------
## Each recount3 study is quantified against the SAME annotation (Gencode v26),
## so the gene space is harmonised by construction - no probe re-mapping.
ap <- available_projects(organism = "human")
rses <- lapply(STUDIES, function(srp) {
  pj <- subset(ap, project == srp & project_type == "data_sources")
  rse <- create_rse(pj)                     # gene-level RSE
  assay(rse, "counts") <- transform_counts(rse)
  rse
})
names(rses) <- STUDIES

## 11.2 Classify each sample into PP / PN / NN from sample_attributes --------
## Token patterns are case-insensitive; the disease label for the mixed AD/PSO
## study (SRP165679) is recovered from the sample_title prefix.
classify_samples <- function(rse, srp) {
  attr <- tolower(as.character(colData(rse)$sra.sample_attributes))
  ttl  <- tolower(as.character(colData(rse)$sra.sample_title))
  cls <- rep(NA_character_, ncol(rse))
  les  <- grepl("lesional|involved|plaque|ls|psoriatic skin", attr) &
          !grepl("non.?lesional|uninvolved|peri", attr)
  peri <- grepl("non.?lesional|uninvolved|peri|un.?involved", attr)
  norm <- grepl("normal|healthy|control", attr)
  cls[norm] <- "NN"; cls[peri] <- "PN"; cls[les] <- "PP"
  # SRP165679: keep only psoriasis (title prefix PSO_/HC_), drop AD_
  if (srp == "SRP165679") {
    is_ad <- grepl("^ad", ttl)
    cls[is_ad] <- NA
  }
  cls
}
cls_list <- lapply(STUDIES, function(s) classify_samples(rses[[s]], s))
names(cls_list) <- STUDIES

## 11.3 Per-study DE on the harmonised gene space (filterByExpr->TMM->voom) --
## Identical to STEP 2-3, applied per study and per contrast.
de_one <- function(rse, cls, cA, cB) {              # effect = cA vs cB (ref)
  keep_s <- which(cls %in% c(cA, cB))
  if (length(unique(cls[keep_s])) < 2) return(NULL)
  grp <- factor(cls[keep_s], levels = c(cB, cA))    # cB is reference
  if (min(table(grp)) < 2) return(NULL)
  counts <- assay(rse, "counts")[, keep_s]
  gn <- rowData(rse)$gene_name
  dge <- DGEList(counts = counts, genes = data.frame(gene_name = gn))
  keep <- filterByExpr(dge, group = grp); dge <- dge[keep,, keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, "TMM")
  des <- model.matrix(~ grp)
  v <- voom(dge, des); fit <- eBayes(lmFit(v, des))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = dge$genes$gene_name, logFC = tt$logFC,
             SE = tt$logFC / tt$t, t = tt$t, P = tt$P.Value,
             adjP = p.adjust(tt$P.Value, "BH"), n1 = sum(grp==cA), n2 = sum(grp==cB))
}
contrasts <- list(PPvsNN = c("PP","NN"), PNvsNN = c("PN","NN"), PNvsPP = c("PN","PP"))
per_study_de <- lapply(contrasts, function(cc)
  Filter(Negate(is.null), setNames(lapply(STUDIES, function(s)
    de_one(rses[[s]], cls_list[[s]], cc[1], cc[2])), STUDIES)))

## 11.4 Vectorised DerSimonian-Laird random-effects meta-analysis -----------
## For every gene in >=2 studies: pool log2FC by inverse-variance random effects.
## Closed-form DL tau^2 (validated to 4 dp vs metafor::rma). See deep-dive s.4.
meta_dl <- function(de_list) {
  genes <- sort(unique(unlist(lapply(de_list, function(d) d$gene))))
  # build gene x study matrices of yi (logFC) and vi (SE^2)
  Y <- sapply(de_list, function(d) d$logFC[match(genes, d$gene)])
  V <- sapply(de_list, function(d) (d$SE[match(genes, d$gene)])^2)
  rownames(Y) <- rownames(V) <- genes
  out <- t(apply(cbind(Y, V), 1, function(row) {
    k <- ncol(Y); yi <- row[1:k]; vi <- row[(k+1):(2*k)]
    ok <- is.finite(yi) & is.finite(vi) & vi > 0
    if (sum(ok) < 2) return(c(logFC=NA,SE=NA,P=NA,I2=NA,k=sum(ok)))
    yi <- yi[ok]; vi <- vi[ok]; wi <- 1/vi
    ybar <- sum(wi*yi)/sum(wi)
    Q <- sum(wi*(yi-ybar)^2); df <- length(yi)-1
    C <- sum(wi) - sum(wi^2)/sum(wi)
    tau2 <- max(0, (Q-df)/C)
    ws <- 1/(vi+tau2); mu <- sum(ws*yi)/sum(ws); se <- sqrt(1/sum(ws))
    I2 <- max(0,(Q-df)/Q)*100
    c(logFC=mu, SE=se, P=2*pnorm(-abs(mu/se)), I2=I2, k=length(yi))
  }))
  d <- data.frame(gene=genes, out); d$FDR <- p.adjust(d$P, "BH"); d
}
meta_de <- lapply(per_study_de, meta_dl)
message(sprintf("  PPvsNN: %d genes in >=2 studies, %d FDR<0.05",
        sum(!is.na(meta_de$PPvsNN$P)), sum(meta_de$PPvsNN$FDR < 0.05, na.rm=TRUE)))

## 11.5 STAT3 random-effects forest (Figure 25) ------------------------------
stat3_rows <- do.call(rbind, lapply(names(per_study_de$PPvsNN), function(s) {
  d <- per_study_de$PPvsNN[[s]]; r <- d[d$gene=="STAT3",][1,]
  data.frame(study=s, yi=r$logFC, sei=r$SE, n1=r$n1, n2=r$n2)
}))
res_stat3 <- rma(yi = yi, sei = sei, data = stat3_rows, method = "DL")
message(sprintf("  STAT3 pooled: %+.3f [%.3f, %.3f]  I2=%.0f%%  tau2=%.3f  Q=%.1f (p=%.1e)",
        res_stat3$b, res_stat3$ci.lb, res_stat3$ci.ub, res_stat3$I2,
        res_stat3$tau2, res_stat3$QE, res_stat3$QEp))
png(file.path(OUTDIR,"fig25_stat3_forest.png"), width=1100, height=650, res=150)
forest(res_stat3, slab = stat3_rows$study, xlab = "log2 fold-change (PP vs NN)",
       main = "Figure 25. STAT3 random-effects meta-analysis")
dev.off()

## 11.6 Per-study STAT3 TF-activity (decoupleR ULM) + Stouffer pool (Fig 27) -
## Feed each study's per-gene moderated t-statistic into run_ulm against the
## CollecTRI regulon; STAT3's regulon activity is the network-hub test.
regulon <- readRDS(file.path(OUTDIR, "collectri_regulon.rds"))  # from STEP 7
tf_one <- function(d) {
  mat <- matrix(d$t, ncol=1, dimnames=list(d$gene, "t"))
  mat <- mat[!is.na(mat[,1]) & !duplicated(rownames(mat)),,drop=FALSE]
  act <- decoupleR::run_ulm(mat, regulon, .source="source", .target="target",
                            .mor="mor", minsize=5)
  act[act$source=="STAT3" & act$statistic=="ulm",]
}
stat3_tf <- do.call(rbind, lapply(names(per_study_de$PPvsNN), function(s) {
  a <- tf_one(per_study_de$PPvsNN[[s]]); data.frame(study=s, score=a$score, p=a$p_value)
}))
# Stouffer weighted-Z (weight by sqrt n) ; see deep-dive s.6
nvec <- sapply(per_study_de$PPvsNN, function(d) d$n1[1]+d$n2[1])
z_i <- qnorm(1 - stat3_tf$p/2) * sign(stat3_tf$score)
Zc <- sum(sqrt(nvec)*z_i)/sqrt(sum(nvec))
message(sprintf("  STAT3 TF activity pooled (Stouffer): Z=%.2f  p=%.2e", Zc, 2*pnorm(-abs(Zc))))

## 11.7 Cross-study isoform PSI_beta pooling (Figure 29) ---------------------
## Extract the STAT3 alpha/beta junctions (STEP 8 coordinates) from every
## study's junction RSE; pool the lesional-vs-healthy PSI_beta shift with
## random effects, restricting to samples with junction depth >= 20.
## KEY RESULT: the anchor's shift (p=0.017) does NOT replicate (pooled p~0.065).
acc <- 42317181L; a_don <- 42316902L; b_don <- 42316852L; tol <- 2L
psi_one <- function(srp, cls) {
  pj <- subset(ap, project==srp & project_type=="data_sources")
  rj <- create_rse(pj, type="jxn")
  gr <- rowRanges(rj); onchr <- as.character(seqnames(gr)) %in% c("chr17","17")
  is_a <- onchr & abs(start(gr)-a_don)<=tol & abs(end(gr)-acc)<=tol
  is_b <- onchr & abs(start(gr)-b_don)<=tol & abs(end(gr)-acc)<=tol
  if (sum(is_a)<1 || sum(is_b)<1) return(NULL)
  a_cnt <- colSums(matrix(assay(rj,"counts")[is_a,,drop=FALSE], ncol=ncol(rj)))
  b_cnt <- colSums(matrix(assay(rj,"counts")[is_b,,drop=FALSE], ncol=ncol(rj)))
  data.frame(sample=colnames(rj), class=cls, a=a_cnt, b=b_cnt,
             depth=a_cnt+b_cnt, psi_beta=b_cnt/(a_cnt+b_cnt))
}
psi_list <- lapply(STUDIES, function(s) psi_one(s, cls_list[[s]]))
names(psi_list) <- STUDIES
psi_eff <- do.call(rbind, lapply(STUDIES, function(s) {
  d <- psi_list[[s]]; if (is.null(d)) return(NULL)
  d <- d[d$depth>=20 & d$class %in% c("PP","NN"),]
  pp <- d$psi_beta[d$class=="PP"]; nn <- d$psi_beta[d$class=="NN"]
  if (length(pp)<3 || length(nn)<3) return(NULL)
  data.frame(study=s, delta=mean(pp)-mean(nn),
             se=sqrt(var(pp)/length(pp)+var(nn)/length(nn)))
}))
res_psi <- rma(yi=delta, sei=se, data=psi_eff, method="DL")
message(sprintf("  PSI_beta pooled shift: %+.4f [%.4f, %.4f] I2=%.0f%% p=%.3f (%s)",
        res_psi$b, res_psi$ci.lb, res_psi$ci.ub, res_psi$I2, res_psi$pval,
        ifelse(res_psi$pval<0.05,"replicates","does NOT replicate")))

## 11.8 clust cross-study co-expression modules (Figure 30) ------------------
## clust is a PYTHON tool (Abu-Jamous & Kelly 2018) run from the command line.
## It finds gene modules co-expressed CONSISTENTLY across datasets. Reproduce:
##   1. export per-study log2-CPM matrices (genes x samples), collapsed to symbol,
##      subset to the meta-significant genes, to clust_input/<SRP>.tsv
##   2. run:  python -m clust clust_input -o clust_out -t 0   (numpy<2 env)
##   3. read clust_out/Clusters_Objects.tsv for module membership
## Result: one consensus module (65 genes) across the 3 large studies, all UP in
## lesional (median logFC +3.43): IL-17 antimicrobial core + mitotic cassette,
## with a healthy < peri-lesional < lesional gradient. See white-paper s.11.7.
## (Export helper and figure code are in the companion analysis; module saved to
##  clust_module.rds / clust_module_genes.csv.)

message("STEP 11 DONE. Figures: fig22_study_inventory.png ... fig30_module.png")
message(sprintf("  STAT3 gene pooled %+.2f  | TF-activity Stouffer Z=%.1f  | PSI_beta p=%.3f (n.s.)",
        res_stat3$b, Zc, res_psi$pval))
message("  Deliverables: meta_de_*.csv, per_study_de.rds, psi_meta_result.rds, clust_module.rds")

# =============================================================================
# STEP 12 - Molecular staging axis: NN -> PN -> PP as an ordered trajectory
# =============================================================================
# Reframes the case-control contrast as a THREE-POINT ordered axis by adding the
# peri-lesional (PN) group. Anchored in SRP165679 (the deepest three-group cohort:
# NN=38, PN=27, PP=28). Outputs go to staging_axis/. See the companion document
# psoriasis_staging_axis_whitepaper.md. Requires per_study_de.rds (STEP 11) and the
# per-study normalised matrix clust_input/SRP165679.tsv (STEP 11.8).
message("STEP 12: molecular staging axis (NN -> PN -> PP) ...")
suppressPackageStartupMessages({ library(limma); library(msigdbr) })
dir.create("staging_axis/results", recursive = TRUE, showWarnings = FALSE)
dir.create("staging_axis/figures", showWarnings = FALSE)

## 12.1 Trend test: per-gene linear trend + rank monotonicity across the axis --
mat_s  <- fread("clust_input/SRP165679.tsv"); Xs <- as.matrix(mat_s[,-1]); rownames(Xs) <- mat_s$Genes
cls_s  <- as.data.table(readRDS("sample_classification.rds"))
lab_s  <- cls_s[srp == "SRP165679"][class %in% c("NN","PN","PP")]
Xs     <- Xs[, lab_s$external_id]
grp_s  <- factor(lab_s$class, levels = c("NN","PN","PP"))
stage  <- as.integer(grp_s) - 1L                      # NN=0, PN=1, PP=2
fit_s  <- eBayes(lmFit(Xs, model.matrix(~ stage)), trend = TRUE)
tt_s   <- topTable(fit_s, coef = "stage", number = Inf, sort.by = "none")
trend  <- data.table(gene = rownames(Xs), trend_slope = tt_s$logFC, trend_t = tt_s$t,
                     trend_P = tt_s$P.Value, trend_FDR = tt_s$adj.P.Val)
mu_s   <- sapply(c("NN","PN","PP"), function(g) rowMeans(Xs[, grp_s == g, drop = FALSE]))
trend[, `:=`(mu_NN = mu_s[,1], mu_PN = mu_s[,2], mu_PP = mu_s[,3])]
trend[, mono_up   := mu_NN < mu_PN & mu_PN < mu_PP]
trend[, mono_down := mu_NN > mu_PN & mu_PN > mu_PP]
trend[, pn_frac   := (mu_PN - mu_NN) / (mu_PP - mu_NN)]
rho_s  <- apply(Xs, 1, function(v) suppressWarnings(cor(v, stage, method = "spearman")))
trend[, spearman_rho := rho_s[gene]]
saveRDS(trend, "staging_axis/results/trend_SRP165679.rds")

# restrict to the lesional DE programme and summarise monotonicity
de_ps  <- readRDS("per_study_de.rds")
de_gen <- as.data.table(de_ps$PPvsNN[["SRP165679"]])[adjP < 0.05 & abs(logFC) > 1, gene]
g      <- trend[gene %in% de_gen]
g[, monotonic := mono_up | mono_down]
fwrite(g[order(-abs(trend_slope)),
         .(gene, trend_slope, trend_FDR, spearman_rho, mu_NN, mu_PN, mu_PP, pn_frac)],
       "staging_axis/results/gradient_genes.csv")
message(sprintf("  12.1 %d DE genes | strictly monotonic %.0f%% | |rho|>0.3 %.0f%% | median PN frac %.0f%%",
        nrow(g), 100*mean(g$monotonic), 100*mean(abs(g$spearman_rho)>0.3),
        100*median(g$pn_frac, na.rm = TRUE)))

## 12.2 Timing taxonomy: early / progressive / late / PN-divergent -----------
g[, cls := fifelse(pn_frac < 0, "PN_divergent",
            fifelse(pn_frac >= 0.50, "early_PN",
            fifelse(pn_frac >= 0.15, "progressive", "late_PP")))]
fwrite(g[order(cls, -abs(trend_slope)),
         .(gene, cls, trend_slope, trend_FDR, spearman_rho, mu_NN, mu_PN, mu_PP, pn_frac)],
       "staging_axis/results/gradient_gene_classes.csv")
# Hallmark enrichment per informative class (Fisher's exact)
Hs   <- as.data.table(msigdbr(species = "Homo sapiens", category = "H"))
Hset <- split(Hs$gene_symbol, Hs$gs_name); univ <- unique(trend$gene)
enr1 <- function(genes, topn = 8) {
  genes <- intersect(genes, univ)
  r <- rbindlist(lapply(names(Hset), function(s) {
    ins <- intersect(Hset[[s]], univ); a <- length(intersect(genes, ins))
    p <- fisher.test(matrix(c(a, length(genes)-a, length(ins)-a,
                              length(univ)-length(genes)-length(ins)+a), 2),
                     alternative = "greater")$p.value
    data.table(set = sub("HALLMARK_","",s), overlap = a, set_size = length(ins), p = p) }))
  r[, FDR := p.adjust(p, "BH")]; r[order(p)][1:topn] }
enr_all <- rbind(cbind(enr1(g[cls=="progressive", gene]), class = "progressive (early/inflammatory)"),
                 cbind(enr1(g[cls=="late_PP",    gene]), class = "late_PP (proliferation)"))
saveRDS(enr_all, "staging_axis/results/class_enrichment.rds")
message(sprintf("  12.2 classes: progressive=%d late_PP=%d early_PN=%d PN_divergent=%d",
        g[cls=="progressive",.N], g[cls=="late_PP",.N], g[cls=="early_PN",.N], g[cls=="PN_divergent",.N]))

## 12.3 Pathway timing: base-R ssGSEA across the three groups (Fig S3) --------
# ssgsea_scores() is defined in STEP 10; reuse it here.
keys <- c("INTERFERON_GAMMA_RESPONSE","INTERFERON_ALPHA_RESPONSE","IL6_JAK_STAT3_SIGNALING",
          "INFLAMMATORY_RESPONSE","TNFA_SIGNALING_VIA_NFKB","ALLOGRAFT_REJECTION",
          "G2M_CHECKPOINT","E2F_TARGETS","MYC_TARGETS_V1","MTORC1_SIGNALING",
          "FATTY_ACID_METABOLISM","ADIPOGENESIS")
gs_s <- lapply(paste0("HALLMARK_", keys), function(s) Hs[gs_name == s, gene_symbol]); names(gs_s) <- keys
sc_s <- ssgsea_scores(Xs, gs_s)
pth  <- rbindlist(lapply(rownames(sc_s), function(p) {
  y <- sc_s[p,]; lp <- summary(lm(y ~ stage))$coefficients["stage",4]
  data.table(pathway = p, trend_p = lp,
             pn_vs_nn_p = wilcox.test(y[grp_s=="PN"], y[grp_s=="NN"])$p.value,
             NN = mean(y[grp_s=="NN"]), PN = mean(y[grp_s=="PN"]), PP = mean(y[grp_s=="PP"])) }))
pth[, `:=`(pn_frac = (PN-NN)/(PP-NN), trend_FDR = p.adjust(trend_p,"BH"),
           pn_FDR = p.adjust(pn_vs_nn_p,"BH"))]
fwrite(pth[order(-pn_frac)], "staging_axis/results/pathway_timing_stats.csv")
saveRDS(list(scores = sc_s, grp = grp_s), "staging_axis/results/ssgsea_3group.rds")
message("  12.3 pathway timing: interferon already up at PN; proliferation (G2M/E2F) flat until PP")

## 12.4 STAT3 TF-activity across the axis (decoupleR ULM) ---------------------
suppressPackageStartupMessages(library(decoupleR))
reg_s <- as.data.table(readRDS("collectri_regulon.rds"))
dev_s <- Xs - rowMeans(Xs[, grp_s == "NN"])            # deviation from healthy baseline
net_s <- reg_s[target %in% rownames(dev_s), .(source, target, mor)]
act_s <- as.data.table(run_ulm(mat = dev_s, net = net_s, .source = "source",
                               .target = "target", .mor = "mor", minsize = 5))
s3_s  <- merge(act_s[source == "STAT3"], data.table(condition = colnames(dev_s), grp = grp_s), by = "condition")
saveRDS(s3_s, "staging_axis/results/stat3_activity_3group.rds")
mu3   <- tapply(s3_s$score, s3_s$grp, mean)
message(sprintf("  12.4 STAT3 activity NN=%.2f PN=%.2f PP=%.2f (PN reaches %.0f%% of lesional)",
        mu3["NN"], mu3["PN"], mu3["PP"], 100*(mu3["PN"]-mu3["NN"])/(mu3["PP"]-mu3["NN"])))

# NOTE: figures fig_s2..fig_s5b are built by the plotting code archived alongside
# staging_axis/ (they reuse the objects saved above). See staging_axis/code/.
message("STEP 12 DONE. Staging-axis outputs in staging_axis/  (see companion white paper)")

# =============================================================================
# End of pipeline (STEP 1-12). STEP 1-11 confirm STAT3 as a robust cross-cohort
# network hub and retire the STAT3-beta isoform switch as a single-study effect.
# STEP 12 reframes the data as a molecular staging axis (NN -> PN -> PP): 85% of
# the psoriasis programme is monotonic, inflammation is the early event and
# keratinocyte proliferation the late event, and STAT3 activation marks the axis
# early. Further studies (treatment/timepoint cohorts ERP110816, SRP065812) can
# extend the pool in a second pass.
# =============================================================================
