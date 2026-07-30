# Auto-extracted generating script
# Produces: fig6_volcano.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): dge_filt_norm.rds
# Source artifact version: 8b732492-ca89-41ef-8052-ea9a2c482086
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR); library(limma); library(ggplot2); library(ggrepel)})

dge <- readRDS("dge_filt_norm.rds")
grp <- dge$samples$group

design <- model.matrix(~ grp)
colnames(design) <- c("Intercept", "PPvsNN")

v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "PPvsNN", number = Inf, sort.by = "P")
res$gene <- rownames(res)

res$neglog10 <- -log10(res$P.Value)
res$cat <- "Not sig"
res$cat[res$adj.P.Val < 0.05 & res$logFC >  1] <- "Up in psoriasis"
res$cat[res$adj.P.Val < 0.05 & res$logFC < -1] <- "Down in psoriasis"

cap <- max(res$neglog10[is.finite(res$neglog10)])
res$neglog10[!is.finite(res$neglog10)] <- cap

lab_genes <- c("DEFB4A","S100A7A","S100A8","S100A9","S100A12","SERPINB4","PI3","IL36A",
               "SPRR2A","LCE3A","STAT3","IL17A","IL17F","SOCS3","IL6","KRT77","BTC","RORC","IL34")
res$lab <- ifelse(res$gene %in% lab_genes, res$gene, "")

p <- ggplot(res, aes(logFC, neglog10, colour = cat)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(max(res$P.Value[res$adj.P.Val < 0.05])),
             linetype = "dashed", colour = "grey50") +
  geom_text_repel(aes(label = lab), size = 3, max.overlaps = 30,
                  segment.size = 0.2, show.legend = FALSE) +
  scale_colour_manual(values = c("Up in psoriasis" = "#C44E52",
                                 "Down in psoriasis" = "#4C72B0",
                                 "Not sig" = "grey75"), name = NULL) +
  labs(title = "Differential expression: lesional (PP) vs normal (NN) skin",
       subtitle = "limma-voom; 3,477 genes at adj.P<0.05 and |log2FC|>1 (1,511 up, 1,966 down)",
       x = "log2 fold-change (PP vs NN)", y = "-log10 raw p-value") +
  theme_bw(base_size = 12) + theme(legend.position = "top")

ggsave("fig6_volcano.png", p, width = 9, height = 6.5, dpi = 150)