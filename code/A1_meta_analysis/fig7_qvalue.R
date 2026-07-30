# Auto-extracted generating script
# Produces: fig7_qvalue.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: d74fa3df-c6ca-4823-a7b3-669c36d13386
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(qvalue); library(ggplot2); library(patchwork)})

de <- readRDS("de_results_full.rds")
de$qvalue <- NULL
qobj <- qvalue(p = de$P.Value)
de$qvalue <- qobj$qvalues
cat("pi0 =", round(qobj$pi0,4), "| genes:", nrow(de), "\n")

pA <- ggplot(de, aes(x = P.Value)) +
  geom_histogram(boundary = 0, bins = 50, fill = "grey70", colour = "white") +
  geom_hline(yintercept = qobj$pi0 * nrow(de) / 50, colour = "#C44E52",
             linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.6, y = qobj$pi0 * nrow(de) / 50,
           label = sprintf("pi0 = %.3f  (null level)", qobj$pi0),
           vjust = -0.6, colour = "#C44E52", size = 4) +
  labs(title = "A. p-value distribution",
       subtitle = "Spike near 0 = real signal; flat part = null genes",
       x = "p-value (PP vs NN)", y = "gene count") +
  theme_bw(base_size = 13)

pB <- ggplot(de, aes(x = adj.P.Val, y = qvalue)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dotted") +
  geom_point(alpha = 0.25, size = 0.5, colour = "#4C72B0") +
  geom_vline(xintercept = 0.05, colour = "#C44E52", linetype = "dashed") +
  geom_hline(yintercept = 0.05, colour = "#C44E52", linetype = "dashed") +
  labs(title = "B. Storey q-value vs BH adjusted p",
       subtitle = "q = pi0 x BH, so q <= BH (more powerful)",
       x = "BH adjusted p (adj.P.Val)", y = "Storey q-value") +
  theme_bw(base_size = 13)

ggsave("fig7_qvalue.png", pA + pB, width = 11, height = 4.6, dpi = 150)