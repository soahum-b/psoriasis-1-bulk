# Auto-extracted generating script
# Produces: fig_significance.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 394ebf21-1958-401f-92c3-0f77e140910e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)

OUT <- "figures_full"
dir.create(OUT, showWarnings = FALSE)

theme_set(theme_bw(base_size = 12))

REPO <- "/net/dali/home/mscbio/sba50/Soahum/Project/psoriasis/psoriasis-1-bulk"

rel  <- readRDS(file.path(REPO, "results_full/reliability_test.rds"))
pn   <- readRDS(file.path(REPO, "results_full/permutation_null.rds"))

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

ggsave(file.path(OUT, "fig_significance.png"), pD1|pD2, width=11, height=4.6, dpi=150)