# Auto-extracted generating script
# Produces: fig_alpha_tuning.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: af19feef-aa5e-4251-a5f8-686003733caa
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(Seurat)
library(ggplot2)

OUT <- "figures_full"
dir.create(OUT, showWarnings = FALSE)

theme_set(theme_bw(base_size = 12))

tun <- readRDS("results_full/scissor_tuning.rds")

tt <- tun$tuning
pA <- ggplot(tt, aes(alpha, selected_frac*100)) +
  geom_hline(yintercept=20, linetype=2, color="grey50") +
  geom_line(color="#377eb8") + geom_point(size=3, color="#377eb8") +
  geom_point(data=tt[which.min(abs(tt$selected_frac-0.15)),], size=5, shape=21, fill="#d7301f") +
  labs(x="alpha (graph-smoothing strength)", y="% cells selected",
       title="Full-census alpha tuning", subtitle="chosen alpha=0.20, 14.67% selected (<20% cutoff)")
ggsave(file.path(OUT,"fig_alpha_tuning.png"), pA, width=6, height=4.5, dpi=150)