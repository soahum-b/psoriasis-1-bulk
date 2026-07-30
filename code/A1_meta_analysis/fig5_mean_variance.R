# Auto-extracted generating script
# Produces: fig5_mean_variance.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): dge_filt_norm.rds
# Source artifact version: 93567404-999e-4579-ba2f-130f56c7d525
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR)})
dge <- readRDS("dge_filt_norm.rds")

grp <- dge$samples$group
nn  <- dge[, grp == levels(grp)[1]]

cnt <- nn$counts
mu  <- rowMeans(cnt)
v   <- apply(cnt, 1, var)
keep <- mu > 0
mu <- mu[keep]; v <- v[keep]

nn <- nn[keep, ]
nn <- estimateDisp(nn, design = model.matrix(~1, data = nn$samples))

suppressMessages({library(ggplot2)})

df <- data.frame(mean = mu, var = v)
df <- df[df$mean > 0 & df$var > 0, ]

p <- ggplot(df, aes(mean, var)) +
  geom_point(alpha = 0.15, size = 0.5, colour = "#555555") +
  geom_abline(aes(slope = 1, intercept = 0, colour = "Poisson (var = mean)"),
              linewidth = 1) +
  stat_function(fun = function(x) x + 0.1241 * x^2,
                aes(colour = "Neg. binomial (var = mean + phi*mean^2)"),
                linewidth = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("Poisson (var = mean)" = "#4C72B0",
                                 "Neg. binomial (var = mean + phi*mean^2)" = "#C44E52"),
                      name = NULL) +
  labs(title = "Mean-variance relationship in NN (normal) samples",
       subtitle = "Genes lie far above the Poisson line: RNA-seq counts are over-dispersed",
       x = "Mean count across samples (log scale)",
       y = "Variance across samples (log scale)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("fig5_mean_variance.png", p, width = 8, height = 5.5, dpi = 150)