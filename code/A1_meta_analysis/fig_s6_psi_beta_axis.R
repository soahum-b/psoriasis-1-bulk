# Auto-extracted generating script
# Produces: fig_s6_psi_beta_axis.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): psi_beta_allstudies.rds, sample_classification.rds
# Source artifact version: 9407b5f5-0092-40ba-93e8-7bbb77ab1af0
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({library(data.table); library(ggplot2)})

d <- as.data.table(readRDS("psi_beta_allstudies.rds"))

d20 <- d[depth>=20]

s <- d20[srp=="SRP165679"]
s[, class:=factor(class, levels=c("NN","PN","PP"))]

cols <- c(NN="#4C72B0", PN="#DD8452", PP="#C44E52")
p <- ggplot(s, aes(class, 100*psi_beta, fill=class))+
  geom_boxplot(width=.6, outlier.size=.7)+
  geom_jitter(width=.12, size=1, alpha=.5)+
  scale_fill_manual(values=cols)+
  labs(title="STAT3 alpha/beta splice ratio does NOT follow the staging axis",
       subtitle="SRP165679 (depth>=20): PSI_beta flat across stages (Spearman rho=0.02, p=0.86)\ncontrast: 85% of EXPRESSION genes are monotonic along the same axis",
       x="Staging axis (NN -> PN -> PP)", y="PSI_beta = STAT3-beta / (alpha+beta)  (%)")+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(), legend.position="none",
       plot.title=element_text(face="bold",size=13), plot.subtitle=element_text(size=9.5))
ggsave("fig_s6_psi_beta_axis.png", p, width=6.5, height=5, dpi=150)