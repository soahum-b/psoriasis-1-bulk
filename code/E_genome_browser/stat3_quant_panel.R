# Auto-extracted generating script
# Produces: stat3_quant_panel.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): stat3_gradient_isoform_vs_gene.csv, trend_SRP165679.rds
# Source artifact version: 5ae85cc3-8708-49d6-9ac8-7afbc06df8dd
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

ap <- function(vid) host$artifact_path(vid)

tr <- readRDS(ap("d61263b9-9355-4e92-97d8-4987bcb38e86"))
setDT(tr)

psi <- read.csv(ap("087461f6-2a17-486a-bb9b-4f76963e8035"))

grp_cols <- c(NN="#4DAF4A", PN="#FF7F00", PP="#984EA3")

s3 <- tr[gene=="STAT3"]
d1 <- data.frame(arm=factor(c("NN","PN","PP"),levels=c("NN","PN","PP")),
                 expr=c(s3$mu_NN,s3$mu_PN,s3$mu_PP))
d2 <- data.frame(arm=factor(psi$arm,levels=c("NN","PN","PP")),
                 psi=psi$PSI_beta_SRP165679_pct)

p1 <- ggplot(d1, aes(arm, expr, fill=arm)) +
  geom_col(width=.65) +
  geom_text(aes(label=sprintf("%.2f",expr)), vjust=-0.4, size=3.2) +
  scale_fill_manual(values=grp_cols, guide="none") +
  coord_cartesian(ylim=c(0, max(d1$expr)*1.12)) +
  labs(title="STAT3 gene expression rises with lesion stage",
       subtitle="mean log2 expression, SRP165679 gradient cohort",
       x=NULL, y="log2 expression") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(size=10.5,face="bold"),
        plot.subtitle=element_text(size=8.5), axis.text=element_text(color="black"))

p2 <- ggplot(d2, aes(arm, psi, fill=arm)) +
  geom_col(width=.65) +
  geom_text(aes(label=sprintf("%.1f%%",psi)), vjust=-0.4, size=3.2) +
  scale_fill_manual(values=grp_cols, guide="none") +
  coord_cartesian(ylim=c(0, max(d2$psi)*1.18)) +
  labs(title="STAT3\u03b2 isoform fraction does NOT track the gene",
       subtitle="PSI (percent-spliced-in) of \u03b2, per group",
       x=NULL, y="STAT3\u03b2 PSI (%)") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(size=10.5,face="bold"),
        plot.subtitle=element_text(size=8.5), axis.text=element_text(color="black"))

pq <- p1 | p2
ggsave("stat3_quant_panel.png", pq, width=9, height=3.4, dpi=200)
cat("quant panel saved\n")