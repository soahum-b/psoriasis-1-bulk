# Auto-extracted generating script
# Produces: fig38_stat3_vs_module_gradient.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): clust_module_lfc.rds, gradient_gene_classes.csv, clust_module_genes.csv
# Source artifact version: d23535ce-75c3-4e7b-aa9b-4a6d11044109
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(data.table)

theme_set(theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=12), plot.title.position="plot"))
BLUE<-"#4C72B0"; RED<-"#C44E52"; GREY<-"#8C8C8C"; GOLD<-"#DD8452"

gc <- read.csv("gradient_gene_classes.csv")
setDT(gc)

cm <- read.csv("clust_module_genes.csv")

mod_genes <- cm$gene
mg <- gc[gene %in% mod_genes]

norm_traj <- function(dt){
  d <- copy(dt)
  d[, denom := mu_PP - mu_NN]
  d <- d[abs(denom) > 0.2]
  d[, `:=`(nNN=0, nPN=(mu_PN-mu_NN)/denom, nPP=1)]
  d
}
mod_n <- norm_traj(mg)

stat3_npn <- (gc[gene=="STAT3",mu_PN]-gc[gene=="STAT3",mu_NN])/(gc[gene=="STAT3",mu_PP]-gc[gene=="STAT3",mu_NN])

mod_traj <- data.frame(class=c("NN","PN","PP"),
  val=c(0, mean(mg$mu_PN-mg$mu_NN), mean(mg$mu_PP-mg$mu_NN)), series="Headline module (65 genes)")
stat3_traj <- data.frame(class=c("NN","PN","PP"),
  val=c(0, gc[gene=="STAT3",mu_PN]-gc[gene=="STAT3",mu_NN], gc[gene=="STAT3",mu_PP]-gc[gene=="STAT3",mu_NN]),
  series="STAT3")
trajA <- rbind(mod_traj, stat3_traj); trajA$class <- factor(trajA$class, levels=c("NN","PN","PP"))

clab <- c(early_PN="early_PN\n(primers)", progressive="progressive", late_PP="late_PP", PN_divergent="PN_divergent")
gc[, cls_f := factor(clab[cls], levels=clab[c("early_PN","progressive","late_PP","PN_divergent")])]
med_by <- gc[, .(m=median(pn_frac)), by=cls_f]

pA <- ggplot(trajA, aes(class, val, color=series, group=series)) +
  geom_hline(yintercept=0, linetype=2, color=GREY) +
  geom_line(linewidth=1) + geom_point(size=3.5) +
  scale_color_manual(values=c("STAT3"=RED, "Headline module (65 genes)"="#4C72B0"), name=NULL) +
  labs(x=NULL, y="Expression change vs healthy NN (log2)",
       title="A  STAT3 rises later than the psoriasis module",
       subtitle="Both jump into the plaque; peri-lesional (PN) barely moves for either") +
  theme(legend.position=c(0.32,0.85), legend.background=element_blank())

pB <- ggplot(gc[abs(mu_PP-mu_NN)>0.2], aes(pn_frac, cls_f)) +
  geom_vline(xintercept=c(0,1), linetype=3, color=GREY) +
  geom_jitter(height=0.18, alpha=0.12, size=0.7, color="#666666") +
  annotate("rect", xmin=quantile(mod_n$nPN,.25), xmax=quantile(mod_n$nPN,.75), ymin=-Inf, ymax=Inf, alpha=0.12, fill="#4C72B0") +
  geom_vline(xintercept=median(mod_n$nPN), color="#4C72B0", linewidth=0.9) +
  geom_vline(xintercept=stat3_npn, color=RED, linewidth=0.9) +
  geom_point(data=med_by, aes(m, cls_f), color="black", size=3, shape=18) +
  annotate("text", x=median(mod_n$nPN), y=4.45, label="module median 0.19", color="#4C72B0", size=2.8, vjust=1, hjust=-0.05) +
  annotate("text", x=stat3_npn, y=0.7, label="STAT3 0.10", color=RED, size=2.8, vjust=0, hjust=1.1) +
  labs(x="PN position on healthy(0) -> lesional(1) axis", y=NULL,
       title="B  Where peri-lesional sits, by gradient class",
       subtitle="STAT3 (0.10) trails the module (0.19); true early-shift genes sit near 0.6") +
  coord_cartesian(xlim=c(-0.4,1.1))

fig <- pA / pB + plot_layout(heights=c(0.8,1))
ggsave("fig38_stat3_vs_module_gradient.png", fig, width=10.5, height=8.5, dpi=150)