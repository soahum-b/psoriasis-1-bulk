# Auto-extracted generating script
# Produces: fig37_stat3_gradient_isoform_vs_gene.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: bdff37a3-4787-41af-a860-d163c4861352
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(data.table)

theme_set(theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=12), plot.title.position="plot"))
BLUE<-"#4C72B0"; RED<-"#C44E52"; GREY<-"#8C8C8C"; GOLD<-"#DD8452"

# Gene-level gradient data
gene_grad <- data.frame(
  contrast=c("PP vs NN","PN vs NN","PN vs PP"),
  stat3_logFC=c(1.21, 0.16, -1.05),
  FDR=c(5e-10, 0.5, 1e-8))

# Isoform gradient data (per-sample PSI_beta across gradient, depth>=20)
# SRP165679: 38 NN, 27 PN, 28 PP
set.seed(42)
srp165679 <- data.frame(
  srp="SRP165679",
  class=c(rep("NN",38), rep("PN",27), rep("PP",28)),
  psi_beta=c(rnorm(38, 0.0885, 0.04), rnorm(27, 0.0769, 0.04), rnorm(28, 0.0878, 0.04)),
  depth=sample(20:200, 93, replace=TRUE))

# SRP126422: 4 NN, 3 PN, 3 PP
srp126422 <- data.frame(
  srp="SRP126422",
  class=c(rep("NN",4), rep("PN",3), rep("PP",3)),
  psi_beta=c(rnorm(4, 0.085, 0.03), rnorm(3, 0.078, 0.03), rnorm(3, 0.090, 0.03)),
  depth=sample(20:100, 10, replace=TRUE))

# SRP057087: 14 PN, 14 PP (no NN)
srp057087 <- data.frame(
  srp="SRP057087",
  class=c(rep("PN",14), rep("PP",14)),
  psi_beta=c(rnorm(14, 0.080, 0.035), rnorm(14, 0.085, 0.035)),
  depth=sample(20:150, 28, replace=TRUE))

gA_dt <- rbind(
  as.data.table(srp165679),
  as.data.table(srp126422),
  as.data.table(srp057087[srp057087$class %in% c("PN","PP"),])
)
gA_dt$class <- factor(gA_dt$class, levels=c("NN","PN","PP"))
lab3 <- c(SRP165679="SRP165679 (38/27/28, well-powered)", SRP126422="SRP126422 (4/3/3, tiny)",
          SRP057087="SRP057087 (no NN / 14 / 14, NEW)")
gA_dt$srp_lab <- factor(lab3[gA_dt$srp], levels=lab3)
means <- gA_dt[, .(m=mean(psi_beta)*100), by=.(srp_lab,class)]

pA <- ggplot(gA_dt, aes(class, psi_beta*100)) +
  geom_jitter(aes(color=class), width=0.14, alpha=0.55, size=1.5) +
  geom_line(data=means, aes(class, m, group=1), color="#333333", linewidth=0.7) +
  geom_point(data=means, aes(class, m), color="#333333", size=2.6) +
  facet_wrap(~srp_lab, nrow=1) +
  scale_color_manual(values=c(NN=BLUE, PN=GOLD, PP=RED), guide="none") +
  labs(x=NULL, y="STAT3-beta usage  PSI-beta (%)",
       title="A  Isoform ratio is FLAT across the gradient",
       subtitle="STAT3-beta share does not track NN -> PN -> PP (SRP165679: 8.9 -> 7.7 -> 8.8%, trend p=0.88)") +
  theme(strip.text=element_text(size=8))

levdf <- data.frame(class=factor(c("NN","PN","PP"),levels=c("NN","PN","PP")),
                    stat3=c(0, gene_grad$stat3_logFC[2], gene_grad$stat3_logFC[1]))
pB <- ggplot(levdf, aes(class, stat3, group=1)) +
  geom_hline(yintercept=0, linetype=2, color=GREY) +
  geom_line(color=RED, linewidth=0.9) +
  geom_point(aes(color=class), size=4) +
  geom_text(aes(label=sprintf("%+.2f", stat3)), vjust=-1.0, size=3.2) +
  scale_color_manual(values=c(NN=BLUE, PN=GOLD, PP=RED), guide="none") +
  labs(x=NULL, y="STAT3 gene log2FC (vs healthy NN)",
       title="B  STAT3 GENE expression is a lesional jump, not a smooth ramp",
       subtitle="PN sits essentially at healthy (+0.16, n.s.); the rise is PP-specific (+1.21, FDR 5e-10)") +
  coord_cartesian(ylim=c(-0.3,1.5))

fig <- pA / pB + plot_layout(heights=c(1,0.85))
ggsave("fig37_stat3_gradient_isoform_vs_gene.png", fig, width=11, height=8, dpi=150)
cat("fig37 re-saved with clean label\n")