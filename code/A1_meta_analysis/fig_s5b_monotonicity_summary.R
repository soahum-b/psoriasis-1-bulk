# Auto-extracted generating script
# Produces: fig_s5b_monotonicity_summary.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 7cedf0e1-c765-438f-987c-e8f63762b179
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(ggplot2)
library(patchwork)

g <- readRDS("{{artifact:MISSING:art_5a3ef6bb-6f03-4df8-8c3f-aa9bd2e73bd1}}")
g[, monotonic_flag := ifelse(mono_up,"up", ifelse(mono_down,"down","non"))]

sm <- data.table(
  cat=c("Strictly\nmonotonic","Monotonic\nrank (|rho|>0.3)","Non-\nmonotonic"),
  pct=c(100*g[monotonic_flag!="non",.N]/nrow(g),
        100*g[abs(spearman_rho)>0.3,.N]/nrow(g),
        100*g[monotonic_flag=="non",.N]/nrow(g)),
  col=c("#55A868","#4C72B0","#B0B0B0"))
sm[, cat:=factor(cat, levels=cat)]
pA <- ggplot(sm, aes(cat, pct, fill=col))+geom_col(width=.7)+scale_fill_identity()+
  geom_text(aes(label=sprintf("%.0f%%",pct)), vjust=-0.3, size=4.2)+
  labs(title="A. Monotonicity of the psoriasis program",
       subtitle=sprintf("%d lesional DE genes tested across NN<PN<PP", nrow(g)),
       x=NULL, y="% of DE genes")+
  scale_y_continuous(expand=expansion(mult=c(0,.12)), limits=c(0,105))+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
       plot.title=element_text(face="bold",size=13), plot.subtitle=element_text(size=9.5))

pf <- g[monotonic_flag!="non" & pn_frac>=0 & pn_frac<=1]
med <- median(pf$pn_frac)
pB <- ggplot(pf, aes(pn_frac))+
  geom_histogram(bins=40, fill="#DD8452", colour="white", linewidth=.2)+
  geom_vline(xintercept=med, linetype=2, colour="#C44E52", linewidth=.7)+
  annotate("text", x=med+0.02, y=Inf, vjust=1.6, hjust=0,
           label=sprintf("median %.0f%%", 100*med), colour="#C44E52", size=4)+
  labs(title="B. Peri-lesional sits early on the axis",
       subtitle="Fraction of the lesional change already reached at the PN stage",
       x="PN fraction toward lesional  (0 = healthy-like, 1 = lesional-like)", y="genes")+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
       plot.title=element_text(face="bold",size=13), plot.subtitle=element_text(size=9.5))

fig <- pA + pB + plot_layout(widths=c(1,1.3))
ggsave("fig_s5b_monotonicity_summary.png", fig, width=11, height=4.4, dpi=150)