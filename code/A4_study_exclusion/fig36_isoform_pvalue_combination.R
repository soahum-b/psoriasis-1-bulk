# Auto-extracted generating script
# Produces: fig36_isoform_pvalue_combination.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 395b84d2-9286-4516-a733-d0eeb7c35f81
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(ggplot2); library(patchwork); library(data.table)})
theme_set(theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=12), plot.title.position="plot"))
BLUE<-"#4C72B0"; RED<-"#C44E52"; GREY<-"#8C8C8C"; GOLD<-"#DD8452"

stouffer_wz <- function(p, effect, n, zcap=38){
  z <- qnorm(1 - p/2) * sign(effect)
  z[z >  zcap] <-  zcap; z[z < -zcap] <- -zcap
  Z <- sum(sqrt(n) * z) / sqrt(sum(n))
  list(Z=Z, p=2*(1-pnorm(abs(Z))), z=z)
}
fisher_comb <- function(p){
  X <- -2*sum(log(p)); df <- 2*length(p)
  list(X=X, p=pchisq(X, df, lower.tail=FALSE))
}

iso <- data.frame(
  study  = c("SRP035988","SRP165679","SRP126422","SRP154474"),
  shift  = c( 1.167, -0.071,  1.008,  2.209),
  p      = c( 0.017,  0.700,  0.860,  0.370),
  n      = c( 178,     66,     7,      17))
iso$sign <- sign(iso$shift)

sw <- stouffer_wz(iso$p, iso$shift, iso$n)
fi <- fisher_comb(iso$p)
iso$signed_z <- round(sw$z, 3); iso$weight_sqrtn <- round(sqrt(iso$n),2)

sw_unw <- { z <- qnorm(1-iso$p/2)*iso$sign; list(Z=sum(z)/sqrt(length(z)), p=2*(1-pnorm(abs(sum(z)/sqrt(length(z)))))) }

iso$contrib <- sqrt(iso$n)*iso$signed_z/sqrt(sum(iso$n))
iso$lab_st <- sprintf("%s\n(n=%d, p=%.2f)", iso$study, iso$n, iso$p)
iso2 <- iso[order(iso$contrib),]; iso2$lab_st <- factor(iso2$lab_st, levels=iso2$lab_st)
pA <- ggplot(iso2, aes(contrib, lab_st, fill=sign>0)) +
  geom_col(width=0.6) +
  geom_vline(xintercept=0, color=GREY) +
  geom_text(aes(label=sprintf("%+.2f", contrib)), hjust=ifelse(iso2$contrib>0,-0.2,1.2), size=3) +
  scale_fill_manual(values=c(`TRUE`=RED,`FALSE`=BLUE), guide="none") +
  labs(x="contribution to combined Z (sqrt(n)-weighted signed z)", y=NULL,
       title="A  Each cohort's contribution to the combined Z",
       subtitle=sprintf("Stouffer weighted-Z = %+.2f, combined p = %.3f  |  anchor SRP035988 supplies %.2f of %.2f (%.0f%%)",
    sw$Z, sw$p, iso$contrib[iso$study=="SRP035988"], sw$Z,
    100*iso$contrib[iso$study=="SRP035988"]/sw$Z)) +
  coord_cartesian(xlim=c(-0.6,1.6))

loo <- data.frame(
  label=c("Stouffer weighted-Z (all 4)","  drop SRP035988 (anchor)","  drop SRP165679","  drop SRP126422","  drop SRP154474",
          "Unweighted Stouffer","Fisher (unsigned)"),
  p=c(sw$p, 0.9133, 0.0113, 0.0449, 0.0655, sw_unw$p, fi$p),
  grp=c("main","loo","loo","loo","loo","alt","alt"))
loo$label <- factor(loo$label, levels=rev(loo$label))
loo$sig <- loo$p<0.05
pB <- ggplot(loo, aes(-log10(p), label, color=grp, shape=sig)) +
  geom_vline(xintercept=-log10(0.05), linetype=2, color=GREY) +
  geom_point(size=3.5) +
  geom_text(aes(label=sprintf("p=%.3f", p)), hjust=-0.25, size=2.9, show.legend=FALSE) +
  annotate("text", x=-log10(0.05), y=0.4, label="p=0.05", size=2.8, color=GREY, vjust=1) +
  scale_color_manual(values=c(main="black",loo=GOLD,alt="#4C72B0"), guide="none") +
  scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1), name="p<0.05") +
  labs(x="-log10(combined p)", y=NULL,
       title="B  The combined signal is anchor-driven and fragile",
       subtitle="Removing the anchor collapses it (p=0.91); it is significant only when SRP035988 is in") +
  coord_cartesian(xlim=c(0,2.6)) +
  theme(legend.position=c(0.85,0.2), legend.background=element_blank(),
        axis.text.y=element_text(size=8.5, family="mono"))

figP <- pA / pB + plot_layout(heights=c(1,1.25))
ggsave("fig36_isoform_pvalue_combination.png", figP, width=10, height=8, dpi=150)