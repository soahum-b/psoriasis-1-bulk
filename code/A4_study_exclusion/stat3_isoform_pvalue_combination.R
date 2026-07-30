# Auto-extracted generating script
# Produces: stat3_isoform_pvalue_combination.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 8db0b0e8-9347-4418-8935-919c6d22fc19
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(data.table)
library(metafor)

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

comb_out <- data.frame(
  method=c("Stouffer weighted-Z (sqrt n)","Unweighted Stouffer","Fisher (unsigned)",
           "LOO drop anchor","LOO drop SRP165679","LOO drop SRP126422","LOO drop SRP154474"),
  combined_p=c(sw$p, sw_unw$p, fi$p, 0.9133,0.0113,0.0449,0.0655))
write.csv(comb_out, "stat3_isoform_pvalue_combination.csv", row.names=FALSE)