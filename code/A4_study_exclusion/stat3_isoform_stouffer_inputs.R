# Auto-extracted generating script
# Produces: stat3_isoform_stouffer_inputs.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 03c72213-ae46-4460-88b5-b2f7c7776db0
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(metafor)})

# Per-study STAT3-beta PP-vs-NN inputs (depth>=20), from the junction tests
iso <- data.frame(
  study  = c("SRP035988","SRP165679","SRP126422","SRP154474"),
  shift  = c( 1.167, -0.071,  1.008,  2.209),
  p      = c( 0.017,  0.700,  0.860,  0.370),
  n      = c( 178,     66,     7,      17))
iso$sign <- sign(iso$shift)

stouffer_wz <- function(p, effect, n, zcap=38){
  z <- qnorm(1 - p/2) * sign(effect)
  z[z >  zcap] <-  zcap; z[z < -zcap] <- -zcap
  Z <- sum(sqrt(n) * z) / sqrt(sum(n))
  list(Z=Z, p=2*(1-pnorm(abs(Z))), z=z)
}

sw <- stouffer_wz(iso$p, iso$shift, iso$n)
iso$signed_z <- round(sw$z, 3)
iso$contrib <- sqrt(iso$n)*iso$signed_z/sqrt(sum(iso$n))

iso_out <- iso[,c("study","n","shift","p","sign","signed_z","contrib")]
write.csv(iso_out, "stat3_isoform_stouffer_inputs.csv", row.names=FALSE)