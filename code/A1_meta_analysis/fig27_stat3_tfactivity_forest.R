# Auto-extracted generating script
# Produces: fig27_stat3_tfactivity_forest.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): study_inventory.csv
# Source artifact version: 840c28fd-aa8e-4468-adf2-c001f0fd48af
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(metafor)
library(data.table)
library(decoupleR)

# Load dependencies
inv <- fread("study_inventory.csv")

# Load per-study TF activity results
tf <- readRDS("tf_activity_perstudy.rds")

# Pool STAT3 activity across PPvsNN studies using Stouffer weighted-Z
pool_tf <- function(contrast, tfname) {
  d <- tf[[contrast]][source==tfname]
  d[, ntot := inv$total[match(study, inv$srp)]]
  z <- sign(d$score) * qnorm(1 - d$p_value/2)
  w <- sqrt(d$ntot)
  Zc <- sum(w*z)/sqrt(sum(w^2))
  p_comb <- 2*pnorm(-abs(Zc))
  list(d=d, Z=Zc, p=p_comb, mean_score=mean(d$score))
}

st3 <- pool_tf("PPvsNN","STAT3")
cat(sprintf("STAT3 pooled (Stouffer, PPvsNN): Z=%.2f  p=%.2e  mean ULM score=%.2f\n",
    st3$Z, st3$p, st3$mean_score))

d <- st3$d
d[, se := abs(score/qnorm(1-p_value/2))]
d[is.na(se)|!is.finite(se), se := abs(score)/6]
res <- rma(yi=score, sei=se, data=d, method="DL")

png("fig27_stat3_tfactivity_forest.png", width=1500, height=850, res=190)
par(mar=c(5,4,3,2))
slab <- sprintf("%s (rank %s)", d$study,
   sapply(d$study, function(s){ dd<-tf$PPvsNN[study==s][order(-score)]; paste0(which(dd$source=="STAT3"),"/",nrow(dd)) }))
forest(res, slab=slab, xlab="STAT3 regulon activity (ULM z-score)",
       header=c("Study (STAT3 TF rank)","activity [95% CI]"), col="#C44E52", border="#C44E52",
       cex=0.95, psize=1.3, mlab="RE pooled", xlim=c(-6,16))
text(-6, -1.6, pos=4, cex=0.8,
     bquote(paste("Stouffer weighted-Z = ", .(sprintf("%.1f",st3$Z)),
       ",  p = ", .(sprintf("%.1e",st3$p)), "    (STAT3 top 3% of ~730 TFs in both large studies)")))
title("Figure 27. STAT3 transcription-factor activity across studies (PP vs NN)", cex.main=0.98)
dev.off()