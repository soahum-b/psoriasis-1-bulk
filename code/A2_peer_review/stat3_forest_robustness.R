# Auto-extracted generating script
# Produces: stat3_forest_robustness.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): Recount-3/per_study_de.rds
# Source artifact version: 106c7130-cfde-42b8-b8ab-ad9b3944a13a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")
.libPaths(c(".r-libs/psoriasis-r", .libPaths()))
suppressMessages(library(metafor))

per_study_de <- readRDS("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/per_study_de.rds")

stat3_rows <- do.call(rbind, lapply(names(per_study_de$PPvsNN), function(s) {
  d <- per_study_de$PPvsNN[[s]]; r <- d[d$gene=="STAT3",][1,]
  data.frame(study=s, yi=r$logFC, sei=r$SE, n1=r$n1, n2=r$n2)
}))

res_dl <- rma(yi=yi, sei=sei, data=stat3_rows, method="DL")
res_hk  <- rma(yi=yi, sei=sei, data=stat3_rows, method="REML", test="knha")
res_dlhk<- rma(yi=yi, sei=sei, data=stat3_rows, method="DL", test="knha")

loo_tbl <- do.call(rbind, lapply(seq_len(nrow(stat3_rows)), function(i){
  sub <- stat3_rows[-i,]
  dl <- rma(yi=yi, sei=sei, data=sub, method="DL")
  hk <- rma(yi=yi, sei=sei, data=sub, method="REML", test="knha")
  data.frame(dropped = stat3_rows$study[i],
             k = nrow(sub),
             DL_est = as.numeric(dl$b), DL_lb = dl$ci.lb, DL_ub = dl$ci.ub, DL_p = dl$pval,
             HK_est = as.numeric(hk$b), HK_lb = hk$ci.lb, HK_ub = hk$ci.ub, HK_p = hk$pval)
}))

robust_tbl <- data.frame(
  analysis = c("DL + z (original)","DL + Knapp-Hartung","REML + HKSJ (recommended)",
               "LOO drop SRP035988 (anchor)","LOO drop SRP165679","LOO drop SRP126422"),
  k = c(3,3,3,2,2,2),
  est = c(res_dl$b,res_dlhk$b,res_hk$b, loo_tbl$DL_est),
  ci_lb = c(res_dl$ci.lb,res_dlhk$ci.lb,res_hk$ci.lb, loo_tbl$DL_lb),
  ci_ub = c(res_dl$ci.ub,res_dlhk$ci.ub,res_hk$ci.ub, loo_tbl$DL_ub),
  p = c(res_dl$pval,res_dlhk$pval,res_hk$pval, loo_tbl$DL_p),
  tau2 = c(res_dl$tau2,res_dlhk$tau2,res_hk$tau2, NA,NA,NA),
  I2  = c(res_dl$I2,res_dlhk$I2,res_hk$I2, NA,NA,NA))
robust_tbl[,3:8] <- round(robust_tbl[,3:8],4)
write.csv(robust_tbl, "stat3_forest_robustness.csv", row.names=FALSE)