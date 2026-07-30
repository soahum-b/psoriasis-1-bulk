# Auto-extracted generating script
# Produces: fig31_stat3_hksj_loo.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): Recount-3/per_study_de.rds
# Source artifact version: 623559a2-3266-44c0-8426-0b4ba60bbca5
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")
.libPaths(c(".r-libs/psoriasis-r", .libPaths()))
suppressMessages(library(metafor))
library(ggplot2)
library(patchwork)

per_study_de <- readRDS("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/per_study_de.rds")

# Reproduce the exact per-study STAT3 effect table the forest was built from
stat3_rows <- do.call(rbind, lapply(names(per_study_de$PPvsNN), function(s) {
  d <- per_study_de$PPvsNN[[s]]; r <- d[d$gene=="STAT3",][1,]
  data.frame(study=s, yi=r$logFC, sei=r$SE, n1=r$n1, n2=r$n2)
}))

# --- Original DL (reproduce fig25) ---
res_dl <- rma(yi=yi, sei=sei, data=stat3_rows, method="DL")
# --- HKSJ: REML + Knapp-Hartung t-based CI (guideline default for few, heterogeneous studies) ---
res_hk  <- rma(yi=yi, sei=sei, data=stat3_rows, method="REML", test="knha")
# --- Also DL+KH to isolate the CI method from the tau2 estimator ---
res_dlhk<- rma(yi=yi, sei=sei, data=stat3_rows, method="DL", test="knha")

# Leave-one-out under BOTH the original DL and the HKSJ primary
loo_tbl <- do.call(rbind, lapply(seq_len(nrow(stat3_rows)), function(i){
  sub <- stat3_rows[-i,]
  dl <- rma(yi=yi, sei=sei, data=sub, method="DL")
  hk <- rma(yi=yi, sei=sei, data=sub, method="REML", test="knha")
  data.frame(dropped = stat3_rows$study[i],
             k = nrow(sub),
             DL_est = as.numeric(dl$b), DL_lb = dl$ci.lb, DL_ub = dl$ci.ub, DL_p = dl$pval,
             HK_est = as.numeric(hk$b), HK_lb = hk$ci.lb, HK_ub = hk$ci.ub, HK_p = hk$pval)
}))

# ---- Panel A data: per-study + three pooled estimates ----
per <- data.frame(label=stat3_rows$study, est=stat3_rows$yi,
                  lb=stat3_rows$yi-1.96*stat3_rows$sei, ub=stat3_rows$yi+1.96*stat3_rows$sei,
                  grp="Per-study", n=paste0("n=",stat3_rows$n1,"/",stat3_rows$n2))
pooled <- data.frame(
  label=c("Pooled: DL + z  (original fig25)","Pooled: DL + Knapp-Hartung","Pooled: REML + HKSJ  (recommended)"),
  est=c(res_dl$b,res_dlhk$b,res_hk$b),
  lb =c(res_dl$ci.lb,res_dlhk$ci.lb,res_hk$ci.lb),
  ub =c(res_dl$ci.ub,res_dlhk$ci.ub,res_hk$ci.ub),
  grp="Pooled", n=c("p=3e-08","p=0.049","p=0.059"))
pA <- rbind(per,pooled)
pA$label <- factor(pA$label, levels=rev(pA$label))
pA$grp <- factor(pA$grp, levels=c("Per-study","Pooled"))

gA <- ggplot(pA, aes(est, label, color=grp)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=lb, xmax=ub), height=0.28, linewidth=0.8) +
  geom_point(aes(shape=grp), size=3.4) +
  geom_text(aes(label=n), hjust=0, nudge_y=0.30, size=3, color="grey30") +
  scale_color_manual(values=c("Per-study"="#4477AA","Pooled"="#CC3311"), guide="none") +
  scale_shape_manual(values=c("Per-study"=16,"Pooled"=18), guide="none") +
  labs(x="STAT3 log2 fold-change (PP vs NN)", y=NULL,
       title="A  Pooled STAT3 effect: choice of method decides significance",
       subtitle="Point estimate stable (~+1.2); interval width triples once few-study uncertainty is honoured") +
  coord_cartesian(xlim=c(-2.5,3.2)) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold",size=12),
        plot.subtitle=element_text(size=9,color="grey35"),
        axis.text.y=element_text(size=9.5), panel.grid.minor=element_blank())

# ---- Panel B data: leave-one-out (DL point + CI; direction robustness) ----
pB <- data.frame(
  label=paste0("drop ",loo_tbl$dropped," (k=2)"),
  est=loo_tbl$DL_est, lb=loo_tbl$DL_lb, ub=loo_tbl$DL_ub, p=loo_tbl$DL_p)
# add full model
pB <- rbind(data.frame(label="ALL 3 studies", est=res_dl$b, lb=res_dl$ci.lb, ub=res_dl$ci.ub, p=res_dl$pval), pB)
pB$anchor <- ifelse(grepl("SRP035988", pB$label), "anchor removed","")
pB$label <- factor(pB$label, levels=rev(pB$label))

gB <- ggplot(pB, aes(est, label)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_vline(xintercept=res_dl$b, linetype=3, color="#CC3311") +
  geom_errorbarh(aes(xmin=lb, xmax=ub), height=0.25, color="#228833", linewidth=0.8) +
  geom_point(size=3.2, color="#228833") +
  geom_text(aes(label=ifelse(anchor!="","anchor removed - STAT3 still up","")),
            hjust=0.5, nudge_y=-0.28, size=3, color="#AA3377") +
  labs(x="Pooled STAT3 log2FC (DerSimonian-Laird)", y=NULL,
       title="B  Leave-one-out: effect is not anchor-driven",
       subtitle="Every re-pool stays positive (+0.95 to +1.39); dropping the large anchor keeps STAT3 up (+1.12)") +
  coord_cartesian(xlim=c(-0.3,2.6)) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold",size=12),
        plot.subtitle=element_text(size=9,color="grey35"),
        axis.text.y=element_text(size=9.5), panel.grid.minor=element_blank())

fig <- gA / gB + plot_layout(heights=c(1.15,1))
ggsave("fig31_stat3_hksj_loo.png", fig, width=9.2, height=8.4, dpi=150, bg="white")