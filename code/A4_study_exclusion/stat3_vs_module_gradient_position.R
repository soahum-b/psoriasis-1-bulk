# Auto-extracted generating script
# Produces: stat3_vs_module_gradient_position.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): clust_module_lfc.rds, gradient_gene_classes.csv, clust_module_genes.csv
# Source artifact version: 949ad6e0-2a8a-491a-8e36-2d8f10597739
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)

gc <- read.csv("gradient_gene_classes.csv")
setDT(gc)

cm <- read.csv("clust_module_genes.csv")

mod_genes <- cm$gene
mg <- gc[gene %in% mod_genes]

norm_traj <- function(dt){
  d <- copy(dt)
  d[, denom := mu_PP - mu_NN]
  d <- d[abs(denom) > 0.2]
  d[, nPN := (mu_PN - mu_NN) / denom]
  d
}
mod_n <- norm_traj(mg)

stat3_npn <- (gc[gene == "STAT3", mu_PN] - gc[gene == "STAT3", mu_NN]) /
             (gc[gene == "STAT3", mu_PP] - gc[gene == "STAT3", mu_NN])

clab <- c(early_PN = "early_PN\n(primers)", progressive = "progressive",
          late_PP = "late_PP", PN_divergent = "PN_divergent")
gc[, cls_f := factor(clab[cls], levels = clab[c("early_PN", "progressive", "late_PP", "PN_divergent")])]
med_by <- gc[, .(m = median(pn_frac)), by = cls_f]

out <- data.frame(
  series = c("STAT3", "Headline module (median)", "early_PN class (median)", "progressive (median)", "late_PP (median)"),
  pn_position_0healthy_1lesional = round(c(
    stat3_npn,
    median(mod_n$nPN),
    med_by$m[med_by$cls_f == "early_PN\n(primers)"],
    med_by$m[med_by$cls_f == "progressive"],
    med_by$m[med_by$cls_f == "late_PP"]
  ), 3),
  interpretation = c(
    "late / plaque-specific",
    "late-ish",
    "early primer (uninvolved already shifting)",
    "modest early shift",
    "late / plaque-specific"
  )
)

write.csv(out, "stat3_vs_module_gradient_position.csv", row.names = FALSE)