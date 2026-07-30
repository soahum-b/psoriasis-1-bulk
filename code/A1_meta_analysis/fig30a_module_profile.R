# Auto-extracted generating script
# Produces: fig30a_module_profile.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): per_study_de.rds, study_inventory.csv
# Source artifact version: f739aa0c-c10c-4655-b11e-39afa82c04b9
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(edgeR)
library(limma)
library(SummarizedExperiment)
library(data.table)
library(recount3)
library(ggplot2)
library(patchwork)
library(metafor)

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/")

de <- readRDS("per_study_de.rds")

NNblue <- "#4C72B0"
PPred <- "#C44E52"

# Read clust module members
lines <- readLines("clust_out3/Clusters_Objects.tsv")
tab <- fread("clust_out3/Clusters_Objects.tsv", skip = 1, header = TRUE)
mod <- tab[[1]]
mod <- mod[mod != "" & !is.na(mod)]

# Load cls_all
cls_all <- readRDS("/tmp/cls_all.rds")
setDT(cls_all)

studies3 <- c("SRP035988", "SRP165679", "SRP076982")

prof <- rbindlist(lapply(studies3, function(s) {
  m <- fread(sprintf("clust_input2/%s.tsv", s))
  mm <- m[Genes %in% mod]
  mat <- as.matrix(mm[, -1])
  rownames(mat) <- mm$Genes
  z <- t(scale(t(mat)))
  score <- colMeans(z, na.rm = TRUE)
  cmap <- cls_all[cls_all$srp == s]
  cl <- cmap$class[match(colnames(mat), cmap$external_id)]
  data.table(study = s, sample = colnames(mat), class = cl, score = score)
}))

prof <- prof[class %in% c("PP", "PN", "NN")]
prof[, class := factor(class, levels = c("NN", "PN", "PP"))]
cols <- c("NN" = NNblue, "PN" = "#DD8452", "PP" = PPred)

pA <- ggplot(prof, aes(class, score, color = class)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_wrap(~study, nrow = 1, scales = "free_x") +
  scale_color_manual(values = cols, name = NULL) +
  labs(x = NULL, y = "module z-score (65 genes)",
       title = "Figure 30. Cross-study co-expression module (clust C0)",
       subtitle = "65 genes co-expressed across all 3 large studies; up in lesional. IL-17 antimicrobial core + mitotic cassette.") +
  theme_bw(base_size = 9) +
  theme(plot.subtitle = element_text(size = 6.8),
        strip.text = element_text(size = 7.5),
        legend.position = "none")

ggsave("fig30a_module_profile.png", pA, width = 10, height = 3.4, dpi = 150)