# Auto-extracted generating script
# Produces: fig_ref_composition.png
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: fa7106af-da2b-4f14-9a7a-bf4ab170158a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)

so <- readRDS("scissor_repo/results/reference_raw.rds")
md <- so@meta.data

pal <- c(NN="#4C9BD4", PN="#E8A33D", PP="#C43C4E")

tab <- as.data.frame(table(md$condition)); names(tab) <- c("condition","cells")
tab$condition <- factor(tab$condition, levels=c("NN","PN","PP"))
pA <- ggplot(tab, aes(condition, cells, fill=condition)) +
  geom_col(width=.7) +
  geom_text(aes(label=formatC(cells, big.mark=",", format="d")), vjust=-0.3, size=3.2) +
  scale_fill_manual(values=pal, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,.12))) +
  labs(title="Reference spans all three ordinal tiers", x=NULL, y="Cells") +
  theme_classic(base_size=11)

db <- as.data.frame(table(md$donor, md$condition)); names(db) <- c("donor","condition","cells")
db <- db[db$cells>0,]
ord <- names(sort(tapply(md$donor, md$donor, length)))
db$donor <- factor(db$donor, levels=ord)
db$condition <- factor(db$condition, levels=c("NN","PN","PP"))
pB <- ggplot(db, aes(cells, donor, fill=condition)) +
  geom_col() +
  scale_fill_manual(values=pal, name="Tier") +
  scale_x_continuous(expand=expansion(mult=c(0,.05))) +
  labs(title="Cells per donor", x="Cells", y="Donor") +
  theme_classic(base_size=10) +
  theme(legend.position=c(.85,.25), legend.background=element_blank(),
        axis.text.y=element_text(size=6))

fig <- pA + pB + plot_layout(widths=c(1,1.3)) +
  plot_annotation(tag_levels="a")
ggsave("fig_ref_composition.png", fig, width=10, height=4.5, dpi=200)