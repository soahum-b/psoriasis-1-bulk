# Auto-extracted generating script
# Produces: fig21_upset_leadingedge.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_results_all.rds
# Source artifact version: 5b4a613f-1423-4cb8-b101-abdf3a0b1273
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(edgeR)
  library(limma)
  library(msigdbr)
  library(fgsea)
  library(UpSetR)
})

wslib <- "/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/.r-libs/psoriasis-r"
dir.create(wslib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(wslib, .libPaths()))

NNblue <- "#4C72B0"; PPred <- "#C44E52"; grey <- "#B0B0B0"

# Load GSEA results
gsea_all <- readRDS("gsea_results_all.rds")

# Hallmark collection
h_df <- msigdbr(species="Homo sapiens", category="H")
H <- split(h_df$gene_symbol, h_df$gs_name)

# Get top up-regulated Hallmark pathways by NES
gh <- as.data.frame(gsea_all$H)

# Pick a focused panel: the immune/inflammatory pathways that share the IL-17 core
pick <- c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE",
          "HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_COMPLEMENT")
le <- setNames(lapply(pick, function(p) {
  x <- gh$leadingEdge[gh$pathway==p][[1]]
  if (is.list(x)) x <- x[[1]]
  as.character(x)
}), gsub("HALLMARK_","",pick))

m <- fromList(le)

png("fig21_upset_leadingedge.png", width=1700, height=1000, res=150)
print(upset(m, sets=rev(names(le)), keep.order=TRUE,
      order.by="freq", nintersects=25,
      mainbar.y.label="Genes shared across pathways (leading-edge)",
      sets.x.label="Leading-edge size",
      main.bar.color=PPred, sets.bar.color=NNblue,
      matrix.color="#333333", shade.color=grey,
      text.scale=c(1.3,1.2,1.1,1.0,1.2,1.1),
      point.size=2.6, line.size=0.8))
grid::grid.text("Fig 21. Leading-edge gene overlap across 6 top immune Hallmark pathways",
                x=0.05, y=0.97, just="left", gp=grid::gpar(fontsize=12, fontface="bold"))
grid::grid.text("Shared core (IL6, IRF1, CXCL9/10/11, CCL2/5, IL1B...) counted repeatedly => pathway redundancy",
                x=0.05, y=0.94, just="left", gp=grid::gpar(fontsize=9, col="#555555"))
dev.off()