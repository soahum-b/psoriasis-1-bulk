# Auto-extracted generating script
# Produces: voom_fit.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): dge_filt_norm.rds
# Source artifact version: 12fadcf7-d139-42dc-8020-c9eb75feb54c
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(edgeR); library(limma)})

dge <- readRDS("dge_filt_norm.rds")

grp <- factor(dge$samples$group, levels=c("NN","PP"))
design <- model.matrix(~grp)
colnames(design) <- c("Intercept","PPvsNN")

v <- voom(dge, design, plot=FALSE)

saveRDS(v, "voom_fit.rds")