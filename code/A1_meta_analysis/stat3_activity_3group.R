# Auto-extracted generating script
# Produces: stat3_activity_3group.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): collectri_regulon.rds, sample_classification.rds
# Source artifact version: 324c85d2-d264-4a56-94bd-9c97fd46b64b
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({library(decoupleR); library(data.table)})

mat <- fread("clust_input/SRP165679.tsv")
X <- as.matrix(mat[,-1]); rownames(X) <- mat$Genes

cls <- as.data.table(readRDS("sample_classification.rds"))
lab <- cls[srp=="SRP165679"][class %in% c("NN","PN","PP")]
X <- X[, lab$external_id]; grp <- factor(lab$class, levels=c("NN","PN","PP"))

reg <- readRDS("collectri_regulon.rds"); reg <- as.data.table(reg)

nnmean <- rowMeans(X[, grp=="NN"])
dev <- X - nnmean
net <- reg[, .(source, target, mor)]
net <- net[target %in% rownames(dev)]
act <- run_ulm(mat=dev, net=net, .source="source", .target="target", .mor="mor", minsize=5)
setDT(act)
s3 <- act[source=="STAT3"]
s3 <- merge(s3, data.table(condition=colnames(dev), grp=grp), by="condition")

saveRDS(s3, "stat3_activity_3group.rds")