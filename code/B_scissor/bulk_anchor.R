# Auto-extracted generating script
# Produces: bulk_anchor.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): rse_SRP165679.rds, Recount-3/sample_classification.rds
# Source artifact version: a7a850eb-ad16-4c77-91bd-178aadea75f0
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(recount3); library(SummarizedExperiment); library(edgeR)})
rse <- readRDS("rse_SRP165679.rds")
sc  <- readRDS("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3/sample_classification.rds")
sc165 <- sc[sc$srp=="SRP165679" & sc$class %in% c("NN","PN","PP"), ]
cat("classified NN/PN/PP samples in SRP165679:", nrow(sc165), "\n")
print(table(sc165$class))

# match RSE columns (external_id) to classification
ext <- colData(rse)$external_id
if (is.null(ext)) ext <- rse$external_id
keep <- ext %in% sc165$external_id
cat("RSE cols matched:", sum(keep), "\n")
rse2 <- rse[, keep]
cls <- sc165$class[match(ext[keep], sc165$external_id)]
tier <- c(NN=0,PN=1,PP=2)[cls]
cat("matched tiers:\n"); print(table(cls))

# logCPM via edgeR
cnts <- assay(rse2, "counts")
dge <- DGEList(counts=cnts, group=cls)
dge <- calcNormFactors(dge)
# filter low-expressed
keepg <- filterByExpr(dge, group=cls)
dge <- dge[keepg,,keep.lib.sizes=FALSE]
logcpm <- cpm(dge, log=TRUE, prior.count=1)
cat("logCPM dim:", paste(dim(logcpm),collapse=" x "), "\n")

# collapse Ensembl -> symbol (max-expression on duplicates)
gn <- rowData(rse2)$gene_name[keepg]
rownames(logcpm) <- gn
# max-expression collapse
ord <- order(rowMeans(logcpm), decreasing=TRUE)
logcpm <- logcpm[ord,]
logcpm <- logcpm[!duplicated(rownames(logcpm)),]
logcpm <- logcpm[!is.na(rownames(logcpm)) & rownames(logcpm)!="",]
cat("collapsed bulk symbols:", nrow(logcpm), "\n")

saveRDS(list(logcpm=logcpm, tier=tier, cls=cls, samples=ext[keep], srp="SRP165679"),
        "bulk_anchor.rds")
cat("STAT3 in anchor:", "STAT3" %in% rownames(logcpm), "\n")
cat("saved bulk_anchor.rds\n")