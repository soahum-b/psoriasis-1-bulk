# Auto-extracted generating script
# Produces: stat3_isoform_eftud2.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): dge_filt_norm.rds
# Source artifact version: a4d91fd8-9345-4f42-8227-934edd8dda50
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(data.table)
  library(edgeR)
  library(httr)
  library(jsonlite)
})

# Load filtered+normalized gene-level DGEList checkpoint
dge <- readRDS("dge_filt_norm.rds")
logcpm <- cpm(dge, log=TRUE)

# Fetch junction-level RSE for SRP035988
ap <- available_projects()
proj_info <- subset(ap, project == "SRP035988" & project_type == "data_sources")
rse_jxn <- create_rse(proj_info, type = "jxn")

# Derive NN/PP groups from sample_attributes
attr_raw <- colData(rse_jxn)$sra.sample_attributes
grp <- ifelse(grepl("normal", attr_raw, ignore.case=TRUE), "NN",
       ifelse(grepl("Psoriasis", attr_raw, ignore.case=TRUE), "PP", NA))

# Get junctions in STAT3 locus (chr17, minus strand)
gr <- rowRanges(rse_jxn)
stat3_win <- GRanges("chr17", IRanges(42313000, 42389000))
hits <- which(overlapsAny(gr, stat3_win) & as.character(strand(gr)) == "-")

# Identify alpha- and beta-specific junctions
# alpha junction: 42316902-42317181
# beta junction:  42316852-42317181
i_alpha <- which(start(gr)==42316902 & end(gr)==42317181)
i_beta  <- which(start(gr)==42316852 & end(gr)==42317181)

ca <- as.numeric(assay(rse_jxn)[i_alpha, ])
cb <- as.numeric(assay(rse_jxn)[i_beta, ])

# Build per-sample PSI table
dt <- data.table(sample=colnames(rse_jxn), grp=grp, a=ca, b=cb)
dt[, depth := a+b]
dt[, psi_beta := b/depth]

# Add EFTUD2 and STAT3 logCPM to per-sample table
setkey(dt, sample)
dt[, EFTUD2 := logcpm["EFTUD2", sample]]
dt[, STAT3  := logcpm["STAT3",  sample]]

saveRDS(dt, "stat3_isoform_eftud2.rds")