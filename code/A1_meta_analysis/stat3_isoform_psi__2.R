# Auto-extracted generating script
# Produces: stat3_isoform_psi.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: c1a2e236-7429-4be0-8c5a-0ba4612fdaf9
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(data.table)
  library(httr)
  library(jsonlite)
})

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")

# fetch junction-level RSE
ap <- available_projects()
proj_info <- subset(ap, project == "SRP035988" & project_type == "data_sources")
rse_jxn <- create_rse(proj_info, type = "jxn")

# derive NN/PP groups
attr_raw <- colData(rse_jxn)$sra.sample_attributes
grp <- ifelse(grepl("normal", attr_raw, ignore.case=TRUE), "NN",
       ifelse(grepl("Psoriasis", attr_raw, ignore.case=TRUE), "PP", NA))

# STAT3 locus junctions
gr <- rowRanges(rse_jxn)
stat3_win <- GRanges("chr17", IRanges(42313000, 42389000))
hits <- which(overlapsAny(gr, stat3_win) & as.character(strand(gr)) == "-")

# identify alpha- and beta-specific junction indices
i_alpha <- which(start(gr)==42316902 & end(gr)==42317181)
i_beta  <- which(start(gr)==42316852 & end(gr)==42317181)

ca <- as.numeric(assay(rse_jxn)[i_alpha, ])
cb <- as.numeric(assay(rse_jxn)[i_beta, ])

dt <- data.table(sample=colnames(rse_jxn), grp=grp, a=ca, b=cb)
dt[, depth := a + b]
dt[, psi_beta := b / depth]

fwrite(dt, "stat3_isoform_psi.csv")