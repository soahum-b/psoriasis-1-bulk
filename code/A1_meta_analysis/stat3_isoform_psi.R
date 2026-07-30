# Auto-extracted generating script
# Produces: stat3_isoform_psi.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds
# Source artifact version: 118993f9-c40d-4a88-a46c-127c7e560491
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

# --- fetch junction-level RSE ---
ap <- available_projects()
proj_info <- subset(ap, project == "SRP035988" & project_type == "data_sources")
rse_jxn <- create_rse(proj_info, type = "jxn")

# --- re-derive NN/PP groups from sample_attributes ---
attr_raw <- colData(rse_jxn)$sra.sample_attributes
grp <- ifelse(grepl("normal", attr_raw, ignore.case=TRUE), "NN",
       ifelse(grepl("Psoriasis", attr_raw, ignore.case=TRUE), "PP", NA))

# --- STAT3 locus, GRCh38 minus strand ---
gr <- rowRanges(rse_jxn)
stat3_win <- GRanges("chr17", IRanges(42313000, 42389000))
hits <- which(overlapsAny(gr, stat3_win) & as.character(strand(gr)) == "-")

# alpha-specific junction: start=42316902, end=42317181
# beta-specific junction:  start=42316852, end=42317181
i_alpha <- which(start(gr)==42316902 & end(gr)==42317181)
i_beta  <- which(start(gr)==42316852 & end(gr)==42317181)

ca <- as.numeric(assay(rse_jxn)[i_alpha, ])
cb <- as.numeric(assay(rse_jxn)[i_beta, ])

dt <- data.table(sample=colnames(rse_jxn), grp=grp, a=ca, b=cb)
dt[, depth := a + b]
dt[, psi_beta := b / depth]

saveRDS(dt, "stat3_isoform_psi.rds")