# Auto-extracted generating script
# Produces: sample_manifest.csv
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: a8220646-c817-4350-bfe9-f46609d614d3
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(readr)
library(stringr)
library(dplyr)

raw_dir <- "scissor_repo/data/raw"
files <- sort(list.files(raw_dir, pattern = "\\.csv\\.gz$", full.names = FALSE))

rows <- lapply(files, function(f) {
  base <- sub("\\.csv\\.gz$", "", f)
  parts <- strsplit(base, "_", fixed = TRUE)[[1]]
  gsm <- parts[1]
  samp <- paste(parts[-1], collapse = "_")
  pref <- toupper(strsplit(samp, "-")[[1]][1])
  cond <- c("NS" = "NN", "PN" = "PN", "PP" = "PP")[pref]
  tier <- c("NS" = 0L, "PN" = 1L, "PP" = 2L)[pref]
  data.frame(gsm = gsm, sample = samp, condition = unname(cond), tier = unname(tier), file = f, stringsAsFactors = FALSE)
})

manifest <- do.call(rbind, rows)

write.csv(manifest, "scissor_repo/data/sample_manifest.csv", row.names = FALSE)