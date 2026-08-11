
suppressMessages({library(recount3)})
# Full human SRA project list
proj <- available_projects(organism="human")
sra  <- proj[proj$file_source=="sra", ]
cat("human SRA projects:", nrow(sra), "\n")

# Fetch per-project metadata and keyword-screen for psoriasis + blood
# recount3 metadata files are small; pull sample_attributes only.
res <- list()
for (i in seq_len(nrow(sra))) {
  p <- sra$project[i]
  ok <- tryCatch({
    # Only the `sra` metadata file carries study_title/sample_attributes.
    # read_metadata() pulls all 5 files per project (~6s); fetching just this one is ~5x faster.
    urls <- locate_url(p, "data_sources/sra", "metadata", "human")
    u <- grep("\\.sra\\.", urls, value=TRUE); if(!length(u)) u <- urls[1]
    md <- read.delim(file_retrieve(u), sep="\t", header=TRUE, quote="", comment.char="")
    # raw read gives UNPREFIXED names (read_metadata() would add the "sra." prefix)
    gv <- function(k) if (k %in% names(md)) paste(head(unique(md[[k]]),50), collapse=" ") else ""
    txt <- tolower(paste(gv("study_title"), gv("study_abstract"), gv("sample_attributes")))
    psor  <- grepl("psoria", txt)
    blood <- grepl("blood|pbmc|mononuclear|neutrophil|leukocyte|monocyte|platelet|serum|plasma", txt)
    if (psor && blood) {
      res[[p]] <- data.frame(project=p, n=nrow(md),
        title=substr(if ("sra.study_title" %in% names(md)) unique(md$sra.study_title)[1] else "", 1, 180),
        has_healthy=grepl("healthy|normal|control", txt),
        stringsAsFactors=FALSE)
      cat("HIT", p, "n=", nrow(md), "\n")
    }
    TRUE
  }, error=function(e) FALSE)
  if (i %% 500 == 0) cat("  ...scanned", i, "/", nrow(sra), "hits:", length(res), "\n")
}
out <- if (length(res)) do.call(rbind, res) else data.frame()
write.csv(out, "recount3_blood_psoriasis_census.csv", row.names=FALSE)
cat("\nDONE. hits:", nrow(out), "\n")
