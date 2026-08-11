#!/usr/bin/env Rscript
# census_blood_par.R -- parallel recount3 blood-psoriasis census (worker shard)
#
# Serial version ran ~7.8 s/project on the cluster (cold cache; the work is
# network-latency-bound, not CPU-bound), i.e. ~19 h for 8,677 projects.
# Sharding by index across N workers overlaps that latency.
#
# Usage: Rscript census_blood_par.R <worker_id 0..N-1> <N> <outdir>
suppressMessages({library(recount3)})
a  <- commandArgs(TRUE); k <- as.integer(a[1]); N <- as.integer(a[2]); outdir <- a[3]
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# Project list: read ONCE from a shared RDS. Calling available_projects() per
# worker re-downloads ~8 data-source manifests each (~30 s) and hammers the server.
# Worker 0 builds it; the rest wait for it to appear.
shared <- file.path(outdir, "sra_projects.rds")
if (k == 0 && !file.exists(shared)) {
  p <- available_projects(organism="human")
  saveRDS(p[p$file_source=="sra", ], paste0(shared, ".tmp")); file.rename(paste0(shared,".tmp"), shared)
}
for (w in 1:120) { if (file.exists(shared)) break; Sys.sleep(5) }
if (!file.exists(shared)) stop("worker ", k, ": shared project list never appeared")
sra <- readRDS(shared)

# per-worker BiocFileCache for the per-project fetches: concurrent writers to one
# cache corrupt its sqlite index
cache <- file.path(outdir, sprintf("cache_%02d", k))
dir.create(cache, showWarnings=FALSE, recursive=TRUE)
Sys.setenv(RECOUNT3_CACHE = cache, XDG_CACHE_HOME = cache)
idx  <- which(seq_len(nrow(sra)) %% N == k)
cat("worker", k, "of", N, "-> ", length(idx), "projects\n"); flush.console()

res <- list(); done <- 0
for (i in idx) {
  p <- sra$project[i]
  tryCatch({
    urls <- locate_url(p, "data_sources/sra", "metadata", "human")
    u <- grep("\\.sra\\.", urls, value=TRUE); if (!length(u)) u <- urls[1]
    md <- read.delim(file_retrieve(u), sep="\t", header=TRUE, quote="", comment.char="")
    # raw read gives UNPREFIXED column names (read_metadata() adds the "sra." prefix)
    gv <- function(key) if (key %in% names(md))
            paste(head(unique(md[[key]]), 50), collapse=" ") else ""
    txt <- tolower(paste(gv("study_title"), gv("study_abstract"), gv("sample_attributes")))
    if (grepl("psoria", txt) &&
        grepl("blood|pbmc|mononuclear|neutrophil|leukocyte|monocyte|platelet|serum|plasma", txt)) {
      res[[p]] <- data.frame(project=p, n=nrow(md),
                             title=substr(gv("study_title"), 1, 180),
                             has_healthy=grepl("healthy|normal|control", txt),
                             stringsAsFactors=FALSE)
      cat("HIT", p, "n=", nrow(md), "|", substr(gv("study_title"),1,90), "\n"); flush.console()
    }
  }, error=function(e) NULL)
  done <- done + 1
  if (done %% 100 == 0) { cat("  w",k,": ",done,"/",length(idx)," hits:",length(res),"\n",sep=""); flush.console() }
}
out <- if (length(res)) do.call(rbind, res) else
       data.frame(project=character(), n=integer(), title=character(), has_healthy=logical())
write.csv(out, file.path(outdir, sprintf("census_shard_%02d.csv", k)), row.names=FALSE)
cat("WORKER", k, "DONE hits:", nrow(out), "\n")
