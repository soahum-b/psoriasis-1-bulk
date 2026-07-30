# Auto-extracted generating script
# Produces: psi_meta_result.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 8db01554-e157-48f1-878f-94270475f18a
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(data.table)
  library(metafor)
})

# Sample classification
cls_all_raw <- readRDS("/tmp/cls_all.rds")
setDT(cls_all_raw)

ap <- available_projects(organism = "human")

acc <- 42317181L; a_don <- 42316902L; b_don <- 42316852L; tol <- 2L

psi_for_study <- function(srp) {
  proj <- subset(ap, project == srp & project_type == "data_sources")
  rj <- create_rse(proj, type = "jxn")
  gr <- rowRanges(rj)
  onchr <- as.character(seqnames(gr)) %in% c("chr17", "17")
  is_a <- onchr & abs(start(gr) - a_don) <= tol & abs(end(gr) - acc) <= tol
  is_b <- onchr & abs(start(gr) - b_don) <= tol & abs(end(gr) - acc) <= tol
  if (sum(is_a) < 1 || sum(is_b) < 1) return(data.table(srp = srp, note = "junction not found",
                                                          n_a = sum(is_a), n_b = sum(is_b)))
  a_cnt <- colSums(matrix(assay(rj, "counts")[is_a, , drop = FALSE], ncol = ncol(rj)))
  b_cnt <- colSums(matrix(assay(rj, "counts")[is_b, , drop = FALSE], ncol = ncol(rj)))
  cmap <- cls_all_raw[cls_all_raw$srp == srp]
  cl <- cmap$class[match(colnames(rj), cmap$external_id)]
  dt <- data.table(srp = srp, sample = colnames(rj), class = cl,
                   a = a_cnt, b = b_cnt, depth = a_cnt + b_cnt)
  dt[, psi_beta := b / depth]
  dt[class %in% c("PP", "PN", "NN")]
}

studies <- c("SRP035988", "SRP165679", "SRP076982", "SRP126422", "SRP016583")
psi_list <- list()
for (s in studies) {
  message("jxn ", s)
  psi_list[[s]] <- tryCatch(psi_for_study(s), error = function(e) data.table(srp = s, note = conditionMessage(e)))
}
psi_all <- rbindlist(psi_list, fill = TRUE)
saveRDS(psi_all, "psi_beta_allstudies.rds")

MIN_DEPTH <- 20
psi_f <- psi_all[!is.na(psi_beta) & depth >= MIN_DEPTH]

pp_nn_effect <- function(d) {
  pp <- d[class == "PP", psi_beta]; nn <- d[class == "NN", psi_beta]
  if (length(pp) < 3 || length(nn) < 3) return(NULL)
  m <- mean(pp) - mean(nn)
  se <- sqrt(var(pp) / length(pp) + var(nn) / length(nn))
  data.table(delta = m, se = se, n_pp = length(pp), n_nn = length(nn),
             mean_pp = mean(pp), mean_nn = mean(nn),
             w_p = wilcox.test(pp, nn)$p.value)
}
eff <- psi_f[, pp_nn_effect(.SD), by = srp]

if (nrow(eff) >= 2) {
  res <- rma(yi = delta, sei = se, data = eff, method = "DL")
  saveRDS(list(eff = eff, res = res, min_depth = MIN_DEPTH), "psi_meta_result.rds")
} else {
  cat("\nOnly", nrow(eff), "study passes depth filter for PP-vs-NN PSI pooling.\n")
}