# Auto-extracted generating script
# Produces: sample_classification.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): email.json
# Source artifact version: 5e70760d-efdc-41f1-81d9-b3533488838e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(recount3)
library(data.table)
library(jsonlite)

studies <- c("SRP035988", "SRP165679", "SRP076982", "SRP126422", "SRP016583")

load_md <- function(srp) read_metadata(file_retrieve(locate_url(srp, "data_sources/sra", "metadata", organism="human")))

classify_samples <- function(srp, md) {
  a  <- md$sra.sample_attributes
  st <- md$sra.sample_title
  ext <- md$external_id
  cls <- rep("DROP", nrow(md))
  if (srp=="SRP035988") {
    cls[grepl("lesional psoriatic", a)] <- "PP"
    cls[grepl("normal skin", a)]        <- "NN"
  } else if (srp=="SRP165679") {
    dis  <- ifelse(grepl("^AD",st),"AD", ifelse(grepl("^PSO",st),"PSO","HC"))
    skin <- ifelse(grepl("non-lesional",st),"PN", ifelse(grepl("lesional",st),"PP","NN"))
    cls[dis=="PSO" & skin=="PP"] <- "PP"
    cls[dis=="PSO" & skin=="PN"] <- "PN"
    cls[dis=="HC"]               <- "NN"
  } else if (srp=="SRP076982") {
    cls[grepl("tissue type;;psoriasis", a)]   <- "PP"
    cls[grepl("tissue type;;uninvolved", a)]  <- "PN"
  } else if (srp=="SRP126422") {
    skinbx <- grepl("skin biopsy", a, ignore.case=TRUE)
    cls[skinbx & grepl("\\(PP\\)", a)] <- "PP"
    cls[skinbx & grepl("\\(NN\\)", a)] <- "NN"
    cls[skinbx & grepl("\\(PN\\)", a)] <- "PN"
  } else if (srp=="SRP016583") {
    cls[grepl("group;;lesional", a, ignore.case=TRUE)]     <- "PP"
    cls[grepl("group;;non-lesional", a, ignore.case=TRUE)] <- "PN"
  }
  data.frame(external_id=ext, srp=srp, class=cls, stringsAsFactors=FALSE)
}

mds <- lapply(studies, load_md); names(mds) <- studies
cls_all <- rbindlist(lapply(studies, function(s) classify_samples(s, mds[[s]])))

saveRDS(cls_all, "sample_classification.rds")