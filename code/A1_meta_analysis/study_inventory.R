# Auto-extracted generating script
# Produces: study_inventory.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): human_studies.json, email.json
# Source artifact version: 93b5cdec-4729-4ca7-b482-72f1ee8373f5
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(recount3); library(data.table); library(jsonlite)})

hs <- fromJSON("human_studies.json")
projs <- hs$project

get_meta <- function(srp) {
  tryCatch({
    md <- read_metadata(file_retrieve(locate_url(srp, "data_sources/sra",
             "metadata", organism="human")))
    sa <- if ("sra.sample_attributes" %in% colnames(md)) md$sra.sample_attributes else rep(NA, nrow(md))
    list(srp=srp, n=nrow(md), attrs=sa,
         src=if("sra.library_source"%in%colnames(md)) unique(md$sra.library_source) else NA)
  }, error=function(e) list(srp=srp, n=NA, attrs=NA, src=paste("ERR",conditionMessage(e))))
}
meta <- lapply(projs, get_meta)
names(meta) <- projs

studies <- c("SRP035988","SRP165679","SRP076982","SRP126422","SRP016583")
mds <- lapply(studies, function(srp) {
  read_metadata(file_retrieve(locate_url(srp,"data_sources/sra","metadata",organism="human")))
})
names(mds) <- studies

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

cls_all <- rbindlist(lapply(studies, function(s) classify_samples(s, mds[[s]])))

depth_of <- function(md) {
  cc <- if ("recount_qc.star.number_of_input_reads_both" %in% colnames(md))
          "recount_qc.star.number_of_input_reads_both" else "recount_qc.star.number_of_input_reads"
  suppressWarnings(as.numeric(md[[cc]]))
}
qc <- rbindlist(lapply(names(mds), function(s){
  v <- depth_of(mds[[s]])
  data.frame(srp=s, median_M_reads=round(median(v,na.rm=TRUE)/1e6,1),
             min_M_reads=round(min(v,na.rm=TRUE)/1e6,1),
             pct_below10M=round(100*mean(v<10e6,na.rm=TRUE)))
}))

tab <- dcast(cls_all[class!="DROP"], srp ~ class, fun.aggregate=length, value.var="external_id")
for (cc in c("PP","PN","NN")) if (!cc %in% colnames(tab)) tab[[cc]] <- 0L
setDT(tab); tab[, total := PP+PN+NN]

labs <- c(SRP035988="Li 2014 (case-control)", SRP165679="Tsoi 2019 (AD/PSO cohort)",
          SRP076982="Anatomic plaques", SRP126422="Keratinocyte+biopsy", SRP016583="Early profiling")
inv <- merge(tab, qc, by="srp")
inv[, label := labs[srp]]
inv[, depth_flag := ifelse(median_M_reads < 10, "SHALLOW", "ok")]
inv <- inv[order(-total)]
setcolorder(inv, c("srp","label","PP","PN","NN","total","median_M_reads","min_M_reads","pct_below10M","depth_flag"))

fwrite(inv, "study_inventory.csv")