# Auto-extracted generating script
# Produces: coexpr_lnc_target_unadjusted.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): lncRNA_significant_ranked.csv, rebuild_logcpm.R
# Source artifact version: 4d7dc584-7d8b-474a-bfd9-1e81f198eb71
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(recount3); library(edgeR); library(limma); library(data.table)})
options(recount3_url = "http://duffel.rail.bio/recount3")
studies <- c("SRP035988","SRP165679","SRP126422","SRP065812")

classify <- function(srp, md){
  a <- md$sra.sample_attributes; st <- md$sra.sample_title; ext <- md$external_id
  cls <- rep("DROP", nrow(md))
  if(srp=="SRP035988"){ cls[grepl("lesional psoriatic",a)]<-"PP"; cls[grepl("normal skin",a)]<-"NN" }
  else if(srp=="SRP165679"){
    dis<-ifelse(grepl("^PSO",st),"PSO",ifelse(grepl("^AD",st),"AD","HC"))
    skin<-ifelse(grepl("non-lesional",st),"PN",ifelse(grepl("lesional",st),"PP","NN"))
    cls[dis=="PSO"&skin=="PP"]<-"PP"; cls[dis=="HC"]<-"NN"
  } else if(srp=="SRP126422"){
    sb<-grepl("skin biopsy",a,ignore.case=TRUE)
    cls[sb&grepl("\\(PP\\)",a)]<-"PP"; cls[sb&grepl("\\(NN\\)",a)]<-"NN"
  } else if(srp=="SRP065812"){
    cls[grepl("before adalimumab",a,ignore.case=TRUE) & grepl("lesional",a,ignore.case=TRUE)]<-"PP"
    cls[grepl("healthy control",a,ignore.case=TRUE) | grepl("normal skin",a,ignore.case=TRUE)]<-"NN"
    cls[grepl("after|post",a,ignore.case=TRUE)]<-"DROP"
  }
  data.table(external_id=ext, srp=srp, class=cls)
}
load_md <- function(srp) read_metadata(file_retrieve(locate_url(srp,"data_sources/sra","metadata",organism="human")))
get_rse <- function(srp){
  rse <- create_rse(subset(available_projects(organism="human"),
                    project==srp & project_type=="data_sources"), type="gene")
  assay(rse,"counts") <- transform_counts(rse); rse
}

logcpm_one <- function(srp){
  md <- load_md(srp); cl <- classify(srp, md); rse <- get_rse(srp)
  cvec <- cl$class[match(colnames(rse), cl$external_id)]
  keep_s <- which(cvec %in% c("NN","PP"))
  grp <- factor(cvec[keep_s], levels=c("NN","PP"))
  dge <- DGEList(counts=assay(rse,"counts")[,keep_s],
                 genes=data.frame(gene_name=rowData(rse)$gene_name))
  keep <- filterByExpr(dge, group=grp); dge<-dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge,"TMM")
  lcpm <- cpm(dge, log=TRUE, prior.count=1)
  rownames(lcpm) <- dge$genes$gene_name
  o <- order(-rowMeans(lcpm)); lcpm <- lcpm[o,]; lcpm <- lcpm[!duplicated(rownames(lcpm)),]
  colnames(lcpm) <- colnames(dge)
  saveRDS(list(logcpm=lcpm, group=grp, study=srp), sprintf("logcpm_%s.rds", srp))
  data.table(study=srp, n_samples=ncol(lcpm), n_genes=nrow(lcpm),
             n_NN=sum(grp=="NN"), n_PP=sum(grp=="PP"))
}

set.seed(1)
summ <- rbindlist(lapply(studies, function(s){ cat("== ", s, " ==\n"); logcpm_one(s) }))
fwrite(summ, "logcpm_rebuild_summary.csv")
print(summ)
cat("DONE\n")

# target panels (from IL1B_STAT3 section)
panels <- list(
  STAT3_JAK_STAT   = c("STAT3","STAT1","JAK1","JAK2","JAK3","TYK2","SOCS3","IL6","IL6R","IL6ST","OSMR"),
  IL1_inflammasome = c("IL1B","IL1A","IL1RN","IL18","CASP1","CASP5","PYCARD","AIM2","NLRP3","GSDMD","IL1R1","IL1R2","IL1RAP"),
  IL36             = c("IL36A","IL36G","IL36B","IL36RN","IL1F10"),
  Th17_output      = c("IL17A","RORC","CCL20","S100A7","S100A8","S100A9","DEFB4A","LCN2")
)
target2panel <- unlist(lapply(names(panels), function(p) setNames(rep(p,length(panels[[p]])), panels[[p]])))
targets <- unique(unlist(panels))

# significant lncRNAs to test
rk <- fread("lncRNA_significant_ranked.csv")
lnc_genes <- unique(rk$gene)

# load logCPM, residualize on group, keep only rows we need
load_resid <- function(srp){
  o <- readRDS(sprintf("logcpm_%s.rds", srp))
  m <- o$logcpm; grp <- o$group
  X <- model.matrix(~grp)
  fitted <- t(X %*% solve(t(X)%*%X) %*% t(X) %*% t(m))
  resid <- m - fitted
  list(resid=resid, raw=m, grp=grp, genes=rownames(m), n=ncol(m))
}
D <- lapply(setNames(studies,studies), load_resid)

# per study: correlation of each lnc vs each target, on residuals (adj) and raw (unadj)
cor_study <- function(d, adj=TRUE){
  M <- if(adj) d$resid else d$raw
  g <- d$genes
  lnc_present <- intersect(lnc_genes, g)
  tgt_present <- intersect(targets, g)
  if(length(lnc_present)==0 || length(tgt_present)==0) return(NULL)
  cl <- M[lnc_present,,drop=FALSE]; ct <- M[tgt_present,,drop=FALSE]
  cc <- cor(t(cl), t(ct))
  dt <- as.data.table(as.table(cc))
  setnames(dt, c("lncRNA","target","r"))
  dt[, n := d$n]
  dt
}

# Fisher z meta across studies
zmeta <- function(rs, ns){
  ok <- is.finite(rs) & ns>3 & abs(rs)<1
  if(sum(ok)==0) return(list(r=NA_real_,p=NA_real_,k=0L,zsum=NA_real_))
  z <- atanh(rs[ok]); w <- ns[ok]-3
  zbar <- sum(w*z)/sum(w); se <- sqrt(1/sum(w))
  zstat <- zbar/se; p <- 2*pnorm(-abs(zstat))
  list(r=tanh(zbar), p=p, k=sum(ok), zstat=zstat)
}

build <- function(adj){
  per <- rbindlist(lapply(studies, function(s) cbind(study=s, cor_study(D[[s]], adj=adj))))
  meta <- per[, {
    zz <- zmeta(r, n)
    .(meta_r=zz$r, meta_p=zz$p, k_studies=zz$k, zstat=zz$zstat,
      r_per_study=paste(sprintf("%s:%.2f", study, r), collapse=";"))
  }, by=.(lncRNA,target)]
  meta[, panel := target2panel[target]]
  meta
}

cat("Building unadjusted (pooled) co-expression...\n")
meta_un <- build(FALSE)
meta_un[, meta_fdr := p.adjust(meta_p, "BH")]
fwrite(meta_un, "coexpr_lnc_target_unadjusted.csv")

cat("\nun rows:", nrow(meta_un), " lncRNAs:", uniqueN(meta_un$lncRNA), " targets:", uniqueN(meta_un$target), "\n")
cat("DONE\n")