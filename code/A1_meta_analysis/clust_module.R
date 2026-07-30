# Auto-extracted generating script
# Produces: clust_module.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): meta_de_results.rds
# Source artifact version: f47bd7fd-da30-44ab-a8a9-df3d23eb5ede
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(msigdbr)

meta <- readRDS("meta_de_results.rds")
mp <- as.data.table(meta$PPvsNN)
bg <- mp[FDR<0.05 & abs(logFC)>1, gene]

lines <- readLines("clust_out3/Clusters_Objects.tsv")
tab <- fread("clust_out3/Clusters_Objects.tsv", skip=1, header=TRUE)
mod <- tab[[1]]; mod <- mod[mod!="" & !is.na(mod)]
cat("module size:", length(mod), "\n")

H <- msigdbr(species="Homo sapiens", category="H")
Hl <- split(H$gene_symbol, H$gs_name)
enr <- rbindlist(lapply(names(Hl), function(nm){
  s <- intersect(Hl[[nm]], bg); k <- intersect(mod, s)
  if (length(k)<2) return(NULL)
  p <- phyper(length(k)-1, length(s), length(bg)-length(s), length(mod), lower.tail=FALSE)
  data.table(set=sub("HALLMARK_","",nm), overlap=length(k), setsize=length(s), p=p,
             genes=paste(head(k,8),collapse=","))
}))
enr <- enr[order(p)][1:8]

modlfc <- mp[gene %in% mod, .(gene, logFC)]

saveRDS(list(module=mod, enrichment=enr, logfc=modlfc), "clust_module.rds")
fwrite(data.table(gene=mod, logFC=mp$logFC[match(mod,mp$gene)]), "clust_module_genes.csv")