# Auto-extracted generating script
# Produces: class_enrichment.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gradient_de_genes.rds
# Source artifact version: 0035527f-9bfe-4265-b86b-043351abab4e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(msigdbr)
library(data.table)

g <- readRDS("gradient_de_genes.rds")

H <- as.data.table(msigdbr(species="Homo sapiens", category="H"))
Hsets <- split(H$gene_symbol, H$gs_name)

universe <- unique(g$gene)

enr <- function(genes, topn=8){
  genes <- intersect(genes, universe)
  res <- rbindlist(lapply(names(Hsets), function(s){
    inset <- intersect(Hsets[[s]], universe)
    a <- length(intersect(genes, inset))
    b <- length(genes)-a
    c <- length(inset)-a
    d <- length(universe)-a-b-c
    data.table(set=sub("HALLMARK_","",s), overlap=a, set_size=length(inset),
               p=fisher.test(matrix(c(a,b,c,d),2), alternative="greater")$p.value)
  }))
  res[, FDR := p.adjust(p,"BH")]
  res[order(p)][1:topn]
}

enr_prog <- enr(g[cls=="progressive", gene])
enr_prog[, class := "progressive (early/inflammatory)"]
enr_late <- enr(g[cls=="late_PP", gene])
enr_late[, class := "late_PP (proliferation)"]
enr_all <- rbind(enr_prog, enr_late)

saveRDS(enr_all, "class_enrichment.rds")