# Auto-extracted generating script
# Produces: tf_activity_collectri.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds
# Source artifact version: 1c168e81-a2ac-4f7e-9812-fc12e4a96ecc
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(decoupleR)
library(data.table)

de_tab <- readRDS("de_results_full.rds")

u <- "https://omnipathdb.org/interactions?resources=CollecTRI&datasets=collectri&fields=sources,references&genesymbols=1&organisms=9606"
raw <- fread(u, sep="\t", header=TRUE)

net <- raw[, .(source = source_genesymbol, target = target_genesymbol,
               mor = fifelse(is_inhibition & !is_stimulation, -1, 1))]
net <- unique(net[source!="" & target!=""])
net <- net[, .(mor = fifelse(sum(mor) < 0, -1, 1)), by=.(source,target)]

mat <- matrix(de_tab$t, ncol=1, dimnames=list(de_tab$gene, "PPvsNN"))
mat <- mat[!is.na(mat[,1]), , drop=FALSE]

set.seed(42)
tf_res <- run_ulm(mat = mat, network = net, .source="source", .target="target",
                  .mor="mor", minsize = 5)
setDT(tf_res)
tf_res <- tf_res[order(-score)]
tf_res[, padj := p.adjust(p_value, "BH")]
tf_res[, rank := .I]

axis_tf <- c("STAT3","STAT1","RELA","RELB","REL","NFKB1","NFKB2","RORC","RORA","IRF1","MYC","E2F1","JUNB","FOS","JUN")
tf_res[, axis := source %in% axis_tf]

saveRDS(tf_res, "tf_activity_collectri.rds")