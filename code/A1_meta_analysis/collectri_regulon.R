# Auto-extracted generating script
# Produces: collectri_regulon.rds
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: 34ad206e-28f3-42b3-9de6-c736f9abc2f3
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)

u <- "https://omnipathdb.org/interactions?resources=CollecTRI&datasets=collectri&fields=sources,references&genesymbols=1&organisms=9606"
raw <- fread(u, sep="\t", header=TRUE)

net <- raw[, .(source = source_genesymbol, target = target_genesymbol,
               mor = fifelse(is_inhibition & !is_stimulation, -1, 1))]
net <- unique(net[source!="" & target!=""])
net <- net[, .(mor = fifelse(sum(mor) < 0, -1, 1)), by=.(source,target)]

saveRDS(net, "collectri_regulon.rds")