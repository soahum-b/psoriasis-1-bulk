# Auto-extracted generating script
# Produces: deconv_scissor_concordance.csv
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): deconv_result.rds, reference_subset20k_umap.rds
# Source artifact version: 6058a7de-9ac9-4910-9e57-df11df29876f
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

.libPaths(c(normalizePath("scissor_repo/.rlib"), .libPaths()))
suppressMessages({library(Seurat)})
sub <- readRDS("reference_subset20k_umap.rds")
dr  <- readRDS("deconv_result.rds"); trend <- dr$trend
md <- sub@meta.data
cts <- levels(factor(md$celltype))
tab <- data.frame()
for (ct in cts){
  isct <- md$celltype==ct
  ft <- fisher.test(table(isct, md$scissor=="Scissor+"))
  posOR <- ft$estimate; posP <- ft$p.value
  side <- "ns"
  if (posP<0.05 && posOR>1) side <- "Scissor+ (lesional-tracking)"
  if (posP<0.05 && posOR<1) side <- "Scissor- (normal-tracking)"
  tab <- rbind(tab, data.frame(celltype=ct, pos_enrich_OR=round(posOR,3),
                               pos_enrich_p=signif(posP,2), scissor_direction=side))
}
m <- merge(tab, trend[,c("celltype","slope_per_tier","p","padj")], by="celltype")
m$deconv_dir  <- ifelse(m$slope_per_tier>0,"rises","falls")
m$deconv_sig  <- m$padj < 0.05
m$expected_deconv <- ifelse(grepl("Scissor\\+",m$scissor_direction),"rises",
                       ifelse(grepl("Scissor-",m$scissor_direction),"falls",NA))
m$mappable <- !is.na(m$expected_deconv)
m$concordant <- m$mappable & m$deconv_sig & (m$deconv_dir==m$expected_deconv)
m$status <- with(m, ifelse(!mappable,"no clear Scissor direction",
                    ifelse(concordant,"CONCORDANT",
                     ifelse(!deconv_sig,"deconv trend n.s.","DISCORDANT"))))
m <- m[order(m$pos_enrich_p),]
write.csv(m, "deconv_scissor_concordance.csv", row.names=FALSE)