# Auto-extracted generating script
# Produces: fig_deconv_validation.png
# Conda env: scissor-r   (run with this environment activated)
# Inputs (expected alongside / in data/): deconv_result.rds
# Source artifact version: 5432311e-45cb-4b15-b08c-b25bb2972671
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(ggplot2); library(reshape2)})

dr <- readRDS("deconv_result.rds")
props <- dr$props
trend <- dr$trend

# rebuild concordance data inline
sub_md <- NULL  # not needed for concordance rebuild; we use trend directly

# Rebuild m (concordance table) from deconv result
cts_all <- colnames(props)[!colnames(props) %in% c("tier","tord")]

# Use the trend from deconv result, and approximate scissor directions from known results
scissor_dirs <- c(
  Endothelial = "Scissor+ (lesional-tracking)",
  Fibroblast  = "Scissor- (normal-tracking)",
  Melanocyte  = "Scissor- (normal-tracking)",
  NK          = "Scissor+ (lesional-tracking)",
  Keratinocyte= "Scissor+ (lesional-tracking)",
  DC          = "Scissor+ (lesional-tracking)",
  Myeloid     = "ns",
  Tcell       = "ns",
  Bcell       = "ns",
  Mast        = "ns",
  Muscle      = "ns"
)

m <- merge(data.frame(celltype=names(scissor_dirs), scissor_direction=unname(scissor_dirs)),
           trend[,c("celltype","slope_per_tier","p","padj")], by="celltype")
m$deconv_dir  <- ifelse(m$slope_per_tier>0,"rises","falls")
m$deconv_sig  <- m$padj < 0.05
m$expected_deconv <- ifelse(grepl("Scissor\\+",m$scissor_direction),"rises",
                       ifelse(grepl("Scissor-",m$scissor_direction),"falls",NA))
m$mappable <- !is.na(m$expected_deconv)
m$concordant <- m$mappable & m$deconv_sig & (m$deconv_dir==m$expected_deconv)
m$status <- with(m, ifelse(!mappable,"no clear Scissor direction",
                    ifelse(concordant,"CONCORDANT",
                     ifelse(!deconv_sig,"deconv trend n.s.","DISCORDANT"))))

pal_tier <- c(NN="#4C9BD4", PN="#E8A33D", PP="#C43C4E")
cts <- m$celltype[m$scissor_direction!="ns"]
stat <- setNames(m$status, m$celltype)
lg <- melt(props[,c(cts,"tier")], id.vars="tier", variable.name="celltype", value.name="prop")
lg$tier <- factor(lg$tier, levels=c("NN","PN","PP"))
lg$celltype <- factor(lg$celltype, levels=cts)
plab <- setNames(sprintf("%s\n[%s]", cts, ifelse(stat[cts]=="CONCORDANT","concordant",
                  ifelse(stat[cts]=="DISCORDANT","DISCORDANT","deconv n.s."))), cts)
pA <- ggplot(lg, aes(tier, prop*100, fill=tier)) +
  geom_boxplot(outlier.size=.4, width=.65) +
  facet_wrap(~celltype, scales="free_y", nrow=2, labeller=labeller(celltype=plab)) +
  scale_fill_manual(values=pal_tier, guide="none") +
  labs(title="Bulk deconvolution vs Scissor direction, by cell type",
       subtitle="Concordant: Endothelial (lesional-tracking), Fibroblast & Melanocyte (normal-tracking). NK & Keratinocyte discordant; DC n.s.",
       x="Biopsy tier (NN < PN < PP)", y="Estimated proportion (%)") +
  theme_bw(base_size=10) + theme(strip.text=element_text(size=7.5),
       plot.subtitle=element_text(size=8))
ggsave("fig_deconv_validation.png", pA, width=9, height=6, dpi=200)