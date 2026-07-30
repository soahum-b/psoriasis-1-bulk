# Auto-extracted generating script
# Produces: fig12_tf_activity_rank.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds
# Source artifact version: 563aa009-56e0-4fe2-ab41-55ca4a23704b
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(decoupleR)
library(ggplot2)
library(ggrepel)

de_tab <- readRDS("de_results_full.rds")

mat <- matrix(de_tab$t, ncol=1, dimnames=list(de_tab$gene, "PPvsNN"))
mat <- mat[!is.na(mat[,1]), , drop=FALSE]

u <- "https://omnipathdb.org/interactions?resources=CollecTRI&datasets=collectri&fields=sources,references&genesymbols=1&organisms=9606"
raw <- fread(u, sep="\t", header=TRUE)

net <- raw[, .(source = source_genesymbol, target = target_genesymbol,
               mor = fifelse(is_inhibition & !is_stimulation, -1, 1))]
net <- unique(net[source!="" & target!=""])
net <- net[, .(mor = fifelse(sum(mor) < 0, -1, 1)), by=.(source,target)]

set.seed(42)
tf_res <- run_ulm(mat = mat, network = net, .source="source", .target="target",
                  .mor="mor", minsize = 5)
setDT(tf_res)
tf_res <- tf_res[order(-score)]
tf_res[, padj := p.adjust(p_value, "BH")]
tf_res[, rank := .I]

axis_tf <- c("STAT3","STAT1","RELA","RELB","REL","NFKB1","NFKB2","RORC","RORA","IRF1","MYC","E2F1","JUNB","FOS","JUN")
tf_res[, axis := source %in% axis_tf]

NNblue<-"#4C72B0"; PPred<-"#C44E52"; grey<-"#B0B0B0"

plotdt <- copy(tf_res)[!grepl("_", source)]
plotdt[, rank2 := frank(-score, ties.method="first")]

labdt <- plotdt[source %in% c("STAT3","STAT1","RELA","NFKB1","RELB","MYC","E2F1","IRF1","RORC","JUNB") | rank2<=3]

p12 <- ggplot(plotdt, aes(rank2, score)) +
  geom_hline(yintercept=0, colour="grey70", linewidth=0.3) +
  geom_point(aes(colour=axis, size=axis, alpha=axis)) +
  scale_colour_manual(values=c(`FALSE`=grey, `TRUE`=PPred), guide="none") +
  scale_size_manual(values=c(`FALSE`=1.1, `TRUE`=2.6), guide="none") +
  scale_alpha_manual(values=c(`FALSE`=0.45, `TRUE`=1), guide="none") +
  ggrepel::geom_text_repel(data=labdt, aes(label=source), size=3.2, fontface="bold",
                           max.overlaps=30, box.padding=0.5, min.segment.length=0,
                           colour="grey15", seed=1) +
  annotate("text", x=nrow(plotdt)*0.55, y=9.2,
           label="Inflammatory hubs (red) dominate the\nmost-active transcription factors", hjust=0, size=3.4, colour=PPred) +
  labs(title="Transcription-factor activity in psoriasis - inferred from target behaviour (decoupleR ULM, CollecTRI)",
       subtitle="732 TFs ranked by activity score (t-value of regulon slope). Red = STAT3 / NF-kB / IFN axis.",
       x="TF activity rank (most active -> least)", y="Activity score (ULM t-value)") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold", size=11))

ggsave("fig12_tf_activity_rank.png", p12, width=10, height=6.2, dpi=150)