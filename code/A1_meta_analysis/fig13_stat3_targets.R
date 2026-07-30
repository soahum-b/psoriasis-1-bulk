# Auto-extracted generating script
# Produces: fig13_stat3_targets.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds
# Source artifact version: f8d8dddb-e69b-4237-8621-500e8f2aa0a8
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(ggplot2)
library(ggrepel)

de_tab <- readRDS("de_results_full.rds")

u <- "https://omnipathdb.org/interactions?resources=CollecTRI&datasets=collectri&fields=sources,references&genesymbols=1&organisms=9606"
raw <- fread(u, sep="\t", header=TRUE)

net <- raw[, .(source = source_genesymbol, target = target_genesymbol,
               mor = fifelse(is_inhibition & !is_stimulation, -1, 1))]
net <- unique(net[source!="" & target!=""])
net <- net[, .(mor = fifelse(sum(mor) < 0, -1, 1)), by=.(source,target)]

NNblue<-"#4C72B0"; PPred<-"#C44E52"

s <- net[source=="STAT3"]
dt <- merge(s, data.table(target=de_tab$gene, t=de_tab$t, logFC=de_tab$logFC), by="target")
dt[, dir := fifelse(mor>0, "activated by STAT3 (+)", "repressed by STAT3 (-)")]
cat("STAT3 targets matched:", nrow(dt), "\n")
cat("mean t of activated targets:", round(mean(dt[mor>0]$t),2),
    "| repressed:", round(mean(dt[mor<0]$t),2), "\n")

lab <- rbind(dt[mor>0][order(-t)][1:10], dt[mor<0][order(t)][1:6])

set.seed(1)
p13 <- ggplot(dt, aes(factor(mor), t, colour=dir)) +
  geom_hline(yintercept=0, colour="grey60", linewidth=0.3) +
  geom_jitter(width=0.18, height=0, alpha=0.55, size=1.6) +
  geom_boxplot(width=0.35, outlier.shape=NA, fill=NA, colour="grey25", linewidth=0.5) +
  ggrepel::geom_text_repel(data=lab, aes(label=target), size=2.9, max.overlaps=30,
                           box.padding=0.4, min.segment.length=0, seed=1, show.legend=FALSE) +
  scale_colour_manual(values=c("activated by STAT3 (+)"=PPred, "repressed by STAT3 (-)"=NNblue), name=NULL) +
  scale_x_discrete(labels=c("-1"="repressed\ntargets (n=39)", "1"="activated\ntargets (n=396)")) +
  labs(title="STAT3 is an ACTIVE regulator - its targets move as its regulon predicts",
       subtitle="Each point = one STAT3 target gene. Activated targets shift up, repressed shift down: the signature of an active TF.",
       x=NULL, y="Observed change in lesional skin (limma t-statistic)") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold", size=12),
        legend.position="top")

ggsave("fig13_stat3_targets.png", p13, width=8.5, height=6.3, dpi=150)