# Auto-extracted generating script
# Produces: fig11_unbiased_top_pathways.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): gsea_ora_merged.rds
# Source artifact version: 2ac7281f-efe7-48f1-8af5-ca6b5ab11a01
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(ggplot2); library(data.table)})

M <- readRDS("gsea_ora_merged.rds")

NNblue<-"#4C72B0"; PPred<-"#C44E52"; grey<-"#B0B0B0"

Mdt <- as.data.table(M)
top <- Mdt[gsea_sig==TRUE][order(-gsea_nes)][1:30]

classify <- function(p){
  pl <- toupper(p)
  if (grepl("STAT3|JAK|IL6|IL_17|INTERLEUKIN_17|IL17|SOCS|NF_KB|NFKB|TNFA", pl)) return("STAT3 / IL-17 / NF-kB axis")
  if (grepl("INTERFERON|TYPE_I_INTERFERON|ANTIVIRAL|ISG", pl)) return("Interferon response")
  if (grepl("CELL_CYCLE|MITOTIC|E2F|G2M|G2_M|MYC|DNA_REPLICATION|CHROMATID|CHROMOSOME|SISTER|SYNTHESIS_OF_DNA|CHECKPOINT|SPINDLE|SEPARATION", pl)) return("Proliferation / cell cycle")
  if (grepl("KERATIN|CORNIF|EPIDERM|SKIN", pl)) return("Keratinocyte / epidermis")
  "Other inflammatory / immune"
}
top[, theme := sapply(pathway, classify)]
top[, short := tolower(gsub("_"," ", gsub("^(HALLMARK|REACTOME|KEGG|GOBP)_","",pathway)))]
top[, short := ifelse(nchar(short)>46, paste0(substr(short,1,44),".."), short)]
top[, short := factor(short, levels=rev(short))]

pal <- c("STAT3 / IL-17 / NF-kB axis"=PPred,
         "Proliferation / cell cycle"="#8172B3",
         "Interferon response"="#CCB974",
         "Keratinocyte / epidermis"="#55A868",
         "Other inflammatory / immune"=grey)

p <- ggplot(top, aes(gsea_nes, short, colour=theme)) +
  geom_segment(aes(x=0, xend=gsea_nes, yend=short), linewidth=0.5, colour="grey85") +
  geom_point(aes(size=gsea_size)) +
  scale_colour_manual(values=pal, name="Theme (labeled post-hoc)") +
  scale_size_continuous(name="set size", range=c(2,7)) +
  labs(title="The 30 strongest pathways in psoriasis - ranked by the data, not chosen by hand",
       subtitle="Every significant up-set (padj<0.05) across all 6,505 tested; top 30 by NES. Nothing here was manually selected.",
       x="Normalized Enrichment Score (NES)", y=NULL) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold", size=12),
        legend.position="right", axis.text.y=element_text(size=8.5))

ggsave("fig11_unbiased_top_pathways.png", p, width=11, height=7.5, dpi=150)