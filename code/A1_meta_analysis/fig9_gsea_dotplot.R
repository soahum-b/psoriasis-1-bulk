# Auto-extracted generating script
# Produces: fig9_gsea_dotplot.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): de_results_full.rds, gsea_ranks.rds
# Source artifact version: 93ac3411-16d6-4ff7-ab0b-7bb9d3f1f74e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(data.table)
})

ranks <- readRDS("gsea_ranks.rds")
de    <- readRDS("de_results_full.rds")

set.seed(42)
build_list <- function(coll, sub=NULL) {
  df <- if (is.null(sub)) msigdbr(species="human", collection=coll)
        else msigdbr(species="human", collection=coll, subcollection=sub)
  split(df$gene_symbol, df$gs_name)
}
run_gsea <- function(lst) {
  r <- fgsea(pathways=lst, stats=ranks, minSize=10, maxSize=500, eps=0)
  r[order(r$padj), ]
}

H_list     <- build_list("H")
React_list <- build_list("C2", "CP:REACTOME")
GOBP_list  <- build_list("C5", "GO:BP")
KEGG_list  <- build_list("C2", "CP:KEGG_LEGACY")

fg_H     <- run_gsea(H_list)
fg_React <- run_gsea(React_list)
fg_GOBP  <- run_gsea(GOBP_list)
fg_KEGG  <- run_gsea(KEGG_list)

mk <- function(dt, coll) { d<-as.data.table(dt)[,.(pathway,NES,padj,size)]; d[,collection:=coll]; d }
allp2 <- rbindlist(list(mk(fg_H,"Hallmark"), mk(fg_React,"Reactome"),
                        mk(fg_GOBP,"GO:BP"), mk(fg_KEGG,"KEGG")))

pick <- c(
 "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
 "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
 "HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1",
 "REACTOME_INTERLEUKIN_17_SIGNALING","REACTOME_STAT3_NUCLEAR_EVENTS_DOWNSTREAM_OF_ALK_SIGNALING_AND_PTK6",
 "REACTOME_GENE_AND_PROTEIN_EXPRESSION_BY_JAK_STAT_SIGNALING_AFTER_INTERLEUKIN_12_STIMULATION",
 "KEGG_JAK_STAT_SIGNALING_PATHWAY",
 "GOBP_KERATINIZATION","GOBP_KERATINOCYTE_DIFFERENTIATION",
 "GOBP_RESPONSE_TO_INTERLEUKIN_17","GOBP_INTERLEUKIN_17_PRODUCTION")
dp <- allp2[pathway %in% pick][order(NES)]

pretty <- function(x){
  x <- gsub("^(HALLMARK|REACTOME|KEGG|GOBP)_","",x)
  x <- gsub("_"," ",x); x <- tolower(x)
  x <- gsub("il6","IL6",x); x<-gsub("stat3","STAT3",x); x<-gsub("jak","JAK",x)
  x <- gsub("nfkb","NF-kB",x); x<-gsub("tnfa","TNFa",x); x<-gsub("il 17","IL-17",x)
  x <- gsub("interleukin 17","IL-17",x); x<-gsub("e2f","E2F",x); x<-gsub("g2m","G2M",x); x<-gsub("myc","MYC",x)
  substr(x,1,44)
}
dp[, lab := pretty(pathway)]
dp[, lab := factor(lab, levels=lab)]
dp[, axis := grepl("STAT|IL6|IL-17|JAK|NF-kB|KERATIN", lab, ignore.case=TRUE)]

PPred<-"#C44E52"
p2 <- ggplot(dp, aes(NES, lab)) +
  geom_segment(aes(x=0, xend=NES, y=lab, yend=lab), colour="grey80", linewidth=0.4) +
  geom_point(aes(size=size, colour=-log10(padj))) +
  scale_colour_gradient(low="#F4C4C4", high=PPred, name=expression(-log[10]~adj.P)) +
  scale_size_continuous(name="set size", range=c(2.5,8)) +
  labs(title="GSEA pathway landscape - psoriasis (PP vs NN)",
       subtitle="Positive NES = coordinately up in lesional skin. All shown sets padj < 0.05.",
       x="Normalized Enrichment Score (NES)", y=NULL) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
        axis.text.y=element_text(size=10))

ggsave("fig9_gsea_dotplot.png", p2, width=9.5, height=6, dpi=150)