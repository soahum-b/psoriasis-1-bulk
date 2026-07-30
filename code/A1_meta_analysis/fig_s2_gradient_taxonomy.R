# Auto-extracted generating script
# Produces: fig_s2_gradient_taxonomy.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: d1eaafd3-76cb-46a1-bf6d-a25a6e775973
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

library(data.table)
library(ggplot2)
library(patchwork)
library(msigdbr)
library(limma)

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")

cls <- readRDS("sample_classification.rds")
setDT(cls)

mat <- fread("clust_input/SRP165679.tsv")
genes <- mat$Genes
X <- as.matrix(mat[,-1]); rownames(X) <- genes
lab <- cls[srp=="SRP165679"]
keep <- lab[class %in% c("NN","PN","PP")]
X <- X[, keep$external_id]
grp <- factor(keep$class, levels=c("NN","PN","PP"))
stage <- as.integer(grp)-1L

des <- model.matrix(~ stage)
fit <- lmFit(X, des); fit <- eBayes(fit, trend=TRUE)
tt <- topTable(fit, coef="stage", number=Inf, sort.by="none")
trend <- data.table(gene=rownames(X), trend_slope=tt$logFC, trend_t=tt$t,
                    trend_P=tt$P.Value, trend_FDR=tt$adj.P.Val)

mu <- sapply(c("NN","PN","PP"), function(g) rowMeans(X[, grp==g, drop=FALSE]))
colnames(mu) <- c("mu_NN","mu_PN","mu_PP")
trend <- cbind(trend, mu)
trend[, mono_up   := mu_NN < mu_PN & mu_PN < mu_PP]
trend[, mono_down := mu_NN > mu_PN & mu_PN > mu_PP]
trend[, monotonic := mono_up | mono_down]
trend[, pn_frac := round((mu_PN-mu_NN)/(mu_PP-mu_NN), 3)]

de <- readRDS("per_study_de.rds")
ppnn <- as.data.table(de$PPvsNN[["SRP165679"]])
de_genes <- ppnn[adjP<0.05 & abs(logFC)>1, gene]

g <- trend[gene %in% de_genes]

Xg <- X[g$gene, ]
rho <- apply(Xg, 1, function(v) suppressWarnings(cor(v, stage, method="spearman")))
g[, spearman_rho := round(rho[gene],3)]

g[, cls := fifelse(pn_frac < 0, "PN_divergent",
            fifelse(pn_frac >= 0.50, "early_PN",
            fifelse(pn_frac >= 0.15, "progressive", "late_PP")))]

H <- as.data.table(msigdbr(species="Homo sapiens", category="H"))
Hsets <- split(H$gene_symbol, H$gs_name)
universe <- unique(trend$gene)

enr <- function(genes, topn=8){
  genes <- intersect(genes, universe)
  res <- rbindlist(lapply(names(Hsets), function(s){
    inset <- intersect(Hsets[[s]], universe); a <- length(intersect(genes, inset))
    b <- length(genes)-a; c <- length(inset)-a; d <- length(universe)-a-b-c
    data.table(set=sub("HALLMARK_","",s), overlap=a, set_size=length(inset),
               p=fisher.test(matrix(c(a,b,c,d),2), alternative="greater")$p.value)
  }))
  res[, FDR := p.adjust(p,"BH")]; res[order(p)][1:topn]
}
enr_prog <- enr(g[cls=="progressive", gene]); enr_prog[, class := "progressive (early/inflammatory)"]
enr_late <- enr(g[cls=="late_PP", gene]);      enr_late[, class := "late_PP (proliferation)"]
enr_all <- rbind(enr_prog, enr_late)

NNblue <- "#4C72B0"; PPred <- "#C44E52"; grey <- "#B0B0B0"
th <- theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold",size=13), legend.position="none")

enr_all[, lbl := gsub("_"," ", set)]
enr_all[, mlogFDR := -log10(FDR)]
enr_all[, class := factor(class, levels=c("progressive (early/inflammatory)","late_PP (proliferation)"))]
top <- enr_all[, .SD[order(-mlogFDR)][1:6], by=class]
top[, fillc := ifelse(grepl("progressive",class), "#55A868", "#8172B3")]

cl <- g[,.N,by=cls]
cl[, lab := c(early_PN="Early\n(PN-specific)", progressive="Progressive\n(NN<PN<PP)",
              late_PP="Late\n(PP-specific)", PN_divergent="PN-divergent")[cls]]
cl[, ord := c(early_PN=1, progressive=2, late_PP=3, PN_divergent=4)[cls]]
cl[, col := c(early_PN="#DD8452", progressive="#55A868", late_PP="#8172B3", PN_divergent=grey)[cls]]

pA <- ggplot(cl, aes(reorder(lab,ord), N, fill=col)) + geom_col(width=.7) + scale_fill_identity() +
  geom_text(aes(label=N), vjust=-0.3, size=4) + labs(title="A. Gradient-gene timing classes", x=NULL, y="genes") +
  scale_y_continuous(expand=expansion(mult=c(0,.12))) + th

pB <- ggplot(top, aes(reorder(lbl,mlogFDR), mlogFDR, fill=fillc)) + geom_col(width=.7) + coord_flip() +
  scale_fill_identity() + facet_wrap(~class, scales="free_y", ncol=1) +
  labs(title="B. What activates early vs late",
       subtitle="Progressive = interferon/inflammation (already rising in peri-lesional);\nLate = cell-cycle/proliferation (switches on only at lesional stage)",
       x=NULL, y="-log10 FDR (Hallmark enrichment)") +
  theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),
       plot.title=element_text(face="bold",size=13), plot.subtitle=element_text(size=9.5),
       strip.text=element_text(face="bold",size=10), legend.position="none")

fig <- pA + pB + plot_layout(widths=c(1,1.5))
ggsave("staging_axis/figures/fig_s2_gradient_taxonomy.png", fig, width=12, height=5.2, dpi=150)