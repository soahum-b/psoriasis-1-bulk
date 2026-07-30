#!/usr/bin/env Rscript
# Staging-axis figures (fig_s2 .. fig_s5b). Run AFTER pipeline STEP 12, which
# writes the objects these plots consume into staging_axis/results/.
# Usage:  conda run -n psoriasis-r Rscript staging_axis/code/staging_figures.R
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(pheatmap); library(grid)
})
setwd(Sys.getenv("STAGING_ROOT", "."))   # run from project root
NN<-"#4C72B0"; PN<-"#DD8452"; PP<-"#C44E52"
cols <- c(NN=NN, PN=PN, PP=PP)

trend <- readRDS("staging_axis/results/trend_SRP165679.rds")
g     <- fread("staging_axis/results/gradient_gene_classes.csv")
enr   <- readRDS("staging_axis/results/class_enrichment.rds")
obj   <- readRDS("staging_axis/results/ssgsea_3group.rds"); sc<-obj$scores; grp<-obj$grp
s3    <- readRDS("staging_axis/results/stat3_activity_3group.rds")

# --- fig_s2: taxonomy (class sizes + enrichment) ---
cl <- g[,.N,by=cls]
cl[, lab := c(early_PN="Early\n(PN-specific)",progressive="Progressive\n(NN<PN<PP)",
              late_PP="Late\n(PP-specific)",PN_divergent="PN-divergent")[cls]]
cl[, ord := c(early_PN=1,progressive=2,late_PP=3,PN_divergent=4)[cls]]
cl[, col := c(early_PN="#DD8452",progressive="#55A868",late_PP="#8172B3",PN_divergent="#B0B0B0")[cls]]
th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold",size=13),legend.position="none")
pA <- ggplot(cl,aes(reorder(lab,ord),N,fill=col))+geom_col(width=.7)+scale_fill_identity()+
  geom_text(aes(label=N),vjust=-0.3,size=4)+labs(title="A. Gradient-gene timing classes",x=NULL,y="genes")+
  scale_y_continuous(expand=expansion(mult=c(0,.12)))+th
enr[, lbl:=gsub("_"," ",set)][, mlogFDR:=-log10(FDR)]
enr[, class:=factor(class,levels=c("progressive (early/inflammatory)","late_PP (proliferation)"))]
top <- enr[,.SD[order(-mlogFDR)][1:6],by=class]
top[, fillc:=ifelse(grepl("progressive",class),"#55A868","#8172B3")]
pB <- ggplot(top,aes(reorder(lbl,mlogFDR),mlogFDR,fill=fillc))+geom_col(width=.7)+coord_flip()+
  scale_fill_identity()+facet_wrap(~class,scales="free_y",ncol=1)+
  labs(title="B. What activates early vs late",x=NULL,y="-log10 FDR (Hallmark)")+
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold",size=13),strip.text=element_text(face="bold"),legend.position="none")
ggsave("staging_axis/figures/fig_s2_gradient_taxonomy.png", pA+pB+plot_layout(widths=c(1,1.5)),
       width=12,height=5.2,dpi=150)

# --- fig_s3: pathway timing lineplot ---
plotd <- rbindlist(lapply(rownames(sc),function(p) data.table(pathway=p,
  stage=factor(c("NN","PN","PP"),levels=c("NN","PN","PP")),
  mean=tapply(sc[p,],grp,mean), se=tapply(sc[p,],grp,function(z)sd(z)/sqrt(length(z))))))
immune<-c("INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","IL6_JAK_STAT3_SIGNALING",
          "INFLAMMATORY_RESPONSE","TNFA_SIGNALING_VIA_NFKB","ALLOGRAFT_REJECTION")
prolif<-c("E2F_TARGETS","G2M_CHECKPOINT","MYC_TARGETS_V1","MTORC1_SIGNALING")
plotd<-plotd[pathway %in% c(immune,prolif)]
plotd[, program:=ifelse(pathway %in% immune,"Immune / interferon (early)","Proliferation (late)")]
plotd[, lbl:=gsub("_"," ",pathway)]
ggsave("staging_axis/figures/fig_s3_pathway_timing.png",
  ggplot(plotd,aes(stage,mean,group=lbl,colour=program))+geom_line(linewidth=.8,alpha=.8)+
    geom_point(size=2)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=.12,alpha=.6)+
    facet_wrap(~program)+scale_colour_manual(values=c("Immune / interferon (early)"="#55A868","Proliferation (late)"="#8172B3"))+
    labs(title="Pathway activation timing along the staging axis (SRP165679)",
         x="Staging axis",y="ssGSEA score (mean +/- SE)")+
    theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position="none",
      plot.title=element_text(face="bold",size=13),strip.text=element_text(face="bold")),
  width=10,height=4.8,dpi=150)

# --- fig_s4 / fig_s4b: canonical hubs + STAT3 activity ---
mat<-fread("clust_input/SRP165679.tsv"); X<-as.matrix(mat[,-1]); rownames(X)<-mat$Genes
cls_s<-as.data.table(readRDS("sample_classification.rds"))
lab<-cls_s[srp=="SRP165679"][class %in% c("NN","PN","PP")]; X<-X[,lab$external_id]
grp2<-factor(lab$class,levels=c("NN","PN","PP"))
genes<-intersect(c("STAT3","STAT1","SOCS3","IL17A","IL23A","S100A7","DEFB4A","IL36G"),rownames(X))
d<-rbindlist(lapply(genes,function(gg)data.table(gene=factor(gg,levels=genes),expr=X[gg,],grp=grp2)))
ggsave("staging_axis/figures/fig_s4_stat3_il17_on_axis.png",
  ggplot(d,aes(grp,expr,fill=grp))+geom_boxplot(outlier.size=.6,width=.65)+
    facet_wrap(~gene,scales="free_y",nrow=2)+scale_fill_manual(values=cols)+
    labs(title="A. Canonical psoriasis hubs rise monotonically along the staging axis",
         x="Staging axis (NN -> PN -> PP)",y="normalized log-expression")+
    theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="none",
      plot.title=element_text(face="bold",size=13),strip.text=element_text(face="bold")),
  width=11,height=5.6,dpi=150)
tp<-summary(lm(s3$score~(as.integer(s3$grp)-1L)))$coefficients[2,4]
ggsave("staging_axis/figures/fig_s4b_stat3_activity_axis.png",
  ggplot(s3,aes(grp,score,fill=grp))+geom_boxplot(width=.6,outlier.size=.6)+scale_fill_manual(values=cols)+
    labs(title="B. STAT3 transcription-factor activity on the axis",
         subtitle=sprintf("decoupleR ULM vs healthy baseline; monotonic rise (trend p=%.0e)",tp),
         x="Staging axis (NN -> PN -> PP)",y="STAT3 activity (ULM score)")+
    theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="none",
      plot.title=element_text(face="bold",size=13)),
  width=6,height=4.6,dpi=150)

# --- fig_s5 / fig_s5b: staging heatmap + monotonicity summary ---
g[, flag:=ifelse(mu_NN<mu_PN & mu_PN<mu_PP,"up",ifelse(mu_NN>mu_PN & mu_PN>mu_PP,"down","non"))]
up<-g[flag=="up"][order(-trend_slope)][1:40,gene]; dn<-g[flag=="down"][order(trend_slope)][1:40,gene]
Z<-t(scale(t(X[c(up,dn),]))); ord<-order(grp2,colMeans(Z[1:40,,drop=FALSE]))
Z<-Z[,ord]; Z[Z>2.5]<-2.5; Z[Z< -2.5]<- -2.5
ac<-data.frame(Stage=grp2[ord]); rownames(ac)<-colnames(Z)
ar<-data.frame(Direction=rep(c("up NN<PN<PP","down NN>PN>PP"),c(40,40))); rownames(ar)<-rownames(Z)
png("staging_axis/figures/fig_s5_staging_heatmap.png",width=1500,height=1150,res=150)
pheatmap(Z,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=ac,annotation_row=ar,
  annotation_colors=list(Stage=cols,Direction=c("up NN<PN<PP"="#C44E52","down NN>PN>PP"="#4C72B0")),
  show_colnames=FALSE,fontsize_row=5.5,gaps_row=40,
  color=colorRampPalette(c("#3B4CC0","white","#B40426"))(101),
  main="Molecular staging axis: top monotonic genes ordered NN -> PN -> PP (SRP165679)")
dev.off()
sm<-data.table(cat=factor(c("Strictly\nmonotonic","Monotonic\nrank (|rho|>0.3)","Non-\nmonotonic"),
  levels=c("Strictly\nmonotonic","Monotonic\nrank (|rho|>0.3)","Non-\nmonotonic")),
  pct=c(100*g[flag!="non",.N]/nrow(g),100*g[abs(spearman_rho)>0.3,.N]/nrow(g),100*g[flag=="non",.N]/nrow(g)),
  col=c("#55A868","#4C72B0","#B0B0B0"))
pf<-g[flag!="non" & pn_frac>=0 & pn_frac<=1]; med<-median(pf$pn_frac)
qA<-ggplot(sm,aes(cat,pct,fill=col))+geom_col(width=.7)+scale_fill_identity()+
  geom_text(aes(label=sprintf("%.0f%%",pct)),vjust=-0.3,size=4.2)+
  labs(title="A. Monotonicity of the psoriasis program",x=NULL,y="% of DE genes")+
  scale_y_continuous(expand=expansion(mult=c(0,.12)),limits=c(0,105))+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=13))
qB<-ggplot(pf,aes(pn_frac))+geom_histogram(bins=40,fill="#DD8452",colour="white",linewidth=.2)+
  geom_vline(xintercept=med,linetype=2,colour="#C44E52",linewidth=.7)+
  labs(title="B. Peri-lesional sits early on the axis",
       x="PN fraction toward lesional (0=healthy-like, 1=lesional-like)",y="genes")+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=13))
ggsave("staging_axis/figures/fig_s5b_monotonicity_summary.png",qA+qB+plot_layout(widths=c(1,1.3)),
       width=11,height=4.4,dpi=150)
cat("staging figures written to staging_axis/figures/\n")
