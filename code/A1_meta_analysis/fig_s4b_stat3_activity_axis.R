# Auto-extracted generating script
# Produces: fig_s4b_stat3_activity_axis.png
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): collectri_regulon.rds, sample_classification.rds
# Source artifact version: 3c64831f-6edb-4bb2-93be-25d16b2165aa
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({library(decoupleR); library(ggplot2); library(data.table)})

mat <- fread("clust_input/SRP165679.tsv"); X<-as.matrix(mat[,-1]); rownames(X)<-mat$Genes
cls <- as.data.table(readRDS("sample_classification.rds"))
lab <- cls[srp=="SRP165679"][class %in% c("NN","PN","PP")]
X <- X[, lab$external_id]; grp <- factor(lab$class, levels=c("NN","PN","PP"))
reg <- readRDS("collectri_regulon.rds"); reg <- as.data.table(reg)

nnmean <- rowMeans(X[, grp=="NN"])
dev <- X - nnmean
net <- reg[, .(source, target, mor)]
net <- net[target %in% rownames(dev)]
act <- run_ulm(mat=dev, net=net, .source="source", .target="target", .mor="mor", minsize=5)
setDT(act)
s3 <- act[source=="STAT3"]
s3 <- merge(s3, data.table(condition=colnames(dev), grp=grp), by="condition")

stage <- as.integer(s3$grp)-1L
tp <- summary(lm(s3$score~stage))$coefficients["stage",4]
mu <- tapply(s3$score, s3$grp, mean)

cols <- c(NN="#4C72B0", PN="#DD8452", PP="#C44E52")
pB <- ggplot(s3, aes(grp, score, fill=grp))+geom_boxplot(width=.6, outlier.size=.6)+
  scale_fill_manual(values=cols)+
  labs(title="B. STAT3 transcription-factor activity on the axis",
       subtitle=sprintf("decoupleR ULM (CollecTRI regulon) vs healthy baseline; monotonic rise (trend p=%.0e)",tp),
       x="Staging axis (NN -> PN -> PP)", y="STAT3 activity (ULM score)")+
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(), legend.position="none",
       plot.title=element_text(face="bold",size=13), plot.subtitle=element_text(size=9))
ggsave("fig_s4b_stat3_activity_axis.png", pB, width=6, height=4.6, dpi=150)