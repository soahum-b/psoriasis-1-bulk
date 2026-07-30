# Auto-extracted generating script
# Produces: fig_s4_stat3_il17_on_axis.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): sample_classification.rds
# Source artifact version: 5eb4f901-0e0b-4b46-846d-d237517c97c3
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressPackageStartupMessages({library(ggplot2); library(data.table)})

setwd("/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3")

mat <- fread("clust_input/SRP165679.tsv"); X<-as.matrix(mat[,-1]); rownames(X)<-mat$Genes
cls <- as.data.table(readRDS("sample_classification.rds"))
lab <- cls[srp=="SRP165679"][class %in% c("NN","PN","PP")]
X <- X[, lab$external_id]; grp <- factor(lab$class, levels=c("NN","PN","PP"))
NNblue<-"#4C72B0"; PNorange<-"#DD8452"; PPred<-"#C44E52"
cols <- c(NN=NNblue, PN=PNorange, PP=PPred)

# canonical hubs
genes <- c("STAT3","STAT1","SOCS3","IL17A","IL23A","S100A7","DEFB4A","IL36G")
genes <- genes[genes %in% rownames(X)]
d <- rbindlist(lapply(genes, function(g)
  data.table(gene=g, expr=X[g,], grp=grp)))
d[, gene := factor(gene, levels=genes)]

th <- theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold",size=13), strip.text=element_text(face="bold",size=10),
      legend.position="none")
pA <- ggplot(d, aes(grp, expr, fill=grp))+
  geom_boxplot(outlier.size=.6, width=.65)+
  facet_wrap(~gene, scales="free_y", nrow=2)+
  scale_fill_manual(values=cols)+
  labs(title="A. Canonical psoriasis hubs rise monotonically along the staging axis",
       subtitle="Each validated gene sits intermediate in peri-lesional (PN) skin - the axis recovers the known biology",
       x="Staging axis (NN -> PN -> PP)", y="normalized log-expression")+
  th+theme(plot.subtitle=element_text(size=9.5))
ggsave("fig_s4_stat3_il17_on_axis.png", pA, width=11, height=5.6, dpi=150)
cat("saved fig_s4\n")