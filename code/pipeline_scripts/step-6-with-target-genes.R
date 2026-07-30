# 7. Functional Enrichment Analysis: Streamlined & Corrected ----

# 1. Setup & Load Libraries -----
library(tidyverse)
library(limma)
library(gprofiler2)
library(clusterProfiler)
library(enrichplot)
library(ggplotify) # Critical for fixing the alignment bug
library(GSVA)
library(gplots)

# 2. Gene Ontology (GO) Enrichment (gProfiler2) -----
# Using the significant module from Section 5
gost.res <- gost(query = rownames(myModule_2), 
                 organism = "hsapiens", 
                 sources = c("GO:BP", "KEGG", "REAC"),
                 correction_method = "fdr")

if (!is.null(gost.res)) {
  # Use standard plot, wrapped in as.ggplot to avoid alignment errors
  p_man <- gostplot(gost.res, interactive = FALSE) + 
           ggtitle("Manhattan Plot: GO & Pathway Enrichment")
  print(as.ggplot(p_man))
}

# 3. Gene Set Enrichment Analysis (GSEA) -----
# Load GMT (Immunologic C7)
c7_gmt <- read.gmt("/net/dali/home/mscbio/sba50/Soahum/Recount3_project/msig/c7.all.v2026.1.Hs.symbols.gmt")

# Prepare Ranked List
all_genes <- topTable(ebFit, coef = 1, number = Inf, sort.by = "none")
gene_list <- setNames(all_genes$t, rownames(all_genes))
gene_list <- sort(gene_list, decreasing = TRUE)

# Run GSEA
set.seed(123)
gsea_res <- GSEA(gene_list, TERM2GENE = c7_gmt, pvalueCutoff = 0.05, verbose = FALSE)
gsea_df <- as_tibble(gsea_res@result)

# 4. Targeted Exploration (STAT3, IL-17, NF-kB, EFTUD2) -----

# Identify pathways containing your key drivers
target_pathways <- gsea_df %>%
  filter(str_detect(ID, "IL17|STAT3|NFKB|INFLAM")) %>%
  slice_max(abs(NES), n = 3)

# visualization: GSEA Running Score
# We use 'gseaplot' (standard) wrapped in as.ggplot for stability
if(nrow(target_pathways) > 0) {
  p_gsea <- gseaplot(gsea_res, 
                     geneSetID = target_pathways$ID[1], 
                     title = target_pathway$ID[1])
  print(as.ggplot(p_gsea))
}

# visualization: Category Net Plot (cnetplot)
# This shows how STAT3/IL-17 connect across different signatures
p_cnet <- cnetplot(gsea_res, 
                   showCategory = 5, 
                   foldChange = gene_list) + 
          ggtitle("Network Analysis: Core Pathway Hubs")
print(as.ggplot(p_cnet))

# 5. Global Visual Summary -----

# Ridge Plot: Distribution of fold changes in top pathways
print(ridgeplot(gsea_res, showCategory = 15) + theme_bw())

# Dot Plot: Summary of enrichment
print(dotplot(gsea_res, showCategory = 10, split = ".sign") + facet_grid(.~.sign))

# 6. Save to PDF (Clean Output) -----
pdf("Psoriasis_Enrichment_Summary.pdf", width = 12, height = 10)
if (exists("p_man")) print(as.ggplot(p_man))
if (exists("p_gsea")) print(as.ggplot(p_gsea))
if (exists("p_cnet")) print(as.ggplot(p_cnet))
dev.off()