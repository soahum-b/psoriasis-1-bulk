# 7. Functional Enrichment Analysis ----

# Load packages ----
library(tidyverse)
library(limma)
library(DT)
library(gplots)
library(gprofiler2)

library(GSEABase)
library(GSVA)

library(enrichplot)
library(clusterProfiler)


# 1. Gene Ontology (GO) Enrichment via gProfiler2 (Over-Representation)-----

## 1. Use the specific UP-regulated module from Code 5
# We use the rownames (gene symbols) from myModule_2
gost.res_module <- gost(rownames(myModule_2), organism = "hsapiens", correction_method = "fdr")

# 2. Produce an interactive Manhattan plot
gostplot(gost.res_module, interactive = TRUE, capped = FALSE)

# 
# 2 Gene Set Enrichment Analysis (GSEA) via clusterProfiler-----
# 1. Load the MSigDB C7 Collection (Immunologic Signatures)
c7cp_file <- read.gmt("/net/dali/home/mscbio/sba50/Soahum/Recount3_project/msig/c7.all.v2026.1.Hs.symbols.gmt")

# 2. Prepare the Ranked List for GSEA
# We extract ALL genes, unfiltered, to build the ranked profile.
# We rank by the t-statistic (accounts for LogFC magnitude AND variance/p-value)
all_genes_for_gsea <- topTable(ebFit, coef = 1, number = Inf, sort.by = "none")

# Construct the named vector
mydata.gsea <- all_genes_for_gsea$t
names(mydata.gsea) <- rownames(all_genes_for_gsea)

# Sort strictly in decreasing order (Required for GSEA)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# 3. Run GSEA
set.seed(123) 
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE = c7cp_file, verbose = FALSE, pvalueCutoff = 0.05)
myGSEA.df <- as_tibble(myGSEA.res@result)

# 4. Visualize GSEA Results
# Bubble plot of top 20 pathways
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "PP (Lesional)",
    NES < 0 ~ "NN (Normal)"
  ))

ggplot(myGSEA.df[1:20, ], aes(x = phenotype, y = ID)) +
  geom_point(aes(size = setSize, color = NES, alpha = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  labs(title = "Top 20 Enriched Immunologic Signatures")

#3 Competitive GSEA via CAMERA (limma) -----
# 1. Read GMT file specifically for CAMERA (Requires a GeneSetCollection / List)
broadSet.C7.all <- getGmt("/net/dali/home/mscbio/sba50/Soahum/Recount3_project/msig/c7.all.v2026.1.Hs.symbols.gmt", geneIdType = SymbolIdentifier())
broadSet.C7.all <- geneIds(broadSet.C7.all) # Extract as list

# 2. Run CAMERA
# Uses the expression matrix (v.DEGList$E) directly
camera.res <- camera(v.DEGList$E, broadSet.C7.all, design, contrast.matrix[,1])
camera.df <- as_tibble(camera.res, rownames = "setName") %>%
  dplyr::filter(FDR <= 0.01) %>%
  mutate(phenotype = case_when(
    Direction == "Up" ~ "PP (Lesional)",
    Direction == "Down" ~ "NN (Normal)"
  ))


# 4. Single Sample GSEA via GSVA-------

# 1. Prepare the Signature List (Fixing the GSVA API requirement)
# Split the clusterProfiler dataframe into a named list of character vectors
C7CP_list <- split(c7cp_file$gene, c7cp_file$term)
C7CP_final <- lapply(C7CP_list, as.character)

# 2. Build the GSVA Parameter Object
GSVA_param_obj <- gsvaParam(exprData = v.DEGList$E, 
                            geneSets = C7CP_final, 
                            minSize = 5, 
                            maxSize = 500, 
                            maxDiff = FALSE)

# 3. Run the GSVA Algorithm (Transforms genes to pathways)
GSVA_results_matrix <- gsva(GSVA_param_obj)

# 4. Apply Linear Model (Limma) to the GSVA output
fit.C7CP <- lmFit(GSVA_results_matrix, design)
fits.C7CP_contrast <- contrasts.fit(fit.C7CP, contrast.matrix)
ebFit.C7CP <- eBayes(fits.C7CP_contrast)

# 5. Extract Significant Pathways
# We drop the lfc threshold to 0 here to capture any significant pathway shift
res.C7CP <- decideTests(ebFit.C7CP, method = "global", adjust.method = "BH", p.value = 0.05, lfc = 0)
diffSets.C7CP <- GSVA_results_matrix[res.C7CP[,1] != 0, ]

# 6. Heatmap of Differentially Enriched Pathways
# Cluster the significant GSVA scores
hr.C7CP <- hclust(as.dist(1 - cor(t(diffSets.C7CP), method = "pearson")), method = "complete")
hc.C7CP <- hclust(as.dist(1 - cor(diffSets.C7CP, method = "spearman")), method = "complete")

myheatcol <- colorRampPalette(colors = c("blue", "white", "red"))(100)

heatmap.2(diffSets.C7CP,
          Rowv = as.dendrogram(hr.C7CP),
          Colv = as.dendrogram(hc.C7CP),
          col = myheatcol, scale = "row",
          density.info = "none", trace = "none",
          cexRow = 0.5, cexCol = 0.5, margins = c(8,14),
          main = "GSVA: Enriched Pathways per Sample")


