# Makefile for the psoriasis meta-analysis.
# Run  `make help`  to see available targets.
#
# Recipe lines are TAB-indented (a Make requirement). If you edit this file,
# keep the tabs - spaces will produce "missing separator" errors.

R_ENV     := psoriasis-r
CLUST_ENV := clust-env
PIPELINE  := psoriasis_pipeline.R

.PHONY: help setup pipeline clust all verify clean

## help: list available targets
help:
	@echo "Psoriasis meta-analysis - make targets:"
	@echo "  make setup      Create both conda environments (R + clust)"
	@echo "  make pipeline   Run the full R analysis pipeline"
	@echo "  make clust      Run the clust cross-study co-expression step"
	@echo "  make all        setup -> pipeline -> clust (everything)"
	@echo "  make verify     Check that R packages and clust are importable"
	@echo "  make clean      Remove regenerable outputs (figures, clust dirs)"

## setup: create the psoriasis-r and clust-env conda environments
setup:
	bash setup.sh

## pipeline: run the full R analysis (STEP 1-11)
pipeline:
	conda run -n $(R_ENV) Rscript $(PIPELINE)

## clust: cross-study co-expression modules (needs clust_input/ from the pipeline)
clust:
	conda run -n $(CLUST_ENV) python -m clust clust_input -o clust_out -t 0

## all: full run from scratch
all: setup pipeline clust

## verify: confirm the environments have what they need
verify:
	@conda run -n $(R_ENV) Rscript -e 'pkgs <- c("recount3","edgeR","limma","fgsea","clusterProfiler","decoupleR","qvalue","org.Hs.eg.db","msigdbr","metafor","UpSetR","data.table","ggplot2","pheatmap"); miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(miss)) cat("MISSING:", paste(miss, collapse=", "), "\n") else cat("All", length(pkgs), "R packages OK.\n")'
	@conda run -n $(CLUST_ENV) python -c "import clust; print('clust import OK')"

## clean: delete regenerable outputs (keeps code, docs, checkpoints)
clean:
	rm -rf clust_input clust_input2 clust_input3 clust_out clust_out2 clust_out3
	rm -f fig*.png
	@echo "Removed regenerable figures and clust directories."
