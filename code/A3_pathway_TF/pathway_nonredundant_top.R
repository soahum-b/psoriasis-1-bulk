# Auto-extracted generating script
# Produces: pathway_nonredundant_top.csv
# Conda env: psoriasis-r   (run with this environment activated)
# Inputs (expected alongside / in data/): collapse_primary.rds
# Source artifact version: 34744df6-93a9-49e9-8bee-23b367399056
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

suppressMessages({library(fgsea); library(msigdbr); library(data.table); library(jsonlite)})

obj <- readRDS("collapse_primary.rds"); fg <- obj$fg; main <- obj$main
setDT(fg)
coll_of <- function(pw) fifelse(grepl("^HALLMARK_",pw),"Hallmark",
                    fifelse(grepl("^REACTOME_",pw),"Reactome",
                    fifelse(grepl("^GOBP_",pw),"GO:BP",
                    fifelse(grepl("^KEGG_",pw),"KEGG","other"))))
nr <- fg[pathway %in% main & NES>0][order(padj_joint,-NES)]
nr[, collection := coll_of(pathway)]
nr[, label := sub("^(HALLMARK|REACTOME|GOBP|KEGG)_","",pathway)]

# coarse theme by keyword (transparent, rule-based)
theme_of <- function(l){
  l <- toupper(l)
  fcase(
    grepl("CELL_CYCLE|MITOT|G2M|E2F|CHROMOSOM|SPINDLE|DNA_REPLICAT|SISTER|KINETOCHORE|M_PHASE|CHECKPOINT",l),"Proliferation / cell cycle",
    grepl("INTERFERON|VIRUS|ANTIVIRAL|IFN|ISG|INNATE_IMMUN|PATTERN_RECOG",l),"Interferon / antiviral",
    grepl("IL6|JAK|STAT|IL_?17|IL17|TNF|NFKB|NF_KAPPA|INFLAMMAT|CYTOKINE|CHEMOKINE|INTERLEUKIN",l),"Inflammatory signaling (IL/JAK-STAT/TNF)",
    grepl("NEUTROPHIL|DEFENS|ANTIMICROB|COMPLEMENT|ALLOGRAFT|LYMPHOCYTE|T_CELL|LEUKOCYTE|IMMUNE_RESPONSE",l),"Immune effector / antimicrobial",
    grepl("RRNA|RIBONUCLEO|TRANSLATION|RIBOSOM|MRNA|SPLICEOS|RNA_PROCESS|NONSENSE",l),"RNA processing / translation",
    grepl("MYC_TARGET|MTORC1|OXIDATIVE_PHOS|GLYCOLYS|METABOL|RESPIRAT",l),"Metabolism / MYC-mTOR",
    grepl("KERATIN|CORNIF|EPIDERM|SKIN_DEV|DESMOSOME",l),"Keratinocyte / epidermis",
    default="Other")
}
nr[, theme := theme_of(label)]

fwrite(nr[, .(theme, pathway, label, collection, NES, ES, padj_joint, padj_percollection=padj, size)],
       "pathway_nonredundant_top.csv")