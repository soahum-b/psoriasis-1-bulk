# Auto-extracted generating script
# Produces: lncRNA_curated_STAT3_IL1_linkage.csv
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): meta_de_PPvsNN_4study_biotyped.csv, lncRNA_significant_ranked.csv, coexpr_lnc_target_adjusted.csv
# Source artifact version: 387f62b2-3703-4093-868d-1ae919f99409
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np

bt = pd.read_csv("meta_de_PPvsNN_4study_biotyped.csv")
rk = pd.read_csv("lncRNA_significant_ranked.csv")
adj = pd.read_csv("coexpr_lnc_target_adjusted.csv")

curated = [
 ("SH3PXD2A-AS1","STAT3/SH3PXD2A-AS1/miR-125b/STAT3 positive-feedback loop; enhances STAT3 expression","STAT3","psoriasis","up","Yang 2021 Cytokine 144:155535"),
 ("AGAP2-AS1","ceRNA to AKT3 via miR-424-5p (psoriasis); binds WTAP to stabilize STAT3 mRNA → IL6/STAT3 (cancer)","STAT3 / AKT","psoriasis + cancer","up","Implications review 2023; Zhang 2023 (GC)"),
 ("MALAT1","miR-330-5p/S100A7 and miR-125b/BRD4 (IL-17) axes; modulates Th17; targets NF-κB/STAT3/PI3K-AKT","STAT3 / Th17 / output","psoriasis","up","Zhou 2022 Autoimmunity; Wang 2024 JBN"),
 ("XIST","Upregulated in psoriatic lesions; miRNA sponge, promotes keratinocyte hyperproliferation","STAT3 / NF-kB (review)","psoriasis","up(review)","ncRNA 2025 review"),
 ("HOTAIR","Upregulated; promotes proliferation, inhibits apoptosis, amplifies inflammation","STAT3 / NF-kB (review)","psoriasis","up","ncRNA 2025 review"),
 ("MIR31HG","Promotes keratinocyte hyperproliferation and inflammation","NF-kB","psoriasis","up","review"),
 ("NEAT1","Downregulated; tumor-suppressor-like, pro-apoptotic/anti-inflammatory in epidermis","(homeostasis)","psoriasis","down","ncRNA 2025 review"),
 ("MEG3","MEG3/miR-21/caspase-8; inhibits keratinocyte proliferation (down in psoriasis)","AKT/apoptosis","psoriasis","down","Jia 2019"),
 ("H19","H19/miR-766-3p/S1PR3 via AKT/mTOR; down in psoriasis promotes proliferation","AKT/mTOR (S1PR3-STAT3 upstream)","psoriasis","down","He 2021"),
 ("HOXC13-AS","Skin-specific; regulates epidermal differentiation via COPA/Golgi-ER transport","(differentiation)","psoriasis/skin","up","Zhang 2023 CDD"),
 ("UCA1","Suppresses METTL14, activates HIF-1a/NF-kB; promotes keratinocyte inflammation","NF-kB","psoriasis","up","Hu 2023 CDD"),
 ("LINC00958","Induced in psoriasis epidermis, modulates epidermal proliferation; JAK/STAT3 in cancer","STAT3 (cancer) / proliferation","psoriasis","up","Sobolev 2022 JID; Chen 2020 TSCC"),
 ("LINC01094","IL6-JAK-STAT3 GSEA association; PTEN/AKT via AZGP1 (cancer)","STAT3 / AKT","cancer","up","Chen 2023 Cancers"),
 ("LINC01206","Enhances keratinocyte cell-cycle progression via EHF (overexpressed in psoriasis)","(cell cycle)","psoriasis","up","2024 (LINC01026/01206)"),
 ("MIR181A2HG","miR-223-3p/SOX6 sponge; inhibits keratinocyte proliferation (down in psoriasis)","(proliferation)","psoriasis","down","2024"),
 ("KLHDC7B-DT","ILF2-dependent; hyperproliferation and skin inflammation","(proliferation/inflam)","psoriasis","up","Yin 2022 Front Genet"),
]
cur = pd.DataFrame(curated, columns=["lncRNA","curated_mechanism","curated_pathway","disease_context","direction_reported","source"])

de_small = bt[["gene","logFC","FDR","I2","k","is_lncRNA","biotype"]].rename(columns={"gene":"lncRNA"})
cur = cur.merge(de_small, on="lncRNA", how="left")
stat3 = adj[adj.target=="STAT3"][["lncRNA","meta_r","meta_fdr","k_studies"]].rename(
    columns={"meta_r":"STAT3_coexpr_r","meta_fdr":"STAT3_coexpr_fdr","k_studies":"STAT3_coexpr_k"})
cur = cur.merge(stat3, on="lncRNA", how="left")
cur["in_meta"] = cur["logFC"].notna()
cur["my_direction"] = cur["logFC"].apply(lambda x: "up_in_PP" if pd.notna(x) and x>0 else ("down_in_PP" if pd.notna(x) else "not_quantified"))
cur["my_FDR_sig"] = cur["FDR"]<0.05

def concord(row):
    if not row["in_meta"]: return "not_quantified"
    rep = row["direction_reported"]
    mine = "up" if row["logFC"]>0 else "down"
    rep_simple = "up" if "up" in rep else ("down" if "down" in rep else "?")
    if rep_simple=="?": return "review_unclear"
    return "concordant" if rep_simple==mine else "DISCORDANT"
cur["direction_concordance"] = cur.apply(concord, axis=1)

cur["coexpr_STAT3_confirms"] = (cur["STAT3_coexpr_fdr"]<0.05) & (cur["STAT3_coexpr_r"].abs()>0.3)

def tier(row):
    has_stat3_curated = "STAT3" in str(row["curated_pathway"])
    psor = "psoriasis" in str(row["disease_context"])
    triple = has_stat3_curated and psor and row["my_FDR_sig"]==True and row["coexpr_STAT3_confirms"]==True
    if triple: return "1_triple_convergent"
    if has_stat3_curated and row["my_FDR_sig"] and row["coexpr_STAT3_confirms"]: return "2_curated+DE+coexpr"
    if row["my_FDR_sig"] and (row["coexpr_STAT3_confirms"] or has_stat3_curated): return "3_two_lines"
    if row["my_FDR_sig"]: return "4_DE_only"
    return "5_curated_only_not_DE"
cur["evidence_tier"] = cur.apply(tier, axis=1)

out_cols = ["lncRNA","evidence_tier","curated_pathway","curated_mechanism","disease_context",
            "direction_reported","in_meta","logFC","my_direction","my_FDR_sig","FDR","I2","k",
            "STAT3_coexpr_r","STAT3_coexpr_fdr","STAT3_coexpr_k","coexpr_STAT3_confirms",
            "direction_concordance","source"]
cur_out = cur[out_cols].sort_values("evidence_tier")
cur_out.to_csv("lncRNA_curated_STAT3_IL1_linkage.csv", index=False)