# Auto-extracted generating script
# Produces: fig_target_ranking.png
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: bb27fd14-1926-4771-9337-7929a809429e
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


BASE = "/Users/soahum/R-Projects/Kostka_project/project-kostka/Recount-3"
WP = "/tmp/psor_pathway"
os.makedirs("work", exist_ok=True)

meta = pd.read_csv(f"{BASE}/meta_de_PPvsNN.csv")

le = pd.read_csv(f"{WP}/leading_edge_long.csv")
le_genes = set(le.gene)
le_by_gene = le.groupby("gene")["set"].apply(lambda s: ";".join(sorted(set(s))))
stat3_targets = set(open(f"{WP}/stat3_targets.txt").read().split())
mod = pd.read_csv(f"{BASE}/clust_module_genes.csv")
mod_genes = set(mod.gene)

sig = meta[(meta.FDR < 0.05) & (meta.logFC.abs() > 1)].copy()
sig["in_leading_edge"] = sig.gene.isin(le_genes)
sig["leading_edge_sets"] = sig.gene.map(le_by_gene).fillna("")
sig["stat3_target"] = sig.gene.isin(stat3_targets)
sig["in_coexpr_module"] = sig.gene.isin(mod_genes)
sig["absLFC"] = sig.logFC.abs()
sig["neglog10FDR"] = -np.log10(sig.FDR.clip(lower=1e-300))
sig["robustness"] = sig.k * (1 - sig.I2 / 100.0)
sig["pathway_support"] = (sig.in_leading_edge.astype(int) + sig.stat3_target.astype(int) + sig.in_coexpr_module.astype(int))

def z(x):
    x = x.astype(float)
    return (x - x.mean()) / x.std(ddof=0)

sig["evidence_score"] = (
    0.35 * z(sig.absLFC) + 0.30 * z(sig.neglog10FDR) +
    0.15 * z(sig.robustness) + 0.20 * z(sig.pathway_support)
)
sig = sig.sort_values("evidence_score", ascending=False).reset_index(drop=True)
sig["evidence_rank"] = np.arange(1, len(sig) + 1)
sig.to_csv("work/target_universe.csv", index=False)

import json

sym2ens = json.load(open("handoff/sym2ens.json"))
ot = json.load(open("handoff/ot_targets.json"))
chembl = json.load(open("handoff/chembl_targets.json"))
bio = json.load(open("handoff/chembl_bioactivity.json"))

def tract_flags(t):
    sm = {x["label"]: x["value"] for x in t["tractability"] if x["modality"] == "SM"}
    ab = {x["label"]: x["value"] for x in t["tractability"] if x["modality"] == "AB"}
    if sm.get("Approved Drug"): sm_bucket = "Approved drug"
    elif sm.get("Advanced Clinical"): sm_bucket = "Advanced clinical"
    elif sm.get("Phase 1 Clinical"): sm_bucket = "Phase 1"
    elif sm.get("Structure with Ligand") or sm.get("High-Quality Ligand"): sm_bucket = "Discovery precedence"
    elif sm.get("High-Quality Pocket") or sm.get("Med-Quality Pocket"): sm_bucket = "Druggable pocket"
    elif sm.get("Druggable Family"): sm_bucket = "Druggable family"
    else: sm_bucket = "Unknown"
    ab_clin = ab.get("Approved Drug") or ab.get("Advanced Clinical") or ab.get("Phase 1 Clinical")
    return dict(
        sm_bucket=sm_bucket,
        sm_has_ligand=bool(sm.get("Structure with Ligand")),
        sm_hq_ligand=bool(sm.get("High-Quality Ligand")),
        sm_hq_pocket=bool(sm.get("High-Quality Pocket")),
        sm_mq_pocket=bool(sm.get("Med-Quality Pocket")),
        sm_druggable_family=bool(sm.get("Druggable Family")),
        ab_clinical=bool(ab_clin),
    )

bucket_rank = {"Approved drug": 6, "Advanced clinical": 5, "Phase 1": 4, "Discovery precedence": 3,
               "Druggable pocket": 2, "Druggable family": 1, "Unknown": 0}

rows = []
for sym, t in ot.items():
    tf = tract_flags(t)
    cls = "; ".join(sorted({c["label"] for c in t.get("targetClass", []) if c.get("level") == "l1"})) or "Unclassified"
    cb = chembl.get(sym) or {}
    bb = bio.get(sym) or {}
    de = sig[sig.gene == sym]
    de = de.iloc[0] if len(de) else None
    rows.append(dict(
        gene=sym, ensembl=t["id"], approvedName=t.get("approvedName", ""),
        protein_class=cls,
        logFC=(de.logFC if de is not None else np.nan),
        FDR=(de.FDR if de is not None else np.nan),
        direction=(de.direction if de is not None else ""),
        evidence_rank=(int(de.evidence_rank) if de is not None else np.nan),
        pathway_support=(int(de.pathway_support) if de is not None else 0),
        in_leading_edge=(bool(de.in_leading_edge) if de is not None else False),
        stat3_target=(bool(de.stat3_target) if de is not None else False),
        ot_sm_bucket=tf["sm_bucket"], ot_sm_bucket_rank=bucket_rank[tf["sm_bucket"]],
        sm_has_ligand=tf["sm_has_ligand"], sm_hq_ligand=tf["sm_hq_ligand"],
        sm_hq_pocket=tf["sm_hq_pocket"], sm_mq_pocket=tf["sm_mq_pocket"],
        sm_druggable_family=tf["sm_druggable_family"], ab_clinical=tf["ab_clinical"],
        ot_drug_candidates=t.get("drugAndClinicalCandidates", {}).get("count", 0),
        chembl_target_id=cb.get("target_chembl_id"),
        n_potent_inhib=bb.get("n_potent_pchembl7", 0),
        n_active_inhib=bb.get("n_active_pchembl6", 0),
    ))
ann = pd.DataFrame(rows)
ann.to_csv("work/druggability_annotation.csv", index=False)

def mm(x):
    x = x.astype(float); lo, hi = np.nanmin(x), np.nanmax(x)
    return (x - lo) / (hi - lo) if hi > lo else x * 0

ann["s_effect"] = mm(ann.logFC.abs())
ann["s_signif"] = mm(-np.log10(ann.FDR.clip(lower=1e-300)))
ann["s_pathway"] = mm(ann.pathway_support.astype(float))
ann["s_tract"] = mm(ann.ot_sm_bucket_rank.astype(float))
ann["s_chem"] = mm(np.log10(ann.n_potent_inhib.clip(lower=0) + 1))
ann["s_drugs"] = mm(np.log10(ann.ot_drug_candidates.clip(lower=0) + 1))
ann["s_struct"] = ((ann.sm_has_ligand.astype(int) + ann.sm_hq_ligand.astype(int)
                    + ann.sm_hq_pocket.astype(int) + ann.sm_mq_pocket.astype(int)) / 4.0)
ann["druggability_score"] = (
    0.18 * ann.s_effect + 0.12 * ann.s_signif +
    0.15 * ann.s_pathway +
    0.16 * ann.s_tract + 0.12 * ann.s_chem + 0.12 * ann.s_drugs +
    0.15 * ann.s_struct
)
ann = ann.sort_values("druggability_score", ascending=False).reset_index(drop=True)
ann["druggability_rank"] = np.arange(1, len(ann) + 1)

axis_map = {
    "Th17/IL-17 axis": ["RORC", "RORA", "IL17A", "IL17F", "IL17RA", "CCL20", "CCR6", "IL23A", "IL12B", "IL23R"],
    "JAK-STAT": ["JAK1", "JAK2", "JAK3", "TYK2", "STAT1", "STAT2", "STAT3", "IL6R", "IL22RA1", "IL22", "IL19"],
    "NF-kB/TNF": ["NFKB1", "NFKB2", "RELA", "RELB", "REL", "BIRC3", "TNIP3", "ZC3H12A", "IRAK2", "TNF"],
    "Innate/antimicrobial": ["S100A9", "S100A7", "S100A8", "S100A12", "DEFB4A", "DEFB4B", "IL36A", "IL36G", "IL36RN", "LTF", "PI3", "NOS2"],
    "Chemokine/leukocyte": ["CXCR1", "CXCR2", "CXCL13", "CCL2", "CCL20", "LTB4R"],
    "Tryptophan/immunometab": ["IDO1", "KYNU", "TDO2", "VNN3"],
    "Eicosanoid/protease": ["PTGS2", "ALOX5", "ALOX15", "MMP9", "MMP1", "MMP12", "PLA2G2A", "PLA2G4D", "LTA4H"],
    "Proliferation": ["PLK1", "AURKA", "AURKB", "CDK1", "RRM2", "TYMS", "BIRC5"],
}
g2axis = {}
for ax_, gs in axis_map.items():
    for g in gs: g2axis.setdefault(g, []).append(ax_)
ann["disease_axis"] = ann.gene.map(lambda g: "; ".join(g2axis.get(g, [])))
ann["on_psoriasis_axis"] = ann.disease_axis.str.len() > 0
ann["struct_ready"] = (ann.sm_has_ligand | ann.sm_hq_pocket | ann.sm_mq_pocket)
ann.to_csv("work/druggable_target_ranking.csv", index=False)

LEADS = ["STAT3", "RORC", "JAK3"]
ann["is_lead"] = ann.gene.isin(LEADS)

bucket_col = {"Approved drug": "#1a9850", "Advanced clinical": "#66bd63", "Phase 1": "#a6d96a",
              "Discovery precedence": "#fdae61", "Druggable pocket": "#fee08b", "Druggable family": "#f6d", "Unknown": "#cccccc"}

apply_figure_style()
onaxis_all = ann[ann.on_psoriasis_axis].sort_values("druggability_score", ascending=False)
top17 = onaxis_all.head(17)
leads_df = ann[ann.gene.isin(LEADS)]
show = pd.concat([top17, leads_df]).drop_duplicates("gene").sort_values("druggability_score", ascending=False)
show = show.iloc[::-1]
colors = show.ot_sm_bucket.map(bucket_col).fillna("#cccccc")

fig, ax = plt.subplots(figsize=(7.8, 6.6))
ypos = np.arange(len(show))
ax.barh(ypos, show.druggability_score, color=colors, edgecolor="white", lw=0.5, zorder=3)
for y, (_, r) in zip(ypos, show.iterrows()):
    lab = r.gene + ("  ★" if r.is_lead else "")
    ax.text(-0.008, y, lab, ha="right", va="center", fontsize=7,
            fontweight=("bold" if r.is_lead else "normal"), fontstyle='italic',
            color=("#b30000" if r.is_lead else "black"))
    ax.text(r.druggability_score + 0.005, y, f"{int(r.n_potent_inhib):,} inh", ha="left", va="center", fontsize=5.6, color="0.35")
ax.set_yticks([]); ax.set_ylim(-0.7, len(show) - 0.3)
ax.set_xlabel("Composite druggability score")
ax.set_title("Top psoriasis-axis druggable targets  (★ = selected lead for docking/MD)", fontsize=8.3)
ax.spines['left'].set_visible(False); ax.margins(x=0.02)
handles = [Patch(facecolor=bucket_col[b], label=b) for b in ["Approved drug", "Advanced clinical", "Phase 1", "Discovery precedence", "Druggable pocket"]]
ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=6, title="SM tractability", title_fontsize=6.2)
fig.tight_layout()
fig.savefig("fig_target_ranking.png", dpi=200)
print("regenerated; genes shown:", show.gene.tolist())