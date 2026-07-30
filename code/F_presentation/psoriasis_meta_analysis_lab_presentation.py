# Auto-extracted generating script
# Produces: psoriasis_meta_analysis_lab_presentation.pptx
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): none
# Source artifact version: d98ea726-13dc-475d-8d96-b9534460eb56
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

# [lineage] reconstruction failed validation (max_tokens) — raw producing cell below; see Execution Log for full trace
# ---------- SLIDE 37: Paper plan ----------
s = add_slide()
header(s, "PAPER PLAN — FOR DISCUSSION", "Proposed submission strategy")
# Paper 1 box
rect(s, Inches(0.5), Inches(1.5), Inches(6.05), Inches(4.7), LGREY)
rect(s, Inches(0.5), Inches(1.5), Inches(6.05), Inches(0.62), ACCENT)
txt(s, Inches(0.7), Inches(1.5), Inches(5.6), Inches(0.62), "PAPER 1 — the discovery paper (now)", size=16, color=WHITE, bold=True, anchor=MSO_ANCHOR.MIDDLE)
bullets(s, Inches(0.75), Inches(2.25), Inches(5.55), Inches(3.8), [
    "Title thrust: reproducible psoriasis signature converging on JAK-STAT/STAT3",
    "Core figures: cohort table, meta-summary + I², STAT3 forest, HKSJ/LOO robustness, TF-activity, clust module",
    "Selling points: multi-cohort generalisability, per-gene reproducibility label, retired isoform switch (integrity)",
    "Target: strong disease-genomics / dermatology journal (e.g. JID, JACI, Genome Medicine tier)",
    "Status: analysis complete — writing-ready",
], size=13.5, gap=9, bcolor=ACCENT)
# Paper 2 box
rect(s, Inches(6.8), Inches(1.5), Inches(6.05), Inches(4.7), LGREY)
rect(s, Inches(6.8), Inches(1.5), Inches(6.05), Inches(0.62), TEAL)
txt(s, Inches(7.0), Inches(1.5), Inches(5.6), Inches(0.62), "PAPER 2 — druggability follow-up (later)", size=16, color=WHITE, bold=True, anchor=MSO_ANCHOR.MIDDLE)
bullets(s, Inches(7.05), Inches(2.25), Inches(5.55), Inches(3.8), [
    "Signature → druggable target triage → structure-based screen",
    "RORγt demonstration screen; STAT3 as degrader/MD target",
    "Needs: full-scale docking + co-folding + extended MD on GPU cluster (packaged, ready to run)",
    "Optional wet-lab validation of top hits",
    "Status: preliminary / proof-of-concept done",
], size=13.5, gap=9, bcolor=TEAL)
takeaway(s, "Recommendation: submit the discovery paper first; carry druggability as a distinct follow-up. Open for lab discussion today.")
notes(s, "The decision slide for the lab. Per the agreed framing, Paper 1 is the meta-analysis/STAT3 discovery story — analysis is complete and writing-ready. Paper 2 is the druggability/structure arm, currently proof-of-concept with cluster-ready packages for the full runs. Discussion prompts: (1) target journal for Paper 1? (2) is STAT3 the headline or the reproducible-signature framing? (3) who owns the cluster runs for Paper 2? (4) do we want wet-lab validation before Paper 2? Adjust journal names to lab norms.")

# ---------- SLIDE 38: Discussion / next steps ----------
s = add_slide()
rect(s, 0,0, SW, SH, NAVY)
rect(s, Inches(0.9), Inches(1.35), Inches(2.6), Pt(5), ACCENT)
txt(s, Inches(0.9), Inches(0.7), Inches(11), Inches(0.6), "DISCUSSION", size=16, color=RGBColor(0xBF,0xD3,0xE8), bold=True)
txt(s, Inches(0.9), Inches(1.55), Inches(11.5), Inches(0.8), "Decisions for the lab", size=30, color=WHITE, bold=True)
bullets(s, Inches(0.95), Inches(2.7), Inches(11.4), Inches(4.2), [
    "Paper scope: confirm discovery-first, druggability as follow-up",
    "Headline framing: STAT3 lead vs. reproducible-signature-first — which sells better?",
    "Target journal & format (word/figure limits shape final figure set)",
    "Any additional cohorts to fold in? (treatment/baseline arms were deferred — ~230 more samples available)",
    "Robustness disclosure: how prominently to feature the HKSJ magnitude caveat (recommend: prominently)",
    "Who drives Paper 2 cluster runs (docking / co-folding / MD packages are ready)",
], size=17, gap=13, color=WHITE, bcolor=ACCENT)
notes(s, "Close on decisions, not data. Drive the lab toward concrete choices: scope, framing, journal, cohort expansion, how to disclose the magnitude caveat, and ownership of the follow-up compute. The deferred treatment/baseline arms (ERP110816, SRP065812; ~230 samples) are a ready lever if a reviewer or the lab wants more power.")

# ---------- SLIDE 39: Methods appendix ----------
s = add_slide()
header(s, "APPENDIX", "Methods provenance")
bullets(s, Inches(0.55), Inches(1.5), Inches(12.2), Inches(5.0), [
    "Data: recount3 whole-skin psoriasis RNA-seq, uniform Gencode v26 annotation; 3 well-powered cohorts (n=178, 93, 12) of 5 screened",
    "Preprocessing: per-study filterByExpr + TMM normalization; voom precision weights; limma moderated linear models",
    "Contrasts: PP–NN (lesional vs healthy), PN–NN, PN–PP (healthy → peri-lesional → lesional gradient)",
    "Meta-analysis: random-effects DerSimonian–Laird on per-study log2FC ± SE; vectorised, validated against metafor (agree to 4 dp)",
    "Heterogeneity: Cochran's Q, I², τ²; HKSJ (REML + t-interval) as recommended primary in the k=3 regime; leave-one-out sensitivity",
    "Pathways: GSEA + ORA (concordance); CAMERA + ROAST (inter-gene correlation); cross-study combination via Stouffer signed weighted-Z",
    "Regulators: decoupleR TF-activity inference on the CollecTRI network; STAT3 regulon and activity meta-analysed",
    "Co-expression: clust consensus modules across the 3 cohorts (numpy<2 env), restricted to meta-significant genes",
    "Splicing: junction-level PSI-β from recount3 RSE, depth floor ≥20; single-study vs pooled replication test",
    "Druggability (follow-up): Open Targets tractability + ChEMBL potent inhibitors (pChEMBL≥7); fpocket; AutoDock Vina; Boltz-2; OpenMM",
], size=14, gap=8)
notes(s, "Reference slide for method questions. Keep in the deck as backup; pull specific rows up if a reviewer-style question comes up in lab. All underlying tables, RDS objects, and figures are saved as project artifacts.")
print("total slides:", len(prs.slides._sldIdLst))
prs.save("psoriasis_meta_analysis_lab_presentation.pptx")
import os
print("saved", os.path.getsize("psoriasis_meta_analysis_lab_presentation.pptx"), "bytes")
