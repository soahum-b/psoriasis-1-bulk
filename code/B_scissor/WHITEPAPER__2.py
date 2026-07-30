# Auto-extracted generating script
# Produces: WHITEPAPER.pdf
# Conda env: python   (run with this environment activated)
# Inputs (expected alongside / in data/): fig_ref_composition.png, fig_ref_umap.png, fig_alpha_tuning.png, fig_scissor_umap.png, fig_scissor_composition.png, fig_significance.png, fig_gradient_program.png, fig_stat3.png, fig_deconv_validation.png
# Source artifact version: 9a926824-182f-4872-b282-75ecfd221b9c
# NOTE: absolute artifact paths rewritten to basenames; place named inputs in the same dir.
#------------------------------------------------------------

import os
import base64
import weasyprint

os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = os.path.join(os.environ.get("CONDA_PREFIX", "/Users/soahum/.claude-science/conda/envs/python"), "lib")

figs = {
    "composition": ("fig_ref_composition.png", "Single-cell reference composition (GSE173706, Ma et al. 2023). (A) Cells per biopsy tier after QC. (B) Cells per donor, stacked by tier; 11 donors contribute paired PN+PP biopsies."),
    "ref_umap": ("fig_ref_umap.png", "Reference UMAP (89,058 cells). (A) Nine broad lineages by canonical markers. (B) Cells colored by tier — peri-lesional (PN) cells occupy an intermediate position between normal (NN) and lesional (PP)."),
    "alpha": ("fig_alpha_tuning.png", "Scissor alpha tuning. Selected fraction across the graph-smoothing grid; alpha=0.40 was chosen (14.84% selected, under the 20% cutoff)."),
    "scissor_umap": ("fig_scissor_umap.png", "Scissor selection on the reference UMAP (20k subset). Scissor+ (gradient-tracking, lesional-associated), Scissor- (normal-associated), and background cells."),
    "scissor_comp": ("fig_scissor_composition.png", "Cell-type and tier composition of the Scissor-selected fractions. Scissor+ peaks at the PN/peri-lesional tier."),
    "signif": ("fig_significance.png", "Significance controls. (A) Reliability test: real CV-MSE vs 100 label permutations. (B) Selection permutation null: real pos-minus-neg tier gap vs 30 shuffled-label reruns."),
    "program": ("fig_gradient_program.png", "Gradient-tracking gene program: volcano of Scissor+ vs background differential expression (1,861 genes at padj<0.05)."),
    "stat3": ("fig_stat3.png", "STAT3 expression across Scissor classes: significantly elevated in Scissor+ (log2FC 0.43, padj 0.018)."),
    "deconv": ("fig_deconv_validation.png", "Orthogonal bulk deconvolution vs Scissor direction, by cell type. Endothelial (lesional-tracking), Fibroblast and Melanocyte (normal-tracking) are concordant; NK and Keratinocyte discordant; DC non-significant."),
}

fig_b64 = {}
for k, (path, cap) in figs.items():
    with open(path, "rb") as fh:
        fig_b64[k] = base64.b64encode(fh.read()).decode()

def img(k):
    vid, cap = figs[k]
    return f'''<figure><img src="data:image/png;base64,{fig_b64[k]}" alt="{k}"/>
<figcaption><b>Figure.</b> {cap}</figcaption></figure>'''

css = """
:root{--ink:#1a1a1a;--muted:#555;--rule:#d5d5d5;--accent:#8a1c30;--blue:#4C9BD4;--soft:#f7f6f4;}
*{box-sizing:border-box}
body{font-family:Georgia,'Times New Roman',serif;color:var(--ink);line-height:1.55;
 max-width:820px;margin:0 auto;padding:56px 40px 80px;background:#fff;font-size:16px}
h1{font-size:27px;line-height:1.25;margin:0 0 6px;font-weight:700}
.sub{font-size:16px;color:var(--muted);font-style:italic;margin:0 0 4px}
.meta{font-size:13px;color:var(--muted);margin:14px 0 0}
h2{font-size:19px;margin:34px 0 8px;padding-bottom:5px;border-bottom:2px solid var(--accent);font-weight:700}
h3{font-size:15.5px;margin:22px 0 6px;color:#2a2a2a;font-weight:700}
p{margin:9px 0}
.abstract{background:var(--soft);border-left:3px solid var(--accent);padding:16px 20px;margin:22px 0;font-size:15px}
.abstract h2{border:0;margin:0 0 8px;font-size:15px;letter-spacing:.06em;text-transform:uppercase;color:var(--accent)}
figure{margin:20px 0;text-align:center}
figure img{max-width:100%;border:1px solid var(--rule);border-radius:3px}
figcaption{font-size:12.5px;color:var(--muted);text-align:left;margin-top:7px;line-height:1.45}
table{border-collapse:collapse;width:100%;margin:16px 0;font-size:13.5px;font-family:Helvetica,Arial,sans-serif}
th,td{border:1px solid var(--rule);padding:6px 9px;text-align:left}
th{background:var(--soft);font-weight:700}
td.n{text-align:right;font-variant-numeric:tabular-nums}
code{font-family:'SF Mono',Consolas,monospace;font-size:13px;background:var(--soft);padding:1px 5px;border-radius:3px}
.key{background:#fbfaf7;border:1px solid var(--rule);border-radius:4px;padding:2px 16px;margin:16px 0}
.foot{margin-top:40px;padding-top:14px;border-top:1px solid var(--rule);font-size:12.5px;color:var(--muted)}
ul,ol{margin:8px 0 8px 4px;padding-left:22px}
li{margin:4px 0}
.tag{display:inline-block;background:var(--accent);color:#fff;font-family:Helvetica,Arial,sans-serif;
 font-size:10.5px;letter-spacing:.08em;padding:2px 8px;border-radius:3px;vertical-align:middle}
@media print{body{padding:0 12px}h2{page-break-after:avoid}figure{page-break-inside:avoid}}
"""

html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Scissor-on-gradient: peri-lesional psoriasis</title><style>{css}</style></head><body>

<h1>Cell states that track the psoriasis progression gradient:<br>a Scissor analysis anchored on an ordinal biopsy-site phenotype</h1>
<p class="sub">A gradient-aware application of phenotype-to-cell mapping in normal &rarr; peri-lesional &rarr; lesional skin</p>
<p class="meta"><span class="tag">WHITE PAPER &middot; BACKBONE</span> &nbsp; Working draft &middot; project <code>psoriasis-1-bulk</code> &middot; single-cell reference GSE173706 (Ma et al. 2023) &middot; bulk anchor SRP165679 (Tsoi et al. 2019)</p>

<div class="abstract"><h2>Abstract</h2>
<p>Psoriatic skin is usually studied as a two-state contrast (lesional vs. normal), which discards the clinically meaningful <i>peri-lesional</i> margin where disease is actively expanding. We reframe the problem as an <b>ordinal gradient</b> &mdash; uninvolved (NN) &lt; peri-lesional (PN) &lt; lesional (PP) &mdash; and ask which single cells co-vary with that gradient. Using Scissor, we project an ordinal bulk RNA-seq phenotype (93 biopsies spanning all three tiers, F=254 for the tier axis on bulk PC1) onto a 89,058-cell single-cell reference that contains the peri-lesional compartment. Because the phenotype label is a clinical biopsy site rather than a molecular signature derived from the same cells, the mapping is <b>non-circular by construction</b>. Selected cells are monotonic on the gradient (mean tier: Scissor&minus; 0.79 &lt; background 1.38 &lt; Scissor+ 1.43) and the gradient-tracking (Scissor+) fraction peaks at the peri-lesional tier. Both a reliability test (p&nbsp;=&nbsp;0.000) and a selection permutation null (p&nbsp;=&nbsp;0.000) confirm the signal is driven by real phenotype structure. Endothelial cells dominate the gradient-tracking population (5.2&times; enriched, OR&nbsp;11.3), and the associated 1,861-gene program is vascular-led; STAT3 is a significant, if modest, member (log2FC&nbsp;0.43, padj&nbsp;0.018). An orthogonal bulk deconvolution independently confirms the compositional trend for the three strongest lineages. This document reports the validated laptop-scale <b>backbone</b>; the full-census run and additional robustness tiers are staged for cluster execution.</p></div>

<h2>1&nbsp;&nbsp;Background and rationale</h2>
<p>Standard psoriasis transcriptomics contrasts lesional (PP) against uninvolved (NN) skin. That design is blind to the <b>peri-lesional</b> zone (PN) &mdash; the advancing margin a few millimetres outside the visible plaque, where the transition from health to disease is presumably underway. If a distinct cellular program initiates psoriatic conversion, the peri-lesional compartment is where it should be visible, and a two-state design cannot see it.</p>
<p>We therefore treat biopsy site as an <b>ordinal phenotype</b> (NN&nbsp;&lt;&nbsp;PN&nbsp;&lt;&nbsp;PP) and use <b>Scissor</b> (Sun, Guan, &hellip; Xia, <i>Nat Biotechnol</i> 2022) to identify the single cells whose expression co-varies with that gradient. Scissor correlates every single cell against a bulk cohort carrying a phenotype label, then solves a network-regularized regression to select the cells most consistent with the phenotype. Our central design choice is that the phenotype is a <b>clinical biopsy-site label</b>, not a molecular score computed from the reference cells &mdash; so a cell being flagged cannot be a restatement of how it was labelled. STAT3, a longstanding candidate in psoriasis, is examined as a member of the resulting program rather than as an assumed driver.</p>

<h2>2&nbsp;&nbsp;Data and design</h2>
<h3>2.1&nbsp;&nbsp;Single-cell reference (GSE173706)</h3>
<p>We assembled the Ma et al. (2023, <i>Nat Commun</i>) psoriasis atlas from 33 per-sample count matrices: <b>96,088 cells &times; 33,538 genes</b> across 22 donors, spanning all three tiers (NN&nbsp;13,534; PN&nbsp;35,518; PP&nbsp;47,036), with <b>11 donors contributing paired PN+PP biopsies</b> &mdash; the design feature that makes the peri-lesional compartment usable.</p>
{img("composition")}
<p>After Ensembl&rarr;symbol mapping (24,185 symbols) and QC (nFeature 200&ndash;6000, nCount&nbsp;&ge;500, percent-MT&nbsp;&lt;20), <b>89,058 cells (92.7%)</b> were retained. Standard normalization, HVG selection (2,000), PCA (30), neighbour graph (dims&nbsp;1:20), Louvain clustering (resolution&nbsp;0.5) and UMAP produced 23 clusters, annotated to nine broad lineages by canonical markers.</p>
{img("ref_umap")}

<h3>2.2&nbsp;&nbsp;Bulk phenotype anchor (SRP165679)</h3>
<p>Among five recount3 psoriasis studies carrying tier labels, only <b>SRP165679 (Tsoi et al. 2019)</b> balances all three tiers at usable depth: <b>NN=38, PN=27, PP=28</b> (93 classified samples). Using a single study removes cross-study confounding. edgeR logCPM (filterByExpr) and Ensembl&rarr;symbol collapse gave 24,533 symbols; the ordinal response is y&nbsp;&isin;&nbsp;{{0,1,2}}. Non-degeneracy checks pass: minimum tier size 27, and bulk PC1 (27.8% variance) tracks the tier axis at <b>F=253.9, p&nbsp;&asymp;&nbsp;10<sup>&minus;38</sup></b>.</p>
<table><tr><th>Check</th><th>Value</th></tr>
<tr><td>Tiers / min samples per tier</td><td class="n">3 / 27</td></tr>
<tr><td>Bulk PC1 variance explained</td><td class="n">27.8%</td></tr>
<tr><td>Bulk PC1 &sim; tier</td><td class="n">F=253.9, p&asymp;10<sup>&minus;38</sup></td></tr>
<tr><td>Cross-study confound</td><td>none (single study)</td></tr>
<tr><td>Phenotype circularity</td><td>none (clinical biopsy-site label)</td></tr></table>

<h2>3&nbsp;&nbsp;Method note: a portable Scissor solver</h2>
<p>Scissor's selection step is a graph-regularized elastic net solved by a compiled routine (APML1). That component could not be built in our environment, so we reimplemented the identical objective in pure R on top of <code>glmnet</code>. The network penalty is the symmetric-normalized graph Laplacian, <code>L_sym = I &minus; D<sup>&minus;1/2</sup> A D<sup>&minus;1/2</sup></code>; at ~10<sup>5</sup> cells a dense Laplacian is intractable, so we encode it as a <b>sparse edge-difference augmentation</b> &mdash; one sparse row per graph edge realizing &beta;<sup>T</sup>L<sub>sym</sub>&beta; = &Sigma;<sub>edges</sub>(&beta;<sub>i</sub>/&radic;d<sub>i</sub> &minus; &beta;<sub>j</sub>/&radic;d<sub>j</sub>)<sup>2</sup>. Sparsity is controlled by walking the glmnet &lambda;-path to a target selected fraction (the operational form of Scissor's cutoff). On synthetic two-community data the port recovered 100/100 true-positive cells (cor(&beta;,&beta;<sub>true</sub>)=0.72). This is a validated equivalent, <b>not</b> the authors' canonical solver; a cross-check against compiled Scissor is staged for the cluster.</p>

<h2>4&nbsp;&nbsp;Results</h2>
<h3>4.1&nbsp;&nbsp;Selection and directionality</h3>
<p>For tractable graph regression the backbone runs on a <b>stratified 20,023-cell subset</b> (preserving cell-type&times;tier proportions; full census staged for cluster). Correlating 93 bulk samples against 20,023 cells over 1,664 shared HVGs and tuning the graph-smoothing parameter gave <b>alpha=0.40, 14.84% of cells selected</b> (1,574 Scissor+, 1,397 Scissor&minus;), under the 20% cutoff.</p>
{img("alpha")}
<p>The selected cells are <b>monotonic on the gradient</b>: mean tier (NN=0, PN=1, PP=2) rises Scissor&minus; <b>0.79</b> &lt; background <b>1.38</b> &lt; Scissor+ <b>1.43</b>. Critically, the gradient-tracking Scissor+ fraction <b>peaks at the peri-lesional tier</b> (9.4% of PN cells vs. 8.0% of PP), the direct signature of an intermediate progression state that a two-state design would miss.</p>
{img("scissor_umap")}
{img("scissor_comp")}

<h3>4.2&nbsp;&nbsp;Significance</h3>
<p>Two independent controls confirm the result. The <b>reliability test</b> (100 label permutations) gives real cross-validated MSE <b>0.147</b> vs. a null mean of <b>0.779</b> (<b>p&nbsp;=&nbsp;0.000</b>, 0/100 permutations lower). The <b>selection permutation null</b> (30 reruns with shuffled tier labels, full re-selection) collapses the positive-minus-negative tier gap from the real <b>0.638</b> to a null mean of <b>&minus;0.305</b> (<b>p&nbsp;=&nbsp;0.000</b>). The selection reflects real phenotype structure, not graph geometry or overfitting.</p>
{img("signif")}

<h3>4.3&nbsp;&nbsp;Which cells track the gradient</h3>
<p>Endothelial cells overwhelmingly dominate the gradient-tracking fraction, with fibroblasts and melanocytes marking the opposite (normal-associated) pole.</p>
<table><tr><th>Cell type</th><th>Fold in Scissor+</th><th>Odds ratio</th><th>p</th><th>Interpretation</th></tr>
<tr><td>Endothelial</td><td class="n">5.18&times;</td><td class="n">11.26</td><td class="n">6&times;10<sup>&minus;245</sup></td><td>lesional-tracking</td></tr>
<tr><td>DC</td><td class="n">1.40&times;</td><td class="n">1.47</td><td class="n">0.014</td><td>weakly lesional</td></tr>
<tr><td>Keratinocyte</td><td class="n">0.87&times;</td><td class="n">0.69</td><td class="n">8&times;10<sup>&minus;12</sup></td><td>near-uniform</td></tr>
<tr><td>NK</td><td class="n">0.40&times;</td><td class="n">0.36</td><td class="n">2&times;10<sup>&minus;12</sup></td><td>depleted from Scissor+</td></tr>
<tr><td>Fibroblast</td><td class="n">0.18&times;</td><td class="n">0.15</td><td class="n">5&times;10<sup>&minus;53</sup></td><td>normal-tracking</td></tr>
<tr><td>Melanocyte</td><td class="n">0.04&times;</td><td class="n">0.04</td><td class="n">3&times;10<sup>&minus;20</sup></td><td>normal-tracking</td></tr></table>

<h3>4.4&nbsp;&nbsp;The gradient-tracking gene program</h3>
<p>Differential expression of Scissor+ vs. background (Wilcoxon) yields <b>1,861 genes at padj&nbsp;&lt;&nbsp;0.05</b>. The top up-regulated genes are vascular/endothelial &mdash; <b>CCL14, ACKR1, RAMP3, PLVAP, APLNR, CYTL1, SPNS2</b> &mdash; consistent with dermal angiogenesis as a hallmark of psoriatic progression.</p>
{img("program")}
<p><b>STAT3</b> is a significant member of the program (<b>log2FC 0.43, padj 0.018</b>; 48.5% vs. 44.7% of cells expressing). The direction matches the STAT3 hypothesis, but the effect is modest and the program is not STAT3-led &mdash; the dominant axis is vascular.</p>
{img("stat3")}

<h3>4.5&nbsp;&nbsp;Orthogonal validation by bulk deconvolution</h3>
<p>As an independent check by a different method, we deconvolved the real SRP165679 bulk (NNLS against a Ma-reference signature) and tested each cell type for a monotonic proportion trend across NN&rarr;PN&rarr;PP. Scissor per-cell direction is taken from the Fisher enrichment among Scissor+ cells; a cell type is concordant only when its bulk proportion trend is significant and matches.</p>
<table><tr><th>Cell type</th><th>Scissor direction</th><th>Bulk proportion trend</th><th>padj</th><th>Status</th></tr>
<tr><td>Endothelial</td><td>lesional-tracking</td><td>rises 0.0&rarr;0.2&rarr;3.8%</td><td class="n">1&times;10<sup>&minus;19</sup></td><td>&#10003; concordant</td></tr>
<tr><td>Fibroblast</td><td>normal-tracking</td><td>falls 12.3&rarr;8.2&rarr;1.4%</td><td class="n">1&times;10<sup>&minus;11</sup></td><td>&#10003; concordant</td></tr>
<tr><td>Melanocyte</td><td>normal-tracking</td><td>falls 7.8&rarr;7.3&rarr;2.2%</td><td class="n">1&times;10<sup>&minus;14</sup></td><td>&#10003; concordant</td></tr>
<tr><td>NK</td><td>normal-tracking</td><td>rises to 7.9% (PP)</td><td class="n">2&times;10<sup>&minus;11</sup></td><td>&#10007; discordant</td></tr>
<tr><td>Keratinocyte</td><td>normal-tracking</td><td>rises 79&rarr;83&rarr;83%</td><td class="n">5&times;10<sup>&minus;3</sup></td><td>&#10007; discordant</td></tr>
<tr><td>DC</td><td>lesional-tracking</td><td>flat</td><td class="n">0.20</td><td>trend n.s.</td></tr></table>
<p>The three strongest-signal lineages are concordant: the endothelial gradient signal reflects a <b>real compositional increase</b> in bulk tissue, not solely a state change. The discordant cases are informative rather than failures &mdash; NK cells <i>infiltrate</i> lesional tissue (bulk proportion rises) yet individual NK cells are depleted from the per-cell gradient-tracking set, and keratinocytes dominate every tier so per-cell direction and bulk fraction are not expected to align. Per-cell tracking and bulk composition are genuinely different measurements; they converge where the biology is a clean compositional shift and diverge otherwise.</p>
{img("deconv")}

<h2>5&nbsp;&nbsp;Discussion</h2>
<p>Reframing psoriasis as an ordinal gradient and mapping it with a non-circular phenotype yields a coherent, reproducible signal: a gradient-tracking cell population that is monotonic on the NN&rarr;PN&rarr;PP axis, peaks at the peri-lesional margin, survives two orthogonal significance controls, and is independently corroborated by bulk deconvolution. The dominant biology is <b>vascular</b> &mdash; endothelial expansion and an angiogenic gene program &mdash; which fits the histology of psoriatic progression and points attention to the dermal vasculature at the advancing edge, not only the epidermal keratinocyte compartment that two-state designs emphasize. STAT3 participates in the program in the expected direction but is not its driver here; its modest single-cell effect suggests any STAT3 role is embedded in a broader vascular/inflammatory circuit rather than acting alone.</p>

<h2>6&nbsp;&nbsp;Limitations</h2>
<ul>
<li><b>Backbone scale.</b> Results are on a 20,023-cell stratified subset; the full 89,058-cell census is staged for cluster execution and may shift point estimates (directionality expected to hold).</li>
<li><b>Solver.</b> The network elastic net is a validated pure-R reimplementation, not the authors' compiled APML1; a cross-check against stock Scissor is a planned step.</li>
<li><b>Single bulk anchor.</b> Tier-1 uses SRP165679 alone; robustness across additional tier-labelled cohorts (Tier-2/3) is future work.</li>
<li><b>Deconvolution method.</b> The validation uses transparent NNLS, not a benchmarked method choice; a benchmarking layer (deconvBenchmarking) would justify the specific deconvolution method for this tissue.</li>
</ul>

<h2>7&nbsp;&nbsp;Next steps</h2>
<ol>
<li><b>Full-census run</b> (89,058 cells) on the cluster; re-estimate all point values.</li>
<li><b>Compiled-Scissor cross-check</b> to confirm the glmnet port against the canonical solver.</li>
<li><b>Benchmarked deconvolution</b> (deconvBenchmarking) to select and re-run the best method for skin.</li>
<li><b>Robustness tiers</b> &mdash; additional bulk anchors and alternative tier definitions.</li>
<li><b>Sequence-level arm</b> (Sei-LLRA / seillra) &mdash; score regulatory variants near gradient-program genes; a distinct data type (DNA sequence) forming a downstream project, with STAT3 at the seam between the expression and sequence arms.</li>
</ol>

<h2>8&nbsp;&nbsp;Data and code availability</h2>
<p>Single-cell reference: <b>GSE173706</b> (Ma et al. 2023). Bulk anchor: <b>SRP165679</b> (Tsoi et al. 2019, via recount3). All code, figures, result tables, the environment specification, and this white paper are in the repository <code>github.com/soahum-b/psoriasis-1-bulk</code>; the full-census cluster script and reproduce/scale instructions are in <code>HANDOFF.md</code>. Heavy inputs are regenerated by <code>code/00_download_data.R</code>.</p>
<p><b>Method reference.</b> Scissor: Sun D, Guan X, Moran AE, <i>et al.</i> Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. <i>Nature Biotechnology</i> <b>40</b>, 527&ndash;538 (2022; published online 11 Nov 2021). doi:10.1038/s41587-021-01091-3.</p>

<p class="foot">Working draft &mdash; internal white paper for the peri-lesional psoriasis Scissor project. All quantitative values are computed from the saved analysis artifacts. Clinical interpretation is preliminary and not intended to guide patient care.</p>
</body></html>"""

weasyprint.HTML(string=html).write_pdf("WHITEPAPER.pdf")
print("PDF:", os.path.getsize("WHITEPAPER.pdf") // 1024, "KB")