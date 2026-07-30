# Network analysis: a deep dive

*Companion to the psoriasis meta-analysis white paper (§8). This document explains what "network
analysis" means in transcriptomics, why we reach for it, and the main families of methods — what
each one assumes, what it computes, what it can and cannot tell you, and when to use it. It is
written to be read on its own; the worked psoriasis result lives in the white paper.*

---

## 1. Why a network at all — the limitation a list cannot escape

Every method up to this point in the pipeline produces a **list**. Differential expression gives a
list of genes with effect sizes. Pathway analysis (GSEA/ORA) gives a list of gene *sets* with
enrichment scores. A list is an unordered bag: it tells you *what* changed, ranked by how much, but
it contains no arrows. It cannot tell you *what caused what*, because causation is a relationship
**between** genes and a list stores no relationships.

This produces a specific, recurring failure of interpretation, which we call the **driver-versus-
output problem**. In psoriasis the single loudest transcriptional signature is proliferation — cell
cycle, E2F targets, MYC targets, mitotic machinery. Rank pathways purely by magnitude and
proliferation dominates the top of the list. But a psoriatic plaque proliferates *because* it is
inflamed; proliferation is the **output** (the visible thickened skin), not the **driver** (the
inflammatory signalling that ordered the thickening). A list ranks by loudness and therefore finds
the symptom. To find the cause you need the arrows — you need to know that inflammatory
transcription factors sit *upstream* of the cell-cycle genes. That is what a network encodes.

A **network** (equivalently, a **graph**) is a set of **nodes** connected by **edges**. In biology
the nodes are usually genes, proteins, or metabolites; the edges are relationships — "regulates,"
"binds," "is co-expressed with," "is in the same pathway as." Edges may be:

- **Directed or undirected.** A→B ("A regulates B") is directed; A—B ("A and B interact") is
  undirected. Direction is what lets you talk about upstream and downstream.
- **Signed or unsigned.** A—|B (A represses B) versus A→B (A activates B). Sign is what lets a
  method predict *which way* a target should move.
- **Weighted or binary.** An edge can carry a confidence or strength, or simply be present/absent.

The whole art of network analysis is (i) choosing where the edges come from — prior knowledge versus
inferred from the data itself — and (ii) choosing what to compute on the resulting graph. The
sections below organise the major method families along exactly those two axes.

---

## 2. Two ways to get the edges: prior-knowledge vs data-driven

This is the single most important distinction, because it determines what a result *means*.

**Prior-knowledge (literature-curated) networks.** The edges come from decades of published
experiments — ChIP-seq, reporter assays, perturbation studies, protein-interaction screens —
aggregated into a database. Examples: CollecTRI and DoRothEA (TF→target regulons), STRING and
IntAct (protein–protein interactions), Reactome and KEGG (pathway/reaction graphs), OmniPath (an
aggregator of many such resources). You bring a database of edges and lay your data on top of it.
The strength: every edge has independent experimental support, so a result is interpretable and
mechanistic. The weakness: coverage is uneven and biased toward well-studied genes, and the edges
are *generic* — a TF→target edge curated in a cancer cell line may not hold in keratinocytes.

**Data-driven (inferred) networks.** The edges are *learned* from the expression data itself:
genes whose expression covaries across samples are connected. Examples: WGCNA (co-expression
modules), GENIE3 and ARACNe (regulatory-network reverse-engineering), Bayesian networks. You bring
no prior graph; the algorithm draws the edges. The strength: the network is specific to *your*
tissue and condition, and can discover relationships no database contains. The weakness: covariation
is not causation (two genes can correlate because a third drives both), edges have no mechanistic
guarantee, and you generally need many samples (tens to hundreds) for the correlations to be stable.

A useful mental model: **prior-knowledge methods test a map you already have; data-driven methods
draw a new map.** They answer different questions and are often used together — infer a module with
WGCNA, then annotate it against curated pathways.

---

## 3. Transcription-factor activity inference (the method used in §8)

**The question it answers:** *Which regulators are active, given the genome-wide expression
changes?* This is the most direct network answer to the driver-versus-output problem, which is why
it is the method we ran.

**The core idea — activity ≠ mRNA.** A transcription factor's own mRNA level is a poor readout of
whether it is *doing* anything. TFs are switched on post-transcriptionally: STAT3 is inactive until
it is phosphorylated, dimerizes, and enters the nucleus; NF-κB is held in the cytoplasm by IκB until
a signal releases it. In all these cases protein and mRNA can be flat while activity swings wildly.
But an active TF cannot hide its footprint: **the genes it activates go up and the genes it
represses go down.** So instead of reading the TF, you read its **targets** and infer the TF's
activity from their collective behaviour. This is sometimes called "footprint" or "regulon-based"
inference.

**The regulon.** A regulon is the curated target set of one TF, each target signed +1 (the TF
activates it) or −1 (the TF represses it). The two dominant human regulon resources:

- **DoRothEA** [3]: TF→target regulons with per-edge **confidence tiers A–E**, assembled from
  literature, ChIP-seq, TF-binding motifs, and inferred co-expression. The convention is to use the
  high-confidence tiers (A–C) for bulk data.
- **CollecTRI** [1]: a newer, larger signed resource (≈1,200 TFs) that expanded coverage while
  keeping the sign of regulation, and outperformed DoRothEA on benchmark perturbation data for
  recovering the perturbed TF. This is what §8 uses.

Both are distributed through **OmniPath** [7], an aggregator that harmonises dozens of prior-
knowledge resources behind one interface.

**The scoring — decoupleR and its models.** Given a regulon and a genome-wide statistic per gene
(we use the limma **t-statistic** from the differential-expression step — one number per gene
encoding direction and confidence of change), a method scores each TF by how coherently its signed
targets moved. **decoupleR** [2] is a framework that implements a family of such scorers behind one
API:

- **ULM (univariate linear model)** — the one we used. For a given TF, regress the per-gene
  statistic on that TF's signed regulon membership (target = +1/−1/0). The regression **slope's
  t-value** is the activity score: strongly positive means the activated targets went up and the
  repressed targets went down, in concert. Simple, fast, robust, and the recommended default for a
  single contrast.
- **MLM (multivariate linear model)** — regress on *all* TFs' memberships simultaneously, so a
  gene shared by two TFs is apportioned rather than double-counted. Better at disentangling TFs with
  overlapping regulons, at the cost of stability when regulons collide heavily.
- **consensus** — run several methods and combine their normalized scores, to avoid over-trusting
  any single model's assumptions. decoupleR's own benchmarking motivates this as a robust default
  when you are not sure which single method to trust [2].
- **ULM/MLM are the modern successors to earlier approaches** such as **VIPER** [4] (an
  enrichment-style regulon method built on the aREA algorithm and ARACNe-derived regulons), and
  simple **GSEA-on-regulons**. The ideas are continuous; decoupleR unifies them.

**What you get out.** A ranked table of TFs by activity score with a p-value (typically permutation-
or model-based, then multiple-testing corrected across TFs). The top of the list is your set of
candidate **master regulators** — the drivers. In §8 this ranking put NF-κB, STAT1, MYC, E2F1, and
STAT3 at the top, directly naming the inflammatory hubs that a list-based method could only imply.

**What it cannot tell you.** (i) Activity is *inferred*, not measured — there is no protein or
phospho data behind it; a footprint is strong evidence of activity but not a direct assay. (ii) The
regulon is generic and possibly incomplete for your tissue. (iii) A regulon with very few or very
one-sided edges (e.g. mostly activating, few repressing) can only report on the arm it covers — the
STAT3 repression-arm caveat in §8 is exactly this. (iv) It ranks regulators by target-coherence,
which correlates with but is not identical to biological importance.

---

## 4. Protein–protein interaction (PPI) and centrality networks

**The question it answers:** *In the network of physical/functional interactions, which proteins are
structurally central — hubs, bottlenecks — among my genes of interest?*

**The edges** come from interaction databases: **STRING** [5] (functional associations combining
experimental interactions, co-expression, text-mining, and conserved co-occurrence, each edge
scored), **IntAct** and **BioGRID** (curated physical interactions), Reactome's functional-
interaction network. Edges here are typically **undirected** (A binds B) and **unsigned** — you know
two proteins interact, not who activates whom.

**What you compute — centrality.** Once you have a graph you rank nodes by graph-topological
importance:

- **Degree** — how many partners a node has. High-degree nodes are **hubs**; hubs tend to be
  essential and are often disease genes.
- **Betweenness** — how often a node lies on the shortest path between other pairs. High-betweenness
  nodes are **bottlenecks** that control flow between modules.
- **Closeness, eigenvector/PageRank centrality** — variants weighting "how central" differently
  (average distance to all nodes; being connected to other central nodes).

A common workflow: take your DE gene list, pull the induced STRING subnetwork, and ask which genes
are hubs/bottlenecks — these are prioritised as candidate key players. Tools: **Cytoscape** with the
**MCODE** (dense-cluster detection) and **cytoHubba** (centrality ranking) plug-ins are the standard
interactive route.

**Strength:** intuitive, well-tooled, mechanistically grounded in physical biology. **Weakness:**
undirected/unsigned edges cannot give you up/downstream or driver/output; topological centrality is
a property of the *database's* wiring, which is study-bias-heavy (well-studied proteins look more
central simply because they are better annotated — the same annotation-bias caveat that dogs pathway
enrichment).

---

## 5. Co-expression networks (data-driven modules)

**The question it answers:** *Which genes move together across my samples, forming modules, and how
do those modules relate to phenotype?* This is the flagship **data-driven** family — the edges are
learned, not looked up.

**WGCNA** [6] (Weighted Gene Co-expression Network Analysis) is the canonical method:

1. Compute pairwise correlation between all gene pairs across samples.
2. Raise it to a **soft-threshold power** β so the network is approximately **scale-free** (a few
   hubs, many leaves) — this softly down-weights weak correlations instead of hard-cutting them.
3. Cluster genes into **modules** of tightly co-expressed genes (assigned colours).
4. Summarise each module by its **eigengene** (first principal component), and correlate eigengenes
   with traits (disease status, severity). Modules whose eigengene tracks the phenotype are the
   interesting ones.
5. Within a module, **intramodular connectivity** identifies **hub genes** — the most-connected,
   most-representative members.

**Strength:** entirely tissue-/condition-specific (it only knows *your* data), reduces thousands of
genes to a handful of interpretable modules, and connects those modules directly to clinical traits.
It is the workhorse for "I have many samples and want unsupervised structure." **Weakness:** needs
many samples (rule of thumb ≥ 15–20, ideally more) for stable correlations; co-expression is
undirected and non-causal (a module is a set of genes that covary, not a regulatory circuit); module
boundaries depend on parameter choices.

**Regulatory-network inference** is the directed cousin. **GENIE3** [8] reframes network inference
as: for each gene, use a random-forest to predict its expression from all other genes; the most
predictive regressors become its regulators — yielding a **directed, weighted** regulatory graph.
**ARACNe** uses mutual information plus a data-processing-inequality step to prune indirect edges.
These learn direction from data (GENIE3 won the DREAM5 network-inference challenge), but still do
not give **sign**, and demand large sample sizes.

---

## 6. Causal and Bayesian networks

**The question it answers:** *What is the most probable directed, causal wiring that generated my
data?* This is the most ambitious — and most assumption-heavy — family.

A **Bayesian network** represents genes as nodes in a **directed acyclic graph** (DAG) with a joint
probability distribution factorised over the edges; structure-learning algorithms search for the DAG
that best explains the observed expression. With **interventional** data (knockouts, perturbations)
the edges can carry genuinely causal meaning; with observational data alone, causal direction is
only partially identifiable and rests on strong assumptions (acyclicity, no hidden confounders,
faithfulness). Related frameworks: **Mendelian randomization** and other causal-inference methods
that use genetic variants as natural experiments to orient edges.

**Strength:** the only family that targets *causality* directly, and it produces a fully generative
model you can reason and simulate over. **Weakness:** computationally hard (DAG search is
super-exponential), data-hungry, and fragile to its assumptions — with pure observational bulk
RNA-seq the identifiability limits are severe. In practice these methods shine with perturbation
data or as a refinement step on a smaller candidate set, not as a first pass on a two-group contrast.

---

## 7. How to choose — a decision guide

| You have… and want to know… | Method family | Representative tools |
|---|---|---|
| A ranked DE statistic; *which regulators are active/driving* | **TF-activity inference** (§3) | decoupleR (ULM/MLM/consensus), VIPER; CollecTRI, DoRothEA regulons |
| A DE gene list; *which genes are structural hubs* | **PPI + centrality** (§4) | STRING, IntAct; Cytoscape (MCODE, cytoHubba) |
| Many samples; *unsupervised modules linked to phenotype* | **Co-expression** (§5) | WGCNA; GENIE3 / ARACNe for directed inference |
| Perturbation/large data; *causal wiring* | **Causal / Bayesian** (§6) | bnlearn, Bayesian structure learning |

Three practical rules that fall out of the above:

1. **Match the method to the question, not the fashion.** "Which driver?" → TF activity. "Which
   hub?" → centrality. "Which module?" → co-expression. Running all three and reporting whichever
   looks best is p-hacking with extra steps.
2. **Prior-knowledge methods work with few samples; data-driven methods need many.** TF-activity and
   PPI-centrality give a sensible answer from a single two-group contrast; WGCNA/GENIE3/Bayesian
   methods need tens-to-hundreds of samples before their edges are trustworthy.
3. **Triangulate.** The strongest claims come from *orthogonal* methods agreeing — as in §8, where
   gene-set overlap (pathway analysis) and regulon-target coherence (TF activity) are different
   statistics on different inputs and named the same STAT3/NF-κB/STAT1 circuit. Agreement across
   data models is worth more than a small p-value from any one of them.

---

## 8. The driver-versus-output concept, stated once more

The reason network analysis earns its place in this pipeline is a single idea worth stating plainly:
**the largest signal is usually the output, and the cause is usually smaller and upstream.** Rank
anything by magnitude — genes by fold-change, pathways by enrichment — and you find the phenotype's
loudest consequence. In psoriasis that is proliferation. To find the *cause* you must impose
directionality, which means imposing a network. TF-activity inference does this by asking not "what
changed most?" but "whose fingerprint explains the changes?", and the answer — inflammatory
regulators sitting upstream of the proliferative output — is the mechanistic hypothesis a list can
never yield. Ranking by magnitude finds the symptom; mapping the network finds the cause.

---

## References

[1] Müller-Dott S, Tsirvouli E, Vazquez M, Ramirez Flores RO, Badia-i-Mompel P, Fallegger R,
Türei D, Lægreid A, Saez-Rodriguez J. **Expanding the coverage of regulons from high-confidence
prior knowledge for accurate estimation of transcription factor activities.** *Nucleic Acids Res.*
2023;51(20):10934–10949. doi:10.1093/nar/gkad841. *(CollecTRI.)*

[2] Badia-i-Mompel P, Vélez Santiago J, Braunger J, Geiss C, Dimitrov D, Müller-Dott S, Taus P,
Dugourd A, Holland CH, Ramirez Flores RO, Saez-Rodriguez J. **decoupleR: ensemble of computational
methods to infer biological activities from omics data.** *Bioinformatics Advances.*
2022;2(1):vbac016. doi:10.1093/bioadv/vbac016.

[3] Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. **Benchmark and integration
of resources for the estimation of human transcription factor activities.** *Genome Res.*
2019;29(8):1363–1375. doi:10.1101/gr.240663.118. *(DoRothEA.)*

[4] Alvarez MJ, Shen Y, Giorgi FM, Lachmann A, Ding BB, Ye BH, Califano A. **Functional
characterization of somatic mutations in cancer using network-based inference of protein activity.**
*Nat Genet.* 2016;48(8):838–847. doi:10.1038/ng.3593. *(VIPER.)*

[5] Szklarczyk D, Kirsch R, Koutrouli M, Nastou K, Mehryary F, Hachilif R, Gable AL, Fang T,
Doncheva NT, Pyysalo S, Bork P, Jensen LJ, von Mering C. **The STRING database in 2023:
protein–protein association networks and functional enrichment analyses for any sequenced genome of
interest.** *Nucleic Acids Res.* 2023;51(D1):D638–D646. doi:10.1093/nar/gkac1000.

[6] Langfelder P, Horvath S. **WGCNA: an R package for weighted correlation network analysis.** *BMC
Bioinformatics.* 2008;9:559. doi:10.1186/1471-2105-9-559.

[7] Türei D, Korcsmáros T, Saez-Rodriguez J. **OmniPath: guidelines and gateway for
literature-curated signaling pathway resources.** *Nat Methods.* 2016;13(12):966–967.
doi:10.1038/nmeth.4077. *(See also Türei et al., Mol Syst Biol 2021, doi:10.15252/msb.20209923, for
the integrated OmniPath resource.)*

[8] Huynh-Thu VA, Irrthum A, Wehenkel L, Geurts P. **Inferring regulatory networks from expression
data using tree-based methods.** *PLoS ONE.* 2010;5(9):e12776. doi:10.1371/journal.pone.0012776.
*(GENIE3.)*

[9] Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, Califano A.
**ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular
context.** *BMC Bioinformatics.* 2006;7(Suppl 1):S7. doi:10.1186/1471-2105-7-S1-S7.

[10] Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T.
**Cytoscape: a software environment for integrated models of biomolecular interaction networks.**
*Genome Res.* 2003;13(11):2498–2504. doi:10.1101/gr.1239303.

---

*Companion documents: `psoriasis_meta_analysis_whitepaper.md` (§8 applies TF-activity inference to
SRP035988), `pathway_analysis_deep_dive.md` (the list-based layer this network layer builds on),
`differential_expression_deep_dive.md`, `normalization_deep_dive.md`.*
