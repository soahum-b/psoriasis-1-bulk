# Pathway & Gene-Set Enrichment Analysis: A Deep-Dive Review

*Companion to the psoriasis meta-analysis white paper (§7). Written to be read on its own — it
assumes you have a ranked list of genes and a differential-expression result, and builds up every
concept needed to interpret pathway enrichment, defend the method choices to a reviewer, and
answer questions about them.*

---

## 0. The problem pathway analysis has to solve

Differential expression (DE) hands you a list — in our psoriasis data, **3,477 significant genes**
out of 24,528 tested. A list is not an explanation. Three problems make a raw gene list hard to
interpret:

1. **Too long to read.** No one can hold 3,477 genes in their head, and the eye fixates on the few
   with the biggest fold-changes (DEFB4A at +13) while ignoring the coordinated behaviour of
   hundreds of moderate movers.
2. **No mechanism.** A list does not tell you *which biological processes* are engaged. "STAT3 is
   up" is a fact; "the IL-6 → JAK → STAT3 → IL-17 signalling circuit is coordinately activated" is
   an explanation.
3. **Not reproducible across studies.** The exact membership of a DE list is noisy — swap a few
   samples and genes near the threshold flip in and out. **Pathways are far more stable**, which is
   exactly why they matter for a meta-analysis.

**Pathway analysis** (also called gene-set analysis, functional enrichment, or over-representation
depending on the flavour) solves all three by asking a different question: not *"which genes
changed?"* but *"which predefined groups of functionally related genes changed together?"* The
groups come from curated databases (§6). The statistics come in two philosophies, which is the
central distinction of this whole document.

---

## 1. Two philosophies: threshold-based vs threshold-free

There are two fundamentally different ways to ask whether a gene set is "enriched," and they
answer subtly different questions. Understanding the difference is the single most important idea
in pathway analysis.

| | **Over-Representation Analysis (ORA)** | **Gene Set Enrichment Analysis (GSEA)** |
|---|---|---|
| Input | A **list** of DE genes (thresholded) | A **ranking** of *all* genes (no threshold) |
| Question | Are set members over-represented **among my hits**? | Do set members cluster at the **extremes** of the ranking? |
| Test | Fisher's exact / hypergeometric on a 2×2 table | Running enrichment score + permutation |
| Needs a cutoff? | **Yes** — you must define "DE" | **No** — uses every gene |
| Blind to | Everything below the threshold | Nothing |
| Classic tools | DAVID, clusterProfiler `enrichGO`/`enrichKEGG` | Broad GSEA, `fgsea` |

Neither is "correct" — they are complementary lenses, and the most defensible analysis (ours)
runs **both**. When they agree, the result is robust. When they disagree, the disagreement is
itself a finding about *how* the pathway is regulated (§5). The rest of this document builds up
each method in full, then shows the contrast on our own IL-17 data.

---

## 2. Over-Representation Analysis (ORA) in detail

### 2.1 The 2×2 contingency table

ORA starts by dividing the world of genes along two axes and counting. Fix a functional category
**X** (a pathway). Fix your DE threshold (ours: BH-FDR < 0.05 **and** |log2FC| > 1 → 3,477 genes).
Every gene in the **universe** (the 24,528 genes we actually tested) falls into exactly one of
four boxes:

|                       | **DE**      | **not DE**   | row total |
|-----------------------|-------------|--------------|-----------|
| **in pathway X**      | a           | b            | a + b     |
| **not in pathway X**  | c           | d            | c + d     |
| **column total**      | a + c       | b + d        | N         |

- **a** = genes that are both in the pathway and differentially expressed (the overlap you care about)
- **a + b** = the size of the pathway (restricted to tested genes)
- **a + c** = the total number of DE genes (3,477)
- **N** = the universe size (24,528)

The question "is X over-represented among my DE genes?" becomes: **is `a` bigger than you'd expect
if DE genes were sprinkled at random across the universe?** Under randomness, the DE rate inside
the pathway (a / (a+b)) should match the background DE rate ((a+c) / N). Enrichment means the
inside rate is meaningfully higher.

### 2.2 Fisher's exact test and the hypergeometric distribution

To turn "bigger than expected" into a p-value, ORA uses **Fisher's exact test** [1]. Its engine is
the **hypergeometric distribution** — the exact probability model for *sampling without
replacement*. The mental picture: an urn holds N genes, of which (a+b) are "in the pathway"
(coloured balls). You draw (a+c) genes without replacement (your DE hits). The probability that
your draw contains **exactly** `a` pathway genes is:

$$P(X = a) = \frac{\binom{a+b}{a}\binom{c+d}{c}}{\binom{N}{a+c}}$$

Fisher's exact test (one-sided, "greater") sums this probability for the observed `a` **and every
more-extreme value** (a+1, a+2, …) to give the p-value: *the probability of seeing an overlap this
large or larger by chance alone.* Because it is exact (not an approximation), it is valid even for
small pathways where chi-squared tests would misbehave. `clusterProfiler`'s `enrichGO`,
`enrichKEGG`, and `enricher` all use this hypergeometric/Fisher machinery under the hood [2].

### 2.3 The background ("universe") is a real choice, not a technicality

The single most common way to get ORA wrong is to pick the wrong **universe** (the N, the
denominator). The temptation is to use "all genes in the genome" (~20,000). But we did not test all
genes — filtering removed lowly expressed ones (§4 of the white paper), leaving 24,528. **The
correct universe is the set of genes that had a chance to be called DE** — i.e., the tested genes.
Using an inflated genome-wide background makes almost everything look falsely enriched, because it
pretends genes that were never testable could have been hits. We use the **24,528 tested genes** as
the universe throughout. This choice is not cosmetic: it can change which pathways pass, and a
reviewer will ask.

### 2.4 Worked example on our data — a cautionary pair

Two pathways, same analysis, opposite outcomes — this is the ORA half of the §5 contrast.

**GO Keratinization** (the psoriatic-plaque output):

|                | DE  | not DE |
|----------------|-----|--------|
| in set         | 46  | 33     |
| not in set     | 3431| 21018  |

DE rate inside = 46/79 = **58%** vs background 3,477/24,528 = **14.2%** → **4.1× enrichment,
Fisher p = 1.0e-19.** Emphatically over-represented. Here ORA and GSEA will agree.

**Reactome IL-17 signaling** (the same pathway GSEA flags strongly):

|                | DE  | not DE |
|----------------|-----|--------|
| in set         | 7   | 63     |
| not in set     | 3470| 20988  |

DE rate inside = 7/70 = **10%** vs 14.2% background → **0.71× (depleted), Fisher p = 0.88.** ORA
says *not enriched.* Only the 7 highest-fold members (IL17A, IL17F, IL17C, IL17RE, NOD2, IRAK2,
MAPKAPK3) survived the |log2FC| > 1 cut; the other 63 signalling-machinery genes were discarded.
GSEA, which keeps them, tells a different story (§3, §5).

**Lesson:** ORA's verdict is only as good as its threshold. A pathway driven by many moderate,
coordinated shifts is invisible to it.

---

## 3. Gene Set Enrichment Analysis (GSEA) in detail

GSEA [3] was invented precisely to fix ORA's blind spot. It never thresholds. Instead it asks:
**do the members of a gene set tend to sit near the top (or bottom) of a ranking of *all* genes?**
Four ingredients.

### 3.1 Ingredient 1 — the ranking metric (which statistic?)

Every gene gets a score, and the genes are sorted from most-up to most-down. The choice of score
matters:

- **log2 fold-change** alone: rewards big effects but ignores reproducibility. A gene that jumps
  +2 in one noisy sample ranks above a rock-solid +1.5.
- **the moderated t-statistic** (our choice): t = (log2FC) / (standard error). It rewards genes
  that are **both large in effect and consistent across samples**. This is the limma empirical-Bayes
  t-statistic from the DE step, so the ranking inherits the variance-shrinkage that makes limma
  robust (§3 of the DE deep-dive). We store it as `gsea_ranks.rds` (24,528 genes, pre-sorted).
- **signed −log10(p)** or the **sign(logFC)·statistic**: common alternatives; all put "big and
  reliable" at the ends and "small or noisy" in the middle.

The key property: the metric must be **signed** (up-genes positive, down-genes negative) so the two
tails of the ranking mean opposite biology. Using an unsigned metric (e.g. |t|) would fold up- and
down-regulation together and destroy directional interpretation.

### 3.2 Ingredient 2 — the walk down the ranked list

Start at the top of the ranking (most up-in-psoriasis) with a running score of 0. Walk downward,
one gene at a time:

- **Hit** (the gene is in the set): the running score steps **up**. The step is *weighted by the
  magnitude of the ranking statistic* — a set member near the very top pushes harder than one in
  the middle. (This is the "weighted" / p=1 Kolmogorov–Smirnov-style enrichment statistic;
  classic GSEA lets you tune the weight exponent, but weighted is standard.)
- **Miss** (not in the set): the running score steps **down** by a small, uniform amount
  (1 / number-of-non-members).

If the set's genes are concentrated at the top, the up-steps come early and cluster, so the running
score **climbs to a peak before the down-steps drag it back**. If the genes are scattered randomly,
ups and downs interleave and the score wanders near zero the whole way.

### 3.3 Ingredient 3 — the Enrichment Score (ES) and its normalization (NES)

- **ES** = the maximum deviation of the running score from zero (the height of the peak). A tall
  positive peak = enrichment at the top; a deep negative trough = enrichment at the bottom.
- **ES depends on set size** — bigger sets can reach higher peaks by chance — so it is not
  comparable across sets as-is. GSEA normalizes it by the mean ES of the permutation null (§4) to
  give the **Normalized Enrichment Score (NES)**, which *is* comparable across sets and collections.
- **Sign of NES is the whole biology:**
  - **NES > 0** → set clusters at the **top** → coordinately **UP** in lesional psoriasis
    (all our STAT3/IL-17/NF-κB sets).
  - **NES < 0** → set clusters at the **bottom** → coordinately **DOWN** (e.g. the
    prostate-development sets that came out negative).

In our data the Hallmark IL6/JAK/STAT3 set peaks at ES = +0.635 near rank 4,761 of 24,528, which
normalizes to **NES = +2.49** (white-paper Fig 8).

### 3.4 Ingredient 4 — the leading edge

The **leading-edge subset** is the group of set members that appear *before* the peak — i.e., the
genes actually responsible for the enrichment signal. These are the biological "core" of the
result and the genes worth reporting. For IL6/JAK/STAT3 the leading edge includes **STAT3 itself**,
STAT1, SOCS3 (its direct negative-feedback target), OSMR, MYD88, IRF1/9, and the interferon /
interleukin receptor chains — a coherent, interpretable core. Leading-edge genes for every set are
saved in `gsea_results_all.csv`.

---

## 4. The p-value: what GSEA permutes, and why it matters

An ES (or NES) is meaningless until you know what a *random* gene set of the same size would score.
GSEA builds that null distribution by **permutation**, and there are two schemes — the difference
is a favourite reviewer question.

### 4.1 Sample (phenotype) permutation — the statistical gold standard

The **original 2005 GSEA** [3] shuffles the **sample labels** (which samples are PP vs NN),
recomputes the *entire per-gene ranking* on the shuffled labels, and recomputes the ES. Repeat
1,000× to get a null distribution of ES for that set.

- **Advantage:** it preserves the **gene–gene correlation structure** — real genes are co-regulated,
  and permuting samples keeps co-expressed genes moving together, so the null is biologically
  honest. This is the statistically correct null.
- **Cost:** you must recompute the ranking for every permutation, which is expensive, and you need
  enough samples per group for label-shuffling to generate a rich null. With our **83 NN / 95 PP**
  we have ample samples to do this for a headline pathway if desired.

### 4.2 Gene-set (gene-rank) permutation — the fast modern default

`fgsea` [4] instead **permutes gene ranks** (equivalently, samples random gene sets of the same
size from the ranking) while holding the observed ranking fixed. It never recomputes the ranking,
so it is orders of magnitude faster, and its **multilevel Monte-Carlo scheme** can estimate the very
small p-values we see (down to ~1e-43 for E2F targets) that ordinary permutation could never reach
with a feasible number of shuffles.

- **Advantage:** speed, and access to tiny accurate p-values; enables testing thousands of sets
  across four collections in seconds.
- **Caveat:** it treats genes as independent, so it **ignores gene–gene correlation** and is
  slightly **anti-conservative** for highly co-regulated sets (real p may be a touch larger than
  reported). This is the accepted, standard trade-off, and it is what we use.

**Rule of thumb:** gene-rank permutation (`fgsea`) for the genome-wide screen across all
collections; sample-label permutation for the one or two headline pathways where you want a
bullet-proof p-value. We do the former; the latter is available given our sample size.

### 4.3 Multiple testing across gene sets

Whichever permutation scheme, you are now testing many sets (50 Hallmark, 1,839 Reactome, 7,538
GO:BP, 186 KEGG), so the **same multiple-testing logic from the DE step applies again** — one level
up. `fgsea` returns a **BH-adjusted `padj`** across all sets in a collection [5]. We report `padj <
0.05`. (The full FWER/FDR/q-value treatment is in the DE deep-dive, §5; nothing new conceptually —
the unit being corrected is now "gene set" instead of "gene.")

---

## 5. When GSEA and ORA disagree — the IL-17 case, in full

This is the payoff of running both methods, worked on our own data.

| Pathway | ORA (Fisher) | GSEA (fgsea) | Verdict |
|---|---|---|---|
| GO Keratinization | 4.1×, p = 1e-19 ✅ | NES +2.54, padj 8.5e-10 ✅ | **Agree** — robust |
| Reactome IL-17 signaling | 0.71×, p = 0.88 ❌ | NES +2.02, padj 0.006 ✅ | **Disagree** — informative |

**Why keratinization agrees:** its genes move *dramatically* — 46 of 79 clear |log2FC| > 1. Whether
you threshold (ORA) or rank (GSEA), the signal survives. Agreement across two different statistical
philosophies is the strongest possible evidence: the finding is not an artefact of method choice.

**Why IL-17 disagrees:** look at *which* genes changed. The set has 70 tested members. The 7 that
clear the strict fold-change bar are the **cytokines and receptors** — IL17A (+5.7), IL17F (+4.7),
IL17C, IL17RE — huge effects. The other ~63 are the **intracellular signal-transduction machinery**
— IRAK1/4, TRAF6, MAP kinases, adaptors — which shift **coordinately but modestly** (perhaps +0.4
to +0.9), not enough to pass |log2FC| > 1. 

- **ORA discards those 63 at the threshold.** It sees only 7 hits out of 70, computes a DE rate
  (10%) *below* background (14.2%), and declares the pathway *depleted*.
- **GSEA keeps all 70.** It sees that essentially the entire set sits in the upper half of the
  ranking — the cytokines at the very top, the machinery pushed consistently upward — and correctly
  reports strong coordinated activation.

**The disagreement is the finding.** It tells you IL-17 signalling in psoriasis is driven by a few
high-amplitude ligands riding on top of a broadly, moderately up-shifted signalling apparatus —
exactly the "coordinated moderate movement" GSEA was built to detect and ORA is structurally blind
to. Reporting both, and explaining the discrepancy, is more informative than either alone. This is
why the analysis is not "GSEA vs ORA" but "GSEA **and** ORA."

### 5.1 The disagreement is systematic: ORA is consistently the more conservative test

The IL-17 case is not a one-off. When we ran ORA across **all** 6,505 gene sets and compared
verdict-by-verdict with GSEA, ORA called far fewer sets significant in *every* collection:

| Collection | GSEA significant | ORA significant | ratio |
|---|---|---|---|
| Hallmark | 31 | 14 | 2.2× |
| Reactome | 545 | 42 | 13× |
| GO:BP | 1,419 | 272 | 5.2× |
| KEGG | 67 | 17 | 3.9× |

Restricting to the 3,584 up-enriched sets both methods evaluate, they agree on **77%**, but the
disagreement is almost entirely one-directional: **1,204 sets are GSEA-significant yet ORA-null,
versus only 64 the reverse** (176 significant by both). ORA is effectively a conservative *subset*
of GSEA — nearly everything ORA flags, GSEA also flags, plus a large tier of coordinated-moderate
programs that never clear ORA's hard threshold. This is the structural consequence of §5's
mechanism, generalized: any pathway whose signal lives in many small, consistent shifts (rather than
a few large ones) is invisible to a threshold test and visible to a rank test. **Neither is "right"
— ORA answers "is this pathway over-represented among my strongest hits?" and GSEA answers "is this
pathway coordinately shifted across the whole transcriptome?" The conservatism is the price ORA pays
for only ever looking at the tail of the distribution.**

### 5.2 "Did we cherry-pick the pathways?" — the unbiased landscape

A fair worry: the dot-plots (Figs 9, main paper) show a *curated* STAT3-axis set. Did we manufacture
the STAT3 story by choosing which pathways to show? The honest test is to rank **every** significant
pathway by the data alone — no manual selection — and see what floats to the top.

![The 30 strongest pathways, ranked purely by NES across all 6,505 tested sets. Nothing manually selected; themes are labeled post-hoc from pathway names. Proliferation/cell-cycle (18/30) and interferon (3/30) dominate the extreme top; the STAT3/IL-17/NF-kB axis appears lower — strong but not the largest-magnitude signal.]({{artifact:art_febc0145-a1fd-416e-9d34-20950a89944c}})

The result is honest and instructive. Of the 1,490 significant up-regulated sets, the 30 with the
**largest** NES are dominated by **proliferation / cell-cycle** (18 of 30 — E2F/G2M/MYC targets, DNA
replication, mitotic checkpoints) and **interferon response** (3 of 30). The STAT3 axis does *not*
top this list: `IL6_JAK_STAT3` ranks **99th**, `TNFA_via_NFKB` 198th, IL-17 signaling ~560th — all
strongly significant (NES 1.8–2.5, padj ≤ 1e-3), but not the largest-magnitude signals in the
transcriptome.

**Does this demote STAT3? No — and understanding why is the key conceptual point.** NES magnitude
measures *how coordinately and strongly a set moved*, not *how causally central it is*. The
cell-cycle sets have enormous NES because keratinocyte **hyperproliferation** is the most florid
downstream *consequence* in a plaque — thousands of cell-cycle genes move in lockstep. But that
proliferation is **driven by** the upstream cytokine circuit: IL-23 → IL-17 → IL-6 → JAK → **STAT3**
→ (SOCS3, cyclins, anti-apoptotic program) → keratinocyte proliferation. STAT3 is the **hub that
causes** the top-ranked signature; it is not itself the largest-amplitude signature. This is a
general lesson in pathway analysis: **the strongest-NES pathway is usually the phenotypic output;
the causal driver sits upstream with a smaller but highly significant footprint.** Ranking by
magnitude finds the *symptom*; mapping the network finds the *cause*. Both matter, and conflating
them is a classic misreading of GSEA output.

So the answer to "should we let the data pick the pathways?" is **yes, and we did** — the analysis
ran genome-wide on all 24,528 genes and all 6,505 sets with no pre-selection (§0, §3). Manual
curation entered only at the *figure* stage, to tell one coherent mechanistic story out of ~1,500
significant sets. The unbiased landscape above is the check that the story was not manufactured: the
STAT3 axis is robustly present, and its position *below* the proliferation signature is exactly what
a driver-versus-output model predicts.

**One caveat that keeps the unbiased ranking honest: the databases are themselves biased.** Counting
"how many significant sets" or "what has the largest NES" is not a neutral census of biology — it is
partly a census of **what has been annotated**. Cell-cycle, DNA-replication, and core-metabolism
processes have been studied for decades and are catalogued in exhaustive, finely-subdivided detail:
Reactome alone splits mitosis into dozens of overlapping sub-sets (checkpoints, sister-chromatid
separation, spindle assembly, …), so a single biological signal — "keratinocytes are dividing" —
gets counted many times over and floods the top of any magnitude ranking. Tissue-specific programs
like keratinocyte terminal differentiation, or the precise wiring of IL-17 signal transduction, are
annotated far more sparsely. So Fig 11's dominance by proliferation reflects **both** real biology
(plaques really are hyperproliferative) **and** an annotation echo (proliferation is over-catalogued
relative to skin-specific immunology). The practical rules this implies: (i) never read "more
significant sets" as "more important biology"; (ii) treat heavily-overlapping sets from one database
as one finding, not many; (iii) weigh a mechanistically specific set (e.g. `IL6_JAK_STAT3`) on par
with a generic high-magnitude one even when its NES is smaller. This is a further reason the curated,
mechanism-driven figure (Fig 9) and the unbiased landscape (Fig 11) are **complementary** — the
first corrects for annotation bias by reasoning about mechanism; the second guards against
cherry-picking. Neither alone is sufficient.

---

## 6. Choosing gene-set databases — each answers a different question

We tested four MSigDB [6] collections. They are not interchangeable; each is a different resolution
and vocabulary.

- **Hallmark (H)** — 50 refined, deliberately **non-redundant** "signatures," each distilled from
  many overlapping sets to reduce redundancy [6]. **Best first look**: minimal overlap, easy to
  read, one clean row per biological theme. This is where *IL6/JAK/STAT3 signaling* and *TNFα
  signaling via NF-κB* live. Start here.
- **KEGG** [7] — classic hand-curated pathway **diagrams** (the wiring maps in textbooks). Best for
  a **canonical named-pathway** story (*JAK-STAT signaling pathway*). Fewer, broader, well known to
  clinicians.
- **Reactome** [8] — expert-curated, **finely granular and hierarchical** (reactions nested inside
  pathways). This is where the **specific mechanistic sets** live: *STAT3 nuclear events*, *IL-17
  signaling*, *formation of the cornified envelope*. Best for drilling into mechanism.
- **GO Biological Process** [9] — the **broadest and most granular** (7,538 sets), but highly
  **redundant** (parent/child terms overlap heavily). Best for **discovery** and for
  biological-process language (*keratinization*, *cornification*, *keratinocyte differentiation*).

**Standard strategy (ours):** lead with **Hallmark** for the clean headline, drill into
**Reactome/KEGG** for named-pathway mechanism, use **GO:BP** for breadth and discovery. Reporting the
same signal across independent collections (as our STAT3 axis does) is itself a robustness argument —
it is not an artefact of one database's idiosyncratic curation.

A note on **redundancy**: GO:BP and Reactome will return many overlapping significant sets (e.g. a
dozen flavours of "NF-κB signalling"). This is expected, not a bug. Tools exist to collapse
redundancy (e.g. clusterProfiler's `simplify()` for GO, or leading-edge similarity clustering /
enrichment maps); for a focused hypothesis like ours, curating a representative set by hand (as in
white-paper Fig 9) is clearer than an automated collapse.

---

## 7. Why the dot-plot is the right summary figure

The **dot-plot** (white-paper Fig 9) encodes four dimensions at once, which is why it is the
standard GSEA centrepiece:

- **x-position = NES** — strength *and direction* (right = up in psoriasis).
- **y-order** — pathways sorted by NES, so the eye reads a ranked landscape top-to-bottom.
- **dot size = gene-set size** — context for how broad each signal is.
- **dot colour = −log10(padj)** — significance, so the deep-red high-NES corner is "strong and
  certain."

One frame communicates "which pathways, how strongly, in which direction, how significantly" — far
more legible than four separate bar charts. The alternative single-pathway figure, the **enrichment
"mountain plot"** (Fig 8), is the right choice when you want to *show the mechanism of the ES
itself* for one headline pathway — the running score, the leading-edge ticks, the peak. Use the
mountain plot to teach/prove one pathway; use the dot-plot to summarize many.

---

## 8. One-paragraph summary defense

*We interpret the 3,477-gene differential-expression result at the pathway level using two
complementary methods. Our primary method is threshold-free **GSEA** (`fgsea`, multilevel
permutation) on all 24,528 genes ranked by the limma moderated t-statistic, tested against four
MSigDB collections (Hallmark, Reactome, KEGG, GO:BP); positive NES denotes coordinate
up-regulation in lesional skin, with BH correction across sets. We cross-check with **ORA**
(Fisher's exact test on 2×2 tables, using the 24,528 tested genes as the background universe) on
the DE hit list. Where the methods agree (keratinization: 4.1× ORA enrichment, GSEA NES +2.54) the
signal is robust to method choice; where they disagree (IL-17 signalling: GSEA-significant but
ORA-null) the discrepancy localizes the regulation to coordinated moderate shifts in the signalling
machinery that a fold-change threshold discards. Across all four collections the IL6→JAK→STAT3→
IL-17/NF-κB→keratinocyte program is coordinately activated, corroborating the STAT3 hypothesis at
the pathway level.*

---

## 9. Beyond competitive GSEA: self-contained tests, correlation correction, and per-sample scores

Everything up to §8 used **competitive** tests: GSEA and ORA both ask *"are the genes in this set
more changed than the rest of the genome?"* Our professor's reading list points to three ideas that
sit outside that frame, and we implemented all three on our data. They are worth understanding
because each fixes a specific blind spot in plain GSEA.

### 9.1 Two null hypotheses, not one — the distinction that reorganizes the whole field

Goeman & Bühlmann [11] made the cleanest statement of a distinction that is easy to miss: a gene-set
test can ask one of **two entirely different questions**, and they have different null hypotheses,
different validity conditions, and different failure modes.

- **Competitive** null: *"the genes in set S are no more associated with the phenotype than the
  genes outside S."* The comparison is **set vs the rest of the genome**. GSEA, ORA, and CAMERA are
  competitive. The unit of replication is the **gene** — so the p-value depends on how many other
  genes exist and how they behave.
- **Self-contained** null: *"no gene in set S is associated with the phenotype."* There is **no
  reference to genes outside S at all**. ROAST is self-contained. The unit of replication is the
  **sample** — so the p-value is bounded by how many samples you have, and it answers "is this set
  changed *at all*?" rather than "is it changed *more than average*?"

Why this matters for us: a self-contained test is the natural question for a **named hypothesis set**
like the STAT3 regulon. We are not really asking "are STAT3 targets more deranged than a random 371
genes?" — we already believe the whole transcriptome is deranged in psoriasis. We are asking "is the
STAT3 target program engaged?" That is the self-contained question, and ROAST answers it directly.

### 9.2 The inter-gene correlation problem — why the plain GSEA p-value is over-confident

The fast, modern GSEA default (gene-rank permutation, §4.2) has a known statistical flaw: it treats
the genes in a set as **independent** draws. They are not. Genes in a pathway are **co-regulated** —
when STAT3 turns on, dozens of its targets move together. Positively correlated genes carry
**redundant** evidence, so counting each as an independent vote inflates the effective sample size
and shrinks the p-value far below what the data support. Wu & Smyth [12] quantified this: inter-gene
correlation of even 0.05 can inflate a gene-set test statistic several-fold.

**CAMERA** [12] is competitive GSEA's correction for exactly this. It estimates a **variance
inflation factor (VIF)** from the residual inter-gene correlation and deflates the test statistic
accordingly. The formula is essentially the classic "design effect" from survey sampling: the
effective number of independent genes is `n / (1 + (n−1)·ρ̄)`, not `n`.

We ran CAMERA on our voom fit (`~group`, NN vs PP contrast, symbol-matched gene sets). The result is
a lesson in calibration. On the **STAT3 regulon** (371 measured targets), the three methods give:

| Test | Type | p-value | What it corrects / ignores |
|---|---|---|---|
| GSEA (`fgsea`) | competitive | **7.1×10⁻¹⁴** | treats 371 targets as independent |
| CAMERA | competitive, correlation-corrected | **4.6×10⁻⁴** | discounts co-regulation (VIF) |
| ROAST | self-contained | **5.0×10⁻⁵** | asks "engaged at all?", rotation-validated |

All three agree the regulon is **up in lesional skin** — the conclusion is robust. But the p-values
span **ten orders of magnitude**, and the gap between GSEA and CAMERA *is* the inter-gene correlation
Goeman & Bühlmann and Wu & Smyth warn about. **For the white paper we quote CAMERA's p, not GSEA's**:
it is the honest competitive statement once co-regulation is accounted for. (White-paper Fig 18.)

Across the 50 Hallmark sets, GSEA and CAMERA **agree on direction and significance for 27/50 sets**
(GSEA-significant: 31; CAMERA-significant: 27) — the correction removes a handful of the least robust
calls but leaves the core immune/proliferation program (E2F, interferon-α/γ, G2M, MYC,
IL6/JAK/STAT3, inflammatory response) intact and up.

### 9.3 ROAST — the self-contained rotation test

ROAST [13] implements the self-contained null by **rotation** rather than permutation. Permutation
needs many samples to generate a null distribution by shuffling phenotype labels; rotation generates
the null by rotating the residual space, so it works with **few samples and complex designs** (it can
handle covariates, which sample-label permutation cannot). It reports directional components
separately. On the STAT3 regulon our ROAST result was **Up p = 5.0×10⁻⁵**, **Mixed p = 1×10⁻⁴**
(the "Mixed" alternative catches sets where genes move in both directions), with **Down p = 1.0** —
i.e. the set is unambiguously engaged in the up direction, exactly as the biology predicts. For a
named hypothesis set this is arguably the most appropriate single test we run. `mroast` extends it to
many sets at once with FDR control.

### 9.4 GSVA / ssGSEA — turning a pathway into a per-sample score

Every method so far collapses the whole cohort into **one** p-value per gene set. **GSVA** [14] and
its close relative **ssGSEA** [15] do something different and complementary: they compute a **score
for each pathway in each sample**, turning a gene × sample expression matrix into a much smaller
pathway × sample matrix. That unlocks everything you can do with a sample-level number — correlate it
with a clinical variable, cluster samples on it, feed it to a downstream model, or pool it across
cohorts in a meta-analysis.

The algorithm (ssGSEA, Barbie et al. [15]): within each sample, rank all genes by expression; walk
the ranked list accumulating a weighted Kolmogorov-Smirnov-like statistic that rises on set members
and falls on non-members; the integral of that walk is the sample's enrichment score. It is GSEA's
running-sum idea (§3.2) applied **one sample at a time** with no phenotype labels.

Because neither the Bioconductor `GSVA` package nor its compiled dependencies would install in our
sandboxed environment, **we implemented ssGSEA directly in ~20 lines of base R** (rank → weighted
ECDF difference → integral → range-normalize). This is a feature, not a workaround: the method is
fully transparent and auditable rather than a black box. We validated it against the known biology —
the `HALLMARK_IL6_JAK_STAT3_SIGNALING` per-sample score separates the groups cleanly (**NN mean
0.505 vs PP mean 0.636, Wilcoxon p = 5.9×10⁻²⁸**), reproducing the cohort-level GSEA verdict from a
completely independent per-sample computation.

Two things the per-sample score buys us that a single p-value cannot (white-paper Figs 19–20):

1. **It tracks STAT3 expression continuously.** The per-sample IL6/JAK/STAT3 score correlates with
   STAT3 log-CPM at **Spearman ρ = 0.81** across all 178 samples — the pathway activity is not a
   binary on/off between groups but a graded axis, and STAT3 sits on it.
2. **It reconnects to the isoform story.** Correlating the pathway score against **PSI_β** (the
   STAT3-β fraction) reproduces the same **Simpson's-paradox** signature we found with EFTUD2:
   overall ρ = 0.09 (p = 0.24), but **within lesional skin ρ = −0.18 (p = 0.074)**. Higher pathway
   activity within PP goes with *less* β (more activating α) — consistent with, not contradicting,
   the model that the β brake is regulated separately from bulk pathway drive.

A pathway-score **heatmap** (Fig 20) orders samples NN→PP and shows the two-block structure at the
per-sample level: the metabolic/lipid programs (myogenesis, adipogenesis, fatty-acid metabolism,
OXPHOS) are high in normal skin, while the immune and proliferation programs (IL6/JAK/STAT3,
inflammatory response, interferon-α/γ, allograft rejection, E2F/G2M/MYC) are high in lesional — the
same biology as our differential-expression and network results, now visible sample by sample. This
per-sample matrix is exactly the object we will pool across studies in the meta-analysis.

### 9.5 Leading-edge redundancy — reading an UpSet plot

GSEA reports many significant immune pathways, and a naive reading counts them as many independent
discoveries. They are not. The **leading edge** (§3.4) of each pathway — the genes actually driving
its score — overlaps heavily across pathways, because psoriasis has essentially **one** dominant
immune program (IL-17 / interferon / NF-κB) that many curated sets each capture a slice of.

An **UpSet plot** [16] is the right way to show this. Where a Venn diagram collapses past 3–4 sets,
UpSet displays intersections as a matrix: each column is a specific combination of sets, its bar is
how many genes fall in exactly that combination. We took the leading edges of the six top immune
Hallmark pathways (inflammatory response, TNFα/NF-κB, IL6/JAK/STAT3, interferon-γ, allograft
rejection, complement) and plotted their overlap (white-paper Fig 21). The numbers make the
redundancy concrete: the six leading edges sum to **537 gene-memberships but only 373 unique genes
(redundancy factor 1.44)**, and a shared core — **IL6, IRF1, CXCL9/10/11, CCL2/5, IL1B, TLR2** —
recurs in four or more of the six. The practical lesson: **do not report six co-significant immune
pathways as six independent findings.** They are largely one signal, and the UpSet plot is the
honest way to say so.

### 9.6 A note on gene-set databases and benchmarking

Two further pointers from the reading list. First, choosing methods is not a matter of taste — it has
been **benchmarked**. Geistlinger et al. [17] built a systematic benchmark (the `GSEABenchmarkeR`
framework) across 42 datasets and found, among other things, that many methods are mis-calibrated
under the competitive null and that set-level tests accounting for correlation (like CAMERA) are
better behaved — empirical support for the choice we made in §9.2. Second, for an **immunology**
phenotype like psoriasis the MSigDB **C7 immunologic collection** [18] (curated from sorted-immune-
cell expression signatures) is a natural specialized vocabulary to layer on top of Hallmark, and the
**GO primer** of du Plessis et al. [19] is the readable entry point for how GO's structure (the DAG,
`is_a`/`part_of` relations, evidence codes, the true-path rule) shapes every GO-based result. Finally,
for the historical record: the GSEA method traces to **Mootha et al. [20]** (2003, *Nature Genetics*)
— the OXPHOS/diabetes paper that introduced the running-enrichment-score idea — before the canonical
2005 PNAS formalization [3].

---

## References

[1] Fisher RA. *Statistical Methods for Research Workers.* Oliver & Boyd, Edinburgh; 1925.
*(Fisher's exact test; hypergeometric tail probability for 2×2 tables.)*

[2] Yu G, Wang L-G, Han Y, He Q-Y. *clusterProfiler: an R package for comparing biological themes
among gene clusters.* OMICS. 2012;16(5):284–287. doi:10.1089/omi.2011.0118. *(hypergeometric ORA
implementation: enrichGO/enrichKEGG/enricher.)*

[3] Subramanian A, Tamayo P, Mootha VK, et al. *Gene set enrichment analysis: a knowledge-based
approach for interpreting genome-wide expression profiles.* Proc Natl Acad Sci USA.
2005;102(43):15545–15550. doi:10.1073/pnas.0506580102. *(original GSEA: running ES, weighted KS
statistic, sample-label permutation, leading edge.)*

[4] Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A. *Fast gene set
enrichment analysis.* bioRxiv. 2021. doi:10.1101/060012. *(fgsea multilevel Monte-Carlo p-value
estimation for arbitrarily small p; gene-rank permutation.)*

[5] Benjamini Y, Hochberg Y. *Controlling the false discovery rate: a practical and powerful
approach to multiple testing.* J R Stat Soc B. 1995;57(1):289–300. *(BH adjustment, applied across
gene sets.)*

[6] Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. *The Molecular
Signatures Database (MSigDB) hallmark gene set collection.* Cell Syst. 2015;1(6):417–425.
doi:10.1016/j.cels.2015.12.004. *(MSigDB; Hallmark refinement to reduce redundancy.)*

[7] Kanehisa M, Goto S. *KEGG: Kyoto Encyclopedia of Genes and Genomes.* Nucleic Acids Res.
2000;28(1):27–30. doi:10.1093/nar/28.1.27.

[8] Milacic M, Beavers D, Conley P, et al. *The Reactome Pathway Knowledgebase 2024.* Nucleic Acids
Res. 2024;52(D1):D672–D678. doi:10.1093/nar/gkad1025.

[9] Ashburner M, Ball CA, Blake JA, et al. *Gene Ontology: tool for the unification of biology.*
Nat Genet. 2000;25(1):25–29. doi:10.1038/75556. *(GO Biological Process.)*

[10] Väremo L, Nielsen J, Nookaew I. *Enriching the gene set analysis of genome-wide data by
incorporating directionality of gene expression.* Nucleic Acids Res. 2013;41(8):4378–4391.
doi:10.1093/nar/gkt111. *(comparison of gene-set methods; directional enrichment.)*

[11] Goeman JJ, Bühlmann P. *Analyzing gene expression data in terms of gene sets: methodological
issues.* Bioinformatics. 2007;23(8):980–987. doi:10.1093/bioinformatics/btm051. *(the
self-contained vs competitive null distinction; choice of background/universe.)*

[12] Wu D, Smyth GK. *Camera: a competitive gene set test accounting for inter-gene correlation.*
Nucleic Acids Res. 2012;40(17):e133. doi:10.1093/nar/gks461. *(inter-gene correlation and its
effect on gene-set test calibration — the caveat behind gene-rank permutation.)*

[13] Wu D, Lim E, Vaillant F, Asselin-Labat M-L, Visvader JE, Smyth GK. *ROAST: rotation gene set
tests for complex microarray experiments.* Bioinformatics. 2010;26(17):2176–2182.
doi:10.1093/bioinformatics/btq401. *(self-contained rotation test; mroast; works with few samples
and covariate-adjusted designs.)*

[14] Hänzelmann S, Castelo R, Guinney J. *GSVA: gene set variation analysis for microarray and
RNA-seq data.* BMC Bioinformatics. 2013;14:7. doi:10.1186/1471-2105-14-7. *(per-sample,
unsupervised gene-set enrichment scores; sample-level pathway matrix.)*

[15] Barbie DA, Tamayo P, Boehm JS, et al. *Systematic RNA interference reveals that oncogenic
KRAS-driven cancers require TBK1.* Nature. 2009;462(7269):108–112. doi:10.1038/nature08460.
*(single-sample GSEA, ssGSEA — the per-sample projection we implemented in base R.)*

[16] Conway JR, Lex A, Gehlenborg N. *UpSetR: an R package for the visualization of intersecting
sets and their properties.* Bioinformatics. 2017;33(18):2938–2940.
doi:10.1093/bioinformatics/btx364. *(UpSet matrix visualization; scales past Venn diagrams. Method:
Lex A, Gehlenborg N, et al. IEEE Trans Vis Comput Graph. 2014;20(12):1983–1992.)*

[17] Geistlinger L, Csaba G, Santarelli M, et al. *Toward a gold standard for benchmarking gene set
enrichment analysis.* Brief Bioinform. 2021;22(1):545–556. doi:10.1093/bib/bbz158. *(systematic
benchmark of enrichment methods across 42 datasets; calibration and correlation caveats.)*

[18] Godec J, Tan Y, Liberzon A, et al. *Compendium of immune signatures identifies conserved and
species-specific biology in response to inflammation.* Immunity. 2016;44(1):194–206.
doi:10.1016/j.immuni.2015.12.006. *(MSigDB C7 immunologic signature collection.)*

[19] du Plessis L, Škunca N, Dessimoz C. *The what, where, how and why of gene ontology — a primer
for bioinformaticians.* Brief Bioinform. 2011;12(6):723–735. doi:10.1093/bib/bbr002. *(GO structure:
the DAG, is_a/part_of relations, evidence codes, the true-path rule.)*

[20] Mootha VK, Lindgren CM, Eriksson K-F, et al. *PGC-1α-responsive genes involved in oxidative
phosphorylation are coordinately downregulated in human diabetes.* Nat Genet. 2003;34(3):267–273.
doi:10.1038/ng1180. *(the original running-enrichment-score idea, predating the 2005 PNAS GSEA
formalization [3].)*

---

*Document status: complete. Companion to `psoriasis_meta_analysis_whitepaper.md` §7 and to
`normalization_deep_dive.md` and `differential_expression_deep_dive.md`. Covers ORA (Fisher/
hypergeometric, 2×2 tables, universe choice), GSEA (ranking metric, running ES, NES, leading edge),
the two permutation schemes, the GSEA-vs-ORA contrast on IL-17, database selection, and figure
choice. §9 adds the competitive-vs-self-contained distinction (Goeman & Bühlmann), CAMERA's
inter-gene-correlation correction, ROAST's self-contained rotation test, GSVA/ssGSEA per-sample
scores (implemented in base R), UpSet leading-edge redundancy, and the benchmarking/database
reading list — all applied to our data and cross-referenced to white-paper Figs 18–21.*
