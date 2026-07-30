# RNA-seq Normalization: A Deep-Dive Review

*Companion note to the psoriasis meta-analysis white paper. Purpose: a self-contained
reference you can study later to defend every normalization decision. Focused on the one
question that determines everything — **what are you trying to compare?***

---

## 0. The one idea that organizes the whole topic

There is no universally "correct" normalization. The right method depends entirely on the
comparison you intend to make. Two fundamentally different questions exist:

| The question you're asking | "Direction" of comparison | Correct method family |
|---|---|---|
| Is gene **A** more abundant than gene **B** *in the same sample*? | **within-sample** | RPKM, FPKM, **TPM** |
| Is gene **A** more abundant *in sample 1 than in sample 2* (or group 1 vs group 2)? | **between-sample** | CPM (depth only), **TMM**, **DESeq2 size factors** |

Our task — **differential expression (DE)** of lesional vs. normal psoriatic skin — is a
**between-sample** comparison. This single fact disqualifies the within-sample methods
(TPM/RPKM) as our DE normalization and points us to the composition-aware between-sample
methods. Everything below is an elaboration of this table.

The classic reviewer trap is using a within-sample unit (TPM/RPKM) for a between-sample
question (DE). We avoid it.

---

## 1. What normalization is actually correcting

Raw RNA-seq counts are corrupted by three nuisance factors that have nothing to do with
biology:

1. **Sequencing depth (library size).** Sample A sequenced to 40M reads will have roughly
   double the counts of sample B at 20M reads, gene-for-gene, purely from depth. Must be
   removed for *any* between-sample comparison.
2. **Gene length.** A 10-kb transcript captures more reads than a 1-kb transcript at the
   same molar abundance, simply because it presents more sequence to fragment and read.
   Matters *only* when comparing *different genes* to each other.
3. **Library composition.** The subtle, dangerous one — the subject of §5. Because a library
   is a fixed budget of reads, a few very highly expressed genes steal reads from all
   others, distorting between-sample comparisons even after depth correction.

A normalization method is essentially a choice about *which* of these three it corrects —
and, crucially, whether it corrects composition. That last property is what separates
DE-appropriate methods from the rest.

---

## 2. Within-sample methods: RPKM, FPKM, TPM

### 2.1 RPKM / FPKM

**RPKM** = Reads Per Kilobase of transcript per Million mapped reads.
**FPKM** = the paired-end equivalent (Fragments instead of Reads; one fragment = one or two
mates). They are the same idea.

**Computation (RPKM):** for each gene, `RPKM = (reads / gene_length_kb) / (total_reads / 1e6)`
— i.e. divide by depth, then by length.

**The defect (this is the citeable result).** Because RPKM divides by depth *before* fixing
length, **the sum of all RPKM values is not the same from sample to sample.** That means an
RPKM value of "20" does not represent the same fraction of the transcriptome in two different
samples — so RPKM is not properly comparable even for the within-sample purpose it was built
for. Wagner, Kin & Lynch (2012) demonstrated this inconsistency directly and argued RPKM
should be replaced by TPM [5]; Bullard et al. (2010) had already shown that total-count/RPKM
style scaling is overly sensitive to a handful of high-count genes [6].

### 2.2 TPM (the better within-sample unit)

**TPM** = Transcripts Per Million. Same two corrections (depth and length) but in the
**opposite order**: divide by length *first* (yielding a per-transcript rate), then rescale
so that **every sample's TPM values sum to exactly 1,000,000**.

**Why the reordering matters.** Because every sample sums to the same total, a TPM of 50
represents the *same fraction of the transcriptome* in every sample. TPM is therefore
internally consistent and is the preferred unit for "what fraction of this sample's
transcripts is gene X" and for comparing genes within a sample. TPM originates in the RSEM
framework of Li & Dewey (2011) [9] and is the community consensus replacement for RPKM
(Wagner 2012) [5].

### 2.3 Why we do **not** use TPM/RPKM for differential expression — two independent reasons

1. **Gene-length correction is unnecessary and adds noise for DE.** DE compares *the same
   gene* across samples. That gene's length is identical on both sides of the comparison, so
   it cancels algebraically. Dividing by length does nothing useful and injects
   length-estimate noise; worse, it lets a few long, highly expressed genes distort the
   per-sample scaling.
2. **The fixed-sum property reintroduces compositional bias (§5).** TPM forcing every sample
   to sum to 1,000,000 is exactly the compositional constraint that causes trouble: if the
   psoriatic keratins and S100 genes consume a large share of the million-unit budget in
   lesional samples, **every other gene's TPM is mechanically pushed down**, producing false
   "down-regulation." TPM *builds in* the artifact that TMM is designed to remove.

**Bottom line:** TPM is a **reporting/abundance unit**, not a DE normalization. It is correct
to quote TPM when describing how much of a sample a gene represents; it is incorrect to feed
TPM into a differential-expression test.

---

## 3. Between-sample, depth-only: CPM / RPM

**CPM** = Counts Per Million; **RPM** = Reads Per Million — the same quantity.
`CPM = count / (total_mapped_reads / 1e6)`. It corrects **sequencing depth only** — no length
term (so it never leaves the "same gene across samples" world) and **no composition term**.

**Its limitation.** CPM assumes the **total** count is a fair scaling factor — that libraries
differ only in depth. But the total is dominated by the most highly expressed genes. When a
few genes balloon in one condition (as in psoriasis), the total is inflated by *them*, and
dividing everything by that inflated total unfairly deflates every other gene. **CPM fixes
depth but is blind to composition.** Bullard et al. (2010) quantified exactly this
sensitivity of total-count scaling to a few high-count genes [6].

**Where CPM is nonetheless correct — and where we use it.** For **visualization, clustering,
and PCA**, depth adjustment is sufficient and log2-CPM is the standard scale. We use
**log2-CPM for Figures 1–4 and the PCA**, but **not** as the normalization inside the DE
model. The honest one-liner: **CPM for display, TMM for inference.**

---

## 4. Between-sample, composition-aware — the DE-correct family

These are the methods actually appropriate for our task. There are two well-validated
options, and they are *siblings*, not competitors.

### 4.1 TMM — Trimmed Mean of M-values (edgeR) — **our choice**

Pick a reference sample. For every gene compute its **M-value**, the log2 ratio of expression
in the test sample vs. the reference, `M = log2(test / reference)`. If there were no
composition bias, most genes should be unchanged (M ≈ 0). TMM then:

1. **Trims the extremes** — discards the top and bottom 30% of genes by M-value (likely
   truly DE) and the most extreme 5% by absolute expression, then
2. takes a **precision-weighted mean of the surviving M-values.**

That trimmed mean is the sample's normalization factor. **Load-bearing assumption:** most
genes are not differentially expressed, so the "typical" log-ratio should be zero; any
systematic departure is composition bias, which TMM scales away. Robust to a few dominant
genes precisely because it trims before averaging. Introduced by Robinson & Oshlack (2010)
[4].

### 4.2 DESeq2 median-of-ratios ("size factors") — the validated sibling

For each gene, form a **pseudo-reference** = the geometric mean of that gene's counts across
all samples. For each sample, take the ratio of every gene to its pseudo-reference; the
sample's **size factor** is the **median** of those ratios. Same robustness idea (a median
resists a few extreme genes) and the **same core assumption** (most genes unchanged).
Introduced for DESeq by Anders & Huber (2010) [7] and carried into DESeq2 by Love, Huber &
Anders (2014) [8].

### 4.3 TMM vs. median-of-ratios — how to defend "why TMM"

- They are the **two validated between-sample normalizations**, and independent benchmarks
  (Dillies et al. 2013 [10]) recommend **both** over total-count/RPKM/TPM approaches for DE,
  finding they give concordant, well-behaved results.
- Practical differences: median-of-ratios drops any gene with a **zero in any sample** (its
  geometric mean is zero), so with many samples/zeros it can discard more genes; TMM is a bit
  more tolerant there. TMM outputs scaling *factors* for library sizes (pairs with
  edgeR/limma-voom); DESeq2 outputs *size factors* (pairs with DESeq2's negative-binomial
  GLM).
- **Our justification:** we run the **edgeR / limma-voom** pipeline, whose native
  normalization is TMM. Choosing TMM keeps the normalization consistent with the test
  framework, and it is known to agree closely with DESeq2's median-of-ratios — so the choice
  is defensible as "the correct partner for our test," not "better than the alternative."

---

## 5. The concept underneath everything: RNA-seq is compositional

**The core fact.** A sequencing run delivers a **fixed budget** of reads. The counts within a
sample necessarily **sum to that library size**. Therefore counts are **not absolute
abundances** — they are **relative proportions of a fixed pie.** This is the definition of
**compositional data** (the same statistical structure as any set of proportions constrained
to sum to a constant).

**The consequence.** If one group massively over-expresses a few genes, those genes take a
larger slice of the pie, so **every other gene's slice mechanically shrinks — even genes
whose true absolute expression never changed.** Uncorrected, this shows up as spurious
**down-regulation** of unrelated genes.

**Why psoriasis is a textbook worst case.** Lesional skin over-expresses antimicrobial and
cornified-envelope genes (DEFB4, S100A7/8/9/12, SPRR2 family, LCE3) to extraordinary levels —
the exact "a few genes eat the pie" scenario. Compositional correction is therefore not
optional in this dataset; it is essential to avoid manufacturing false down-regulation across
the rest of the transcriptome.

**What corrects composition and what does not:**

| Method | Depth | Length | **Composition** | Verdict for our DE |
|---|:---:|:---:|:---:|---|
| RPKM / FPKM | ✓ | ✓ | ✗ (worsens) | ✗ within-sample unit |
| TPM | ✓ | ✓ | ✗ (fixed-sum builds it in) | ✗ within-sample unit |
| CPM / RPM | ✓ | ✗ | ✗ | ✗ for inference (✓ for display) |
| **TMM** (ours) | ✓ | ✗ | **✓** | **✓ DE-correct** |
| DESeq2 size factors | ✓ | ✗ | **✓** | ✓ DE-correct (sibling) |

---

## 6. The 30-second verbal defense

> *"We used **TMM**. RNA-seq is compositional — counts are proportions of a fixed read budget
> — and psoriasis is an extreme case where a few over-expressed skin genes would deflate
> everything else. TPM and RPKM are within-sample units, and TPM even bakes in the
> compositional constraint, so they are wrong for differential expression. Plain CPM corrects
> depth but not composition. TMM and DESeq2's median-of-ratios are the two validated
> composition-aware methods; they agree closely, and TMM is the native normalization of our
> edgeR/limma-voom pipeline. We use log2-CPM only for visualization, never for the test
> itself."*

---

## References

[4] Robinson MD, Oshlack A. *A scaling normalization method for differential expression
analysis of RNA-seq data.* Genome Biology. 2010;11:R25. doi:10.1186/gb-2010-11-3-r25. *(TMM.)*

[5] Wagner GP, Kin K, Lynch VJ. *Measurement of mRNA abundance using RNA-seq data: RPKM
measure is inconsistent among samples.* Theory in Biosciences. 2012;131(4):281–285.
doi:10.1007/s12064-012-0162-3. *(RPKM inconsistent across samples; advocates TPM.)*

[6] Bullard JH, Purdom E, Dudoit S, Speed TP. *Evaluation of statistical methods for
normalization and differential expression in mRNA-Seq experiments.* BMC Bioinformatics.
2010;11:94. doi:10.1186/1471-2105-11-94. *(Total-count/RPKM scaling is sensitive to a few
high-count genes; robust between-sample normalization does better.)*

[7] Anders S, Huber W. *Differential expression analysis for sequence count data.* Genome
Biology. 2010;11:R106. doi:10.1186/gb-2010-11-10-r106. *(DESeq median-of-ratios size
factors.)*

[8] Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for
RNA-seq data with DESeq2.* Genome Biology. 2014;15:550. doi:10.1186/s13059-014-0550-8.
*(DESeq2.)*

[9] Li B, Dewey CN. *RSEM: accurate transcript quantification from RNA-Seq data with or
without a reference genome.* BMC Bioinformatics. 2011;12:323. doi:10.1186/1471-2105-12-323.
*(Origin of the TPM unit.)*

[10] Dillies M-A, Rau A, Aubert J, et al. *A comprehensive evaluation of normalization
methods for Illumina high-throughput RNA sequencing data analysis.* Briefings in
Bioinformatics. 2013;14(6):671–683. doi:10.1093/bib/bbs046. *(Benchmark recommending
TMM/DESeq-type methods over total-count/RPKM for DE.)*
