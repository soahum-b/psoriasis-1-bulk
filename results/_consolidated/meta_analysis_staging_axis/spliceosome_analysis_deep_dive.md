# Spliceosome & alternative-splicing analysis: a deep dive

*Companion to the psoriasis meta-analysis white paper (§9, STAT3 isoform splicing). This document
explains, from first principles, what the spliceosome is, what alternative splicing is, how you
measure it from RNA-seq data, and how its dysregulation shows up in disease. It is written for a
reader new to the topic and goes into the mechanistic and methodological "weeds." The worked STAT3
α/β result lives in the white paper; the concepts behind it live here.*

---

## 1. Why splicing exists at all

A human gene is not a continuous coding sequence. It is a mosaic of **exons** (the segments that
end up in the mature message) interrupted by **introns** (intervening segments that are removed).
The primary transcript made by RNA polymerase II — the **pre-mRNA** — contains both. Before it can
be translated, the introns must be excised and the exons stitched together. That cut-and-paste
operation is **splicing**.

The number of introns is large: a typical human gene has ~8 introns, and some have dozens. Introns
are often far longer than exons (kilobases of intron around ~150-nt exons). So most of a pre-mRNA is
discarded, precisely and repeatedly, every time the gene is transcribed. This is expensive, and
evolution keeps it for a reason: splicing is not merely removal, it is a **decision point**. The
same pre-mRNA can be spliced in different ways to yield different mature mRNAs — and therefore
different proteins — from one gene. That is **alternative splicing**, and it is why ~20,000 human
genes can encode well over 100,000 distinct proteins. STAT3α versus STAT3β (white paper §9) is
exactly one such decision: one gene, one pre-mRNA, two proteins with opposite functions, selected
by which splice site is used.

---

## 2. The chemistry: what defines a splice, and the two transesterification steps

Splicing is defined by short sequence signals at the intron boundaries, recognised with high
precision:

- The **5′ splice site (donor)** — almost always a **GU** at the intron's 5′ end (the "GT–AG rule"
  at the DNA level).
- The **3′ splice site (acceptor)** — an **AG** at the intron's 3′ end.
- The **branch point** — an adenosine ~18–40 nt upstream of the acceptor.
- The **polypyrimidine tract** — a C/U-rich stretch between branch point and acceptor.

The reaction itself is two sequential **transesterifications** [1]:
1. The branch-point adenosine's 2′-OH attacks the 5′ splice site, cutting the exon–intron bond and
   forming a looped **lariat** intermediate (the intron now joined to itself in a lasso shape).
2. The freed 3′-OH of the upstream exon attacks the 3′ splice site, joining the two exons and
   releasing the intron lariat.

No net energy is stored in the bonds (it is transesterification, not hydrolysis), but the **fidelity
and control** require an enormous molecular machine — because choosing *which* GU pairs with *which*
AG, out of thousands of possibilities in a long pre-mRNA, is the whole game.

---

## 3. The spliceosome: a multi-megadalton, dynamically assembled machine

The **spliceosome** is the ribonucleoprotein complex that catalyses splicing [1,2]. It is not a
pre-formed enzyme like a polymerase; it **assembles anew on each intron** from five small nuclear
ribonucleoprotein particles (**snRNPs**) — **U1, U2, U4, U5, U6** — plus well over 100 associated
proteins. Each snRNP is a small nuclear RNA (snRNA) wrapped in proteins.

The **assembly cycle** (major/U2-dependent spliceosome) proceeds through defined states:
- **E complex** — U1 snRNP recognises the 5′ splice site by base-pairing.
- **A complex** — U2 snRNP binds the branch point.
- **B complex** — the pre-assembled U4/U6.U5 **tri-snRNP** joins.
- **B\*** / **activation** — a major remodelling: U1 and U4 leave, U6 replaces U1 at the 5′ site,
  and the catalytic core forms. This step is driven by RNA helicases (e.g. BRR2/SNRNP200) and the
  **U5 snRNP GTPase EFTUD2** (also called SNU114) — the protein at the centre of our psoriasis lead.
- **C complex** — the two transesterifications occur.
- **Disassembly** — the lariat intron and snRNPs are released and recycled.

The machine is therefore **dynamic**: it is built, activated, catalyses, and disassembles for every
single intron, and each of those transitions is a potential control point where regulation — or
dysregulation — can act. **EFTUD2** is a core, essential component of the **U5 snRNP** and the
U4/U6.U5 tri-snRNP; its GTPase activity is required for the activation remodelling that forms the
catalytic centre. Because it is core (not a peripheral regulator), changes in EFTUD2 are expected to
have broad, pleiotropic effects on splicing efficiency rather than gene-specific effects — a point
that matters for interpreting the white-paper §9 result.

**A note on the minor spliceosome.** A small fraction of introns (~0.3%) use a distinct **U12-
dependent (minor) spliceosome** with U11, U12, U4atac, U6atac (U5 is shared). It is mechanistically
analogous but recognises different consensus sequences; most analyses, and ours, concern the major
spliceosome.

---

## 4. The types of alternative splicing (the vocabulary you must know)

Alternative splicing is categorised into a small number of canonical **event types** — this
vocabulary is what every quantification tool reports, so it is worth knowing precisely:

1. **Exon skipping (cassette exon, SE).** An exon is either included or excluded. The most common
   type in mammals. Quantified as the inclusion level of the cassette exon.
2. **Alternative 5′ splice site (A5SS).** Two competing donors at the same exon's 3′ end shift the
   exon boundary. *(This is the STAT3 α/β mechanism in §9: a shared acceptor with two competing
   donors.)*
3. **Alternative 3′ splice site (A3SS).** Two competing acceptors at an exon's 5′ end.
4. **Mutually exclusive exons (MXE).** Exactly one of two adjacent exons is included, never both.
5. **Intron retention (IR).** An intron is left in the mature mRNA. Often a signal of regulation or
   of splicing failure/stress; frequently targets the transcript for nonsense-mediated decay.
6. **Alternative first/last exons and alternative polyadenylation** — promoter- or 3′-end-driven
   isoform diversity; sometimes grouped separately because they are transcription-/cleavage-driven
   rather than purely spliceosomal.

The functional consequence depends on where the change lands: an in-frame swap changes a protein
domain (STAT3α→β loses the transactivation domain); a frameshift or a retained intron usually
introduces a **premature termination codon (PTC)**, routing the transcript to **nonsense-mediated
decay (NMD)** — so alternative splicing is also a way to *down-regulate* a gene, not only to
diversify it.

---

## 5. The core measurement: PSI (percent spliced-in)

The universal quantity in splicing analysis is **PSI (Ψ), percent-spliced-in** — the fraction of a
gene's transcripts that include a particular splicing choice [3]. For a cassette exon:

```
        reads supporting inclusion
Ψ  =  ----------------------------------------------
      reads supporting inclusion + reads supporting exclusion
```

Ψ ranges 0 (never included) to 1 (always included). A **ΔΨ** between conditions is the effect size —
"this exon is included 20% more in disease." The elegance of PSI is that it is a **ratio internal to
the gene**, so it is automatically normalised for the gene's overall expression: a gene can double
its total output with no change in Ψ, or hold output flat while Ψ swings. This is exactly why §9
reports PSI_β rather than raw β counts — it isolates the *splice decision* from the *expression
level*, which move independently.

**Reads that inform Ψ come in two flavours:**
- **Junction reads** — reads that span an exon–exon boundary, mapping with a gap to the genome.
  These are the *direct* evidence of a splice: a read crossing the α junction proves an α splice
  happened. recount3 provides these as precomputed **junction counts**, which is what §9 uses.
- **Exon-body reads** — reads falling within an exon; informative for inclusion but less specific
  than junction reads.

Junction-based Ψ is the cleaner measurement when junction coverage is adequate (in §9, ≥22 junction
reads in every sample), because each read is an unambiguous vote for one splice choice.

---

## 6. Methods and tools for splicing analysis from RNA-seq

There are three broad strategies, in increasing order of assumption:

**(a) Event-level, annotation-based (count the junctions).** Define the competing junctions for a
known event and compute Ψ directly from junction counts, then test between groups. This is what §9
does by hand, and it is the most transparent approach when you have a *specific hypothesis* about a
*specific event* (STAT3 α/β). The standard automated tool in this family is **rMATS** [4], which
enumerates the five canonical event types genome-wide, computes Ψ per replicate, and tests ΔΨ with a
statistical model that accounts for both technical (read-count) and biological variability. Other
tools: **SUPPA2** [5] (fast, transcript-quantification-based Ψ, scales to hundreds of samples),
**MISO** [3] (Bayesian Ψ estimation; the paper that popularised the PSI framework), **JunctionSeq**
and **dexseq** (exon-bin level).

**(b) Annotation-free, local splicing variation.** Instead of pre-defining events, model the local
graph of splice junctions and detect *any* differential usage. **LeafCutter** [6] clusters
overlapping introns and tests for differential **intron excision ratios** — it needs no transcript
annotation and captures novel/unannotated splicing, at the cost of reporting intron clusters rather
than tidy named event types. **MAJIQ** is a similar graph-based method (local splicing variations,
"LSVs").

**(c) Isoform-level (quantify whole transcripts).** Estimate the abundance of each full-length
transcript isoform, then compare. Tools: **Salmon**/**kallisto** (fast pseudo-alignment transcript
quantification) feeding **DRIMSeq** or **DEXSeq** for **differential transcript usage (DTU)**, or
**IsoformSwitchAnalyzeR** for a full DTU + functional-consequence workflow. This is the most
complete view but the most assumption-heavy: assigning short reads to full isoforms is an
inherently under-determined problem, and estimates are noisier than a single junction ratio.

**Which to use.** For a **targeted, specific isoform question** (STAT3 α/β), event-level junction
counting (a) is the most direct and defensible — you can point to the exact reads. For a **discovery
scan** ("what is differentially spliced anywhere?"), rMATS (a) or LeafCutter (b) genome-wide. For
**isoform-resolution biology across a transcriptome**, the DTU route (c). As always, match the
method to the question; a genome-wide scan and a single-event test answer different things.

**Practical cautions specific to splicing:**
- **Junction coverage is everything.** Ψ from few reads is unstable — filter on a minimum junction
  depth per sample (§9 used ≥ ~10–20 and confirmed a floor of 22).
- **Read length and mappability** matter: junction reads must span the boundary with enough
  overhang on each side to map uniquely.
- **Bulk tissue confounds splicing with cell-type composition.** If two conditions differ in their
  cellular mix, apparent ΔΨ can reflect *which cells are present*, not a splicing change within a
  cell type. This is a real caveat for the §9 EFTUD2 correlation, and single-cell/long-read data are
  the remedy.
- **Long-read sequencing (PacBio/Oxford Nanopore)** sidesteps the short-read isoform-assignment
  problem by reading whole transcripts end to end — the emerging gold standard for isoform
  discovery, at higher cost and lower depth.

---

## 7. Spliceosome dysregulation in disease

Splicing is not just a housekeeping step; its disruption causes disease through several distinct
routes [2,7]:

- **Cis mutations that create/destroy splice sites.** A point mutation in a donor/acceptor or a
  branch point can abolish a splice or activate a **cryptic** one, producing an aberrant transcript.
  A large fraction of disease-causing mutations act by disrupting splicing rather than by changing a
  codon directly.
- **Trans mutations in spliceosome components.** Recurrent somatic mutations in core splicing
  factors — **SF3B1, U2AF1, SRSF2, ZRSR2** — are hallmark drivers of **myelodysplastic syndromes**
  and other cancers [7]; they alter splice-site selection genome-wide, producing mis-spliced
  transcripts that drive the malignant phenotype.
- **Spliceosome-component dosage / developmental syndromes.** Germline haploinsufficiency of core
  factors causes developmental disorders — notably **EFTUD2 mutations cause mandibulofacial
  dysostosis with microcephaly (MFDM)** [8], establishing that EFTUD2 dosage has strong phenotypic
  consequences via broad splicing effects.
- **Splicing-factor expression changes in inflammation/proliferation.** Beyond mutation, the
  *expression* of splicing regulators shifts with cell state — proliferating and inflamed cells
  reprogram their splicing. This is the regime relevant to psoriasis: EFTUD2 is *up* in lesional
  skin (white paper §9) as part of a broadly activated, proliferative keratinocyte program, and the
  question is whether that shift bends specific splice ratios such as STAT3 α/β.

**The therapeutic angle.** Splicing is druggable: **splice-switching antisense oligonucleotides**
(e.g. nusinersen for SMA, which redirects SMN2 exon-7 inclusion) and small-molecule spliceosome
modulators (SF3B1-targeting compounds in oncology) are approved or in trials. This is why isoform-
resolved analysis is not merely descriptive — if a specific splice choice drives a phenotype, it can
in principle be corrected.

---

## 8. How this maps onto the psoriasis analysis (§9)

To close the loop between concepts and the worked result:

- **Event type:** the STAT3 α/β choice is an **alternative 5′ splice site (A5SS)** event — a shared
  acceptor (chr17:42,317,181) with two competing donors (α 42,316,902; β 42,316,852), 50 nt apart.
- **Measurement:** **junction-based PSI_β = β/(α+β)**, computed directly from recount3 junction
  counts — strategy (a), event-level, annotation-anchored, the most transparent choice for a
  specific isoform hypothesis.
- **Result:** PSI_β rises modestly but significantly in lesional skin (ΔΨ ≈ +1.2 pp, p = 0.017)
  while both isoforms rise in absolute terms — the splice decision and the expression level moving
  as separable quantities, exactly the separation PSI is designed to expose.
- **Spliceosome link:** EFTUD2, a **core U5 snRNP GTPase** required for catalytic activation
  (§3), is strongly up in lesional skin; but within groups it correlates *negatively* with PSI_β,
  arguing against a simple "EFTUD2 makes the β brake" model and pointing instead toward EFTUD2
  supporting activating-α production.
- **The key caveat, from §6:** this is bulk tissue, so composition confounds splicing; and it is
  correlational. The definitive follow-ups are isoform-resolved protein/phospho data and single-cell
  or long-read sequencing to remove the cell-type confound.

---

## References

[1] Wahl MC, Will CL, Lührmann R. **The spliceosome: design principles of a dynamic RNP machine.**
*Cell.* 2009;136(4):701–718. doi:10.1016/j.cell.2009.02.009.

[2] Scotti MM, Swanson MS. **RNA mis-splicing in disease.** *Nat Rev Genet.* 2016;17(1):19–32.
doi:10.1038/nrg.2015.3.

[3] Katz Y, Wang ET, Airoldi EM, Burge CB. **Analysis and design of RNA sequencing experiments for
identifying isoform regulation (MISO).** *Nat Methods.* 2010;7(12):1009–1015. doi:10.1038/nmeth.1528.
*(Introduced the PSI/Ψ framework for RNA-seq.)*

[4] Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y. **rMATS: robust and flexible
detection of differential alternative splicing from replicate RNA-Seq data.** *Proc Natl Acad Sci
USA.* 2014;111(51):E5593–E5601. doi:10.1073/pnas.1419161111.

[5] Trincado JL, Entizne JC, Hysenaj G, Singh B, Skalic M, Elliott DJ, Eyras E. **SUPPA2: fast,
accurate, and uncertainty-aware differential splicing analysis across multiple conditions.** *Genome
Biol.* 2018;19(1):40. doi:10.1186/s13059-018-1417-1.

[6] Li YI, Knowles DA, Humphrey J, Barbeira AN, Dickinson SP, Im HK, Pritchard JK. **Annotation-free
quantification of RNA splicing using LeafCutter.** *Nat Genet.* 2018;50(1):151–158.
doi:10.1038/s41588-017-0004-9.

[7] Yoshida K, Sanada M, Shiraishi Y, et al. **Frequent pathway mutations of splicing machinery in
myelodysplasia.** *Nature.* 2011;478(7367):64–69. doi:10.1038/nature10496.

[8] Lines MA, Huang L, Schwartzentruber J, et al. **Haploinsufficiency of a spliceosomal GTPase
encoded by EFTUD2 causes mandibulofacial dysostosis with microcephaly.** *Am J Hum Genet.*
2012;90(2):369–377. doi:10.1016/j.ajhg.2011.12.023. *(EFTUD2 dosage and disease.)*

[9] Schaefer TS, Sanders LK, Nathans D. **Cooperative transcriptional activity of Jun and Stat3β, a
short form of Stat3.** *Proc Natl Acad Sci USA.* 1995;92(20):9097–9101. doi:10.1073/pnas.92.20.9097.
*(STAT3β first described.)*

[10] Caldenhoven E, van Dijk TB, Solari R, et al. **STAT3β, a splice variant of transcription factor
STAT3, is a dominant negative regulator of transcription.** *J Biol Chem.* 1996;271(22):13221–13227.
doi:10.1074/jbc.271.22.13221. *(α vs β definition.)*

---

*Companion documents: `psoriasis_meta_analysis_whitepaper.md` (§9 applies junction-based PSI to
STAT3 α/β), `network_analysis_deep_dive.md`, `pathway_analysis_deep_dive.md`,
`differential_expression_deep_dive.md`, `normalization_deep_dive.md`.*
