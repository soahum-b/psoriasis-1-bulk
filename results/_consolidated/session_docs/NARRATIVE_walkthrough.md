# Psoriasis project — narrative walkthrough

*A guided tour of everything generated so far, in the order the story actually unfolds.
Written to be read start-to-finish, or used as a map to jump into any one arm. Every number
here is traceable to a saved artifact; where a claim is soft or a caveat matters, it is flagged
inline rather than buried. Clinical interpretation is preliminary and not for patient care.*

---

## 0. The one-paragraph story

We set out to find druggable targets in psoriasis by moving from **genes → pathways → proteins**,
with STAT3 as the lead candidate. What the data actually delivered is more interesting than a
single-gene story: across three independent analytical arms, the dominant, reproducible signal in
psoriatic progression is **vascular** — an IL-1-responsive endothelial compartment — while **STAT3
turns out to be a passenger, not a driver** at single-cell resolution. IL-1β itself, the other lead,
is the *least* reliable transcript in the whole panel; its pathway acts through the inflammasome and
through IL-36, not through IL1B mRNA. The through-line that emerged is a **paracrine circuit**:
dendritic-cell-made IL-1β acting on IL1R1⁺ endothelium that runs STAT3/inflammasome programs — the
same "primed vasculature" the cardiovascular literature ties to the psoriasis→atherosclerosis link.
That reframing — from "STAT3 driver" to "IL-1-primed vascular program" — is the intellectual arc of
the project so far.

---

## 1. How the project is organized (so the pieces make sense)

There are **three scientific arms**, each with its own writeup, plus a set of living project docs.
They are deliberately kept separate because they answer different questions on different data:

| Arm | Question | Data | Main document |
|---|---|---|---|
| **A — Bulk meta-analysis** | Which genes are up/down in lesion vs normal, reproducibly across studies? | 4 bulk RNA-seq studies (recount3) | `IL1B_STAT3_expanded_analysis.md` |
| **B — Scissor gradient** | Which individual *cells* track the normal→peri-lesional→lesional progression? | 89k-cell single-cell atlas (GSE173706) | `WHITEPAPER.md` |
| **C — Theory 1: vascular IL-1** | Is there an IL-1-responsive vascular population linking psoriasis to atherosclerosis? | Same 89k single-cell object | `theory1_endothelial_IL1_and_blood_arm.md` |

Two side-investigations sit on top of these: a **loss-of-Y** check (`notes/lossofY_skin_check.md`)
and a **blood arm** that was attempted and dropped. The project record lives in `docs/DECISIONS_LOG.md`,
`docs/METHODS.md`, `notes/PLAN.md`, `notes/lab_notebook.md`, and `CHECKPOINT.md`, all under
`github.com/soahum-b/psoriasis-1-bulk`.

---

## 2. Arm A — The bulk meta-analysis (where STAT3 still looks strong)

**The setup.** This is the classic, well-powered question: pool multiple independent
psoriasis studies and ask which genes reliably differ between lesional (PP) and normal (NN) skin. We
deliberately **drop the peri-lesional tier** here — it is underpowered for gene-level effects (only
2 studies contribute it, and only *IL18* reached significance across a 49-gene panel) — and ask the
clean two-group question. Four treatment-free studies qualified (SRP035988, SRP165679, SRP126422,
SRP065812; **141 normal + 145 lesional**). Each study went through an identical limma-voom pipeline,
combined by DerSimonian-Laird random-effects meta-analysis. "k=4" means four independent studies
voting on each gene.

**What it found — three clean messages:**
1. **STAT3/JAK-STAT is the tightest, most reproducible signal.** STAT3 +1.06, STAT1 +1.70,
   JAK3 +1.12, TYK2 +0.23, all significant across 4 studies. This is the arm where the STAT3
   hypothesis looks strongest — and it is already druggable via JAK inhibition clinically.
2. **The IL1B ligand gene is the *least* reliable member of the panel** (+1.09 but non-significant,
   I²=94%; the newest cohort even pulls it *down*). Yet the inflammasome *machinery* around it is
   solid — CASP1, PYCARD, IL1RN, AIM2, CASP5 all significant. Message: IL-1β pathway activity is
   carried by the assembled inflammasome and its output, **not** by a clean IL1B mRNA rise.
3. **IL-36 dwarfs IL-1β** (IL36A +9.0) — within the IL-1 family, IL-36 is the dominant psoriasis
   axis. The strongest signal of all is the shared downstream output (S100A7/8/9, DEFB4A, LCN2, CCL20).

**The honest craft here:** we started at k=3 and added a 4th study mid-analysis. That shrank large
fold-changes toward zero (S100A9 7.9→6.1) and flipped three genes (SOCS3, IL18, CXCL2) to
non-significant. This is **not** signal loss — it is the 4th study exposing genuine between-study
disagreement (all three flips now have I²=95–98%). Every gene's *direction* was preserved. Adding an
independent lab regresses lab-specific magnitudes toward an honest pooled value.

> **The key takeaway for the whole project:** in bulk, STAT3 looks like a real target. Hold that
> thought — Arm B complicates it.

---

## 3. Arm B — The Scissor gradient (where STAT3 becomes a passenger)

**The idea, and why it is clever.** Instead of the usual two-state (lesion vs normal) contrast, we
treat biopsy site as an **ordinal gradient**: normal (NN) < peri-lesional (PN) < lesional (PP). The
peri-lesional margin is the advancing edge of a plaque — where health is *becoming* disease — and a
two-state design is blind to it. We then use **Scissor** to ask, of 89,058 individual cells in a
single-cell atlas, *which ones co-vary with that clinical gradient*. The crucial design safeguard:
the phenotype is a **clinical biopsy-site label**, not a molecular score computed from the same
cells — so a flagged cell cannot merely be restating how it was labelled (no circularity).

**What it found:**
- Selected cells are **monotonic on the gradient** (mean tier: Scissor− 0.79 < background 1.38 <
  Scissor+ 1.43), and the gradient-tracking fraction **peaks at the peri-lesional tier** — exactly
  the signature of an intermediate progression state.
- Two independent significance controls (a reliability test and a permutation null) both give
  p ≈ 0 — the signal is real phenotype structure, not graph artifact.
- **Endothelial cells dominate the gradient-tracking population** — 5.2× enriched, the top lineage,
  and an orthogonal bulk deconvolution independently confirms endothelial proportion rises NN→PN→PP.
  The program is **vascular-led** (top genes CCL14, ACKR1, RAMP3, PLVAP — all endothelial).
- **STAT3 does not survive.** On a 20k-cell "backbone" subset STAT3 looked significant (log2FC 0.43,
  padj 0.018), but at the full 89k-cell census it is **not significant** under proper BH-FDR
  correction (log2FC 0.30, q=0.13; 42.7% vs 44.9% of cells expressing). The backbone signal was a
  small-subset effect that dissolved at full resolution.

> **The pivot of the whole project lives here:** the bulk arm says "STAT3 up"; the single-cell arm
> says "yes, but it's expressed nearly everywhere and doesn't *track progression* — the cells that
> track progression are the blood vessels." STAT3 is embedded in a broader vascular/inflammatory
> circuit, not the driver of it.

**Two methodological honesties worth understanding:**
- **Backbone vs full census.** Early results were on a 20k-cell subset (the "backbone") for speed;
  the full 89k run later confirmed the *structure* (endothelial-led, monotonic) but *revised* STAT3
  to non-significant. The full census also re-tuned a smoothing parameter (alpha 0.40→0.20). The
  paper's figures still show the backbone pending a refresh — a flagged limitation.
- **The solver.** Scissor's official engine (APML1) would not build in our environment, so we
  reimplemented the identical math in pure R on glmnet, validated on synthetic data (recovered
  100/100 true cells). It is a validated equivalent, not the canonical binary — a cross-check is staged.

---

## 4. Arm C — Theory 1: the IL-1-primed vasculature (where it all connects)

This is the clinician's own hypothesis, and it is where the bulk and single-cell arms fuse into one
mechanism. **Original theory:** IL-1β on the vessels/circulation drives IL-1 signaling → STAT3 +
inflammasome → and this vascular inflammation is what predisposes psoriasis patients to
**atherosclerosis** (the "psoriatic march"). The data required two corrections but supported the core:

1. **IL-1β is not made by the endothelium — it is dendritic-cell-restricted** (48% of DCs express
   it, vs 0.2% of endothelial cells). So the vessels don't *produce* IL-1β. This upgrades the theory
   to a **paracrine model**: DC/systemic IL-1β acting *on* the vessels — a stronger, more specific claim.
2. **The endothelium is the *responder*.** IL1R1 (the receptor) is on ~29% of endothelial cells. And
   critically: **IL1R1⁺ endothelial cells are 2.3× enriched in the gradient-tracking Scissor+ set**
   (OR 2.32, p=1.4×10⁻³⁸), and they co-express far more STAT3 (66% vs 40%), IL6, NF-κB, CASP1, GSDMD,
   PYCARD than IL1R1⁻ endothelial cells. **This is the causal chain made visible** — the vessel cells
   bearing the IL-1 receptor are the same ones running STAT3 and the inflammasome.
3. **Direction nuance:** the receptor is *early-high* — IL1R1 is actually *down* in fully lesional
   endothelium (51%→29%→27% across NN→PN→PP), consistent with IL-1 as an **initiating** signal that
   desensitizes as lesions establish, not a maintaining one.

This is well-grounded in literature (psoriatic march; ~50% increased CV risk; IL-1↔STAT3 axis), and
the specific cell-level wiring (IL1R1⁺ endothelium = the STAT3/inflammasome-active, progression-
tracking population) is a less-explored, publishable angle. Figure: `fig_theory1_endothelial_IL1.png`.

> **How this rescues STAT3 and IL-1β without overclaiming:** neither is the headline driver, but both
> sit inside the vascular program. STAT3 is what the IL-1-primed endothelium runs; IL-1β is the
> (paracrine) trigger. The target is not "STAT3" or "IL-1β" in isolation — it is the IL-1-responsive
> vascular circuit.

---

## 5. The two side-investigations (and why both are negatives worth keeping)

**The blood arm — attempted, then dropped.** The theory points at circulation, so we tried to add a
blood dataset. We ingested ERP110814 (10 treatment-naive plaque-psoriasis blood samples, recovered
via the ArrayExpress timepoint annotation). But it has no healthy controls, and the only in-recount3
healthy whole-blood reference (GTEx) is **perfectly study-batch confounded** with it — disease status
is 100% aliased with lab/protocol (globin fraction 0.05% vs 37%, ~10-year age gap). No matching fixes
a confounder collinear with the exposure, so the contrast is uninterpretable and was **dropped**. A
real blood arm would need a *single* study containing both plaque-psoriasis and healthy blood,
processed together — not available in recount3. (Separately: the known recount3 psoriasis-blood
studies are all Generalized Pustular Psoriasis, a different subtype, and excluded — though this was a
screen of 3 pre-selected accessions, not an exhaustive census.)

**Loss-of-Y (theory 2 down payment) — clean negative in skin.** We can measure Y-chromosome gene
*expression* from RNA-seq (we did), but not genomic mLOY (that needs DNA/blood). In 21 male donors /
65,524 cells, the fraction of "Y-silent" cells is **driven by sequencing dropout** (31.5% in shallow
cells → 1.4% in deep cells) and **does not track disease tier**. So there is no loss-of-Y-expression
signal above the technical floor in skin — expected, since mLOY is a *blood*-cell phenomenon. A real
mLOY test needs blood.

> **Why keep negatives?** Both close off tempting-but-invalid directions with evidence, so we don't
> chase them later. The blood arm looks appealing until you see the batch confound; LOY looks testable
> until you see it's all dropout. Recording *why* they fail is as valuable as a positive result.

---

## 6. The intellectual arc, in one view

```
STAT3 as lead target (hypothesis)
        │
   Arm A (bulk):  STAT3/JAK-STAT robustly up  ── looks like a target ✓
        │         IL1B unreliable; IL-36 dominant; inflammasome up
        ▼
   Arm B (single-cell): STAT3 does NOT track progression at full census ✗
        │              the cells that track it are ENDOTHELIAL (vascular-led)
        ▼
   Arm C (theory 1): the vascular program is an IL-1-RESPONSIVE endothelium
                      IL1R1⁺ vessels run STAT3 + inflammasome (paracrine, DC-sourced IL-1β)
                      → ties psoriasis to atherosclerosis (psoriatic march)
        │
        ▼
   Reframed target: not STAT3 alone, not IL-1β alone —
                    the IL-1-primed vascular circuit (STAT3/inflammasome downstream)
```

The story is not "we confirmed STAT3." It is "we tested STAT3 rigorously, found it is a passenger in
a vascular-led program, and identified the vascular IL-1 circuit as the more defensible target — with
a direct line to the cardiovascular comorbidity that matters clinically."

---

## 7. What stands, what's soft, what's open

**Stands (well-supported):**
- Vascular/endothelial dominance of the progression program (Arm B, two significance controls +
  orthogonal deconvolution).
- STAT3/JAK-STAT up in bulk (Arm A, k=4).
- IL1R1⁺ endothelium = the STAT3/inflammasome-active, gradient-tracking population (Arm C, OR 2.32).
- IL1B mRNA is unreliable; IL-36 + inflammasome are the real IL-1-family signal (Arm A).

**Soft / caveated:**
- STAT3 single-cell result depends on backbone-vs-full-census scale (resolved: n.s. at full census).
- Scissor solver is a validated re-implementation, not the canonical binary.
- Single bulk anchor for the gradient (SRP165679); NNLS deconvolution not benchmarked.
- High between-study heterogeneity on several bulk genes.
- WHITEPAPER §4 figures still show the 20k backbone, not the 89k object.

**Open (not yet done):**
- Theory-2 Part A: **sex-stratified analysis of the skin data** (feasible now — the clean next step).
- WHITEPAPER HTML/PDF re-render + §4 figure refresh from the 89k object.
- Compiled-Scissor (APML1) cross-check; benchmarked deconvolution (Stage 1 of the two-stage design);
  sequence-level Sei-LLRA arm.
- A proper blood arm (needs a non-recount3 internally-controlled dataset).

---

## 8. Questions this walkthrough should let you ask

Use these as entry points — each maps to a section and a document:
- *"Why did STAT3 look significant then not?"* → §3 (backbone vs full census, BH-FDR).
- *"If IL1B is unreliable, why do we still care about IL-1?"* → §2 + §4 (inflammasome/IL-36 carry it;
  the receptor on vessels is the story).
- *"What exactly is the druggable target now?"* → §4 + §6 (the IL-1-primed vascular circuit).
- *"How does this connect to heart disease?"* → §4 (psoriatic march, IL1R1⁺ endothelium).
- *"What's the difference between the two-state and gradient designs?"* → §3.
- *"Why did the blood arm fail?"* → §5 (study-batch confounding).
- *"Is the peri-lesional margin actually special?"* → §3 (Scissor+ peaks at PN).
- *"What should we do next?"* → §7 (sex-stratified skin pass).
