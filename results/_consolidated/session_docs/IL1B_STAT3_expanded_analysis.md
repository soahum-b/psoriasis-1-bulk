# IL-1β/inflammasome and STAT3 in psoriatic lesions — a focused two-group meta-analysis

*Companion section to the psoriasis bulk meta-analysis. Designed to slot into `psoriasis_meta_analysis_whitepaper.md`. All quantitative values computed from the saved analysis artifacts; clinical interpretation is preliminary and not intended to guide patient care.*

---

## Motivation

An earlier arm of this project treated biopsy site as an **ordinal gradient** (normal `NN` < peri-lesional `PN` < lesional `PP`) to expose the advancing margin. That framing is valuable for the progression question, but it is underpowered for a **focused mechanistic contrast**: in the ordinal set, only two studies contribute peri-lesional samples, and across a 49-gene IL-1β/STAT3 panel only a single gene (*IL18*) reached significance in the `PN vs NN` contrast. The peri-lesional arm cannot resolve gene-level effects for this panel.

Here we deliberately **drop the peri-lesional tier** and ask the well-powered two-group question — **psoriatic lesion (PP) vs normal skin (NN)** — focusing on two interacting modules: the **IL-1β/inflammasome** axis and the **STAT3/JAK-STAT** axis, with IL-36, Th17, and the shared downstream output (S100/defensin/chemokine) genes as context. These two modules are not on a single linear pathway; they are parallel arms of one **TNF → IL-1β/IL-6 → Th17 (STAT3) → IL-17A** feed-forward circuit, so reading them side by side is the point.

## Design

Dropping the peri-lesional tier relaxes the study-eligibility constraint from "must contain NN **and** PP" (the actual constraint on the ordinal set — one qualifying study, SRP035988, already contributes NN+PP with no PN) to simply "must contain lesional and normal skin." We screened all 23 human psoriasis studies in the recount3 SRA pool against sample-level metadata (`study_eligibility_PPvsNN.csv`). Most large candidates disqualify on defensible grounds — no normal samples, in-vitro/cell-line/explant material, or blood compartment. One clean addition qualified:

- **SRP065812** (16 healthy-donor normal + 18 lesional, median 43 M reads). It is an adalimumab study; we retained **only the 18 pre-treatment lesional samples** and dropped the 18 post-treatment ones.

We excluded **ERP110816** (86 uninvolved + 92 lesional, deep) despite its size: every patient is on **etanercept (anti-TNF)**, which pharmacologically suppresses the exact TNF→IL-1β→STAT3 axis under study, and the recount3-deposited metadata carries **no visit/timepoint label**, so a treatment-naive baseline cannot be reliably isolated (10 patients, 178 longitudinal biopsies, PASI spanning 0.9–17+ with no marker of which biopsy is pre-dose).

The expanded set therefore contains **4 independent studies, all treatment-free** (SRP035988, SRP165679, SRP126422, SRP065812; 141 normal + 145 lesional). Per-study differential expression was fit with the **identical pipeline** used for the existing studies (`filterByExpr` → TMM → limma-voom, effect = log₂FC, SE = log₂FC/t), and the four studies were combined by the same closed-form **DerSimonian-Laird random-effects meta-analysis**. This keeps every effect size comparable across the k=3 → k=4 upgrade.

## Result

![IL-1β/inflammasome and STAT3 modules, lesional vs normal, 4-study meta-analysis]({{artifact:art_6ceb085d-1b9b-4d53-a13d-4fe59cfce924}})

*Expanded 4-study meta-analysis, psoriatic lesion (PP) vs normal (NN), peri-lesional excluded. Points are random-effects pooled log₂ fold-changes with 95% CIs; filled = FDR<0.05, faded = n.s.; ‡ marks high between-study heterogeneity (I²≥75%); ᵏ⁼ⁿ marks genes supported by fewer than 4 studies.*

**The STAT3/JAK-STAT arm is uniformly, robustly activated** and is confirmed by the 4th independent cohort: *STAT3* +1.06, *STAT1* +1.70, *JAK3* +1.12, *TYK2* +0.23 (all FDR<0.05, k=4). *JAK1* is the only member down (−0.23), consistent with a constitutive rather than induced role. This is the tightest, most reproducible signal in the panel.

**The inflammasome machinery is up, but the IL1B ligand gene itself is not.** The processing apparatus is significant — *CASP1* +0.59, *PYCARD* +1.09, *IL1RN* +0.83, *AIM2* +2.84, *CASP5* +3.54, *GSDMD* — but *IL1B* is non-significant (+1.09, FDR n.s., I²=94%), and the new cohort pulls it *down* (SRP065812 log₂FC −0.28). Across four independent studies, IL1B transcript is the single least reliable member of the panel — a well-replicated negative worth stating plainly before any target hypothesis rests on it. *IL1A* and *IL18* are down.

**IL-36 dwarfs IL-1β** (*IL36A* +9.0, *IL36G* +4.7), reinforcing that within the IL-1 family it is IL-36, not IL-1β, that is the dominant psoriasis axis. The **downstream output** (*S100A7/8/9*, *DEFB4A*, *LCN2*, *CCL20*) is the strongest signal overall — the shared readout both arms converge on.

### What adding the 4th study changed

| Gene | Module | k=3 log₂FC | k=4 log₂FC | k=4 sig | I² (k=4) |
|---|---|--:|--:|:--:|--:|
| STAT3 | STAT3/JAK | +1.25 *** | **+1.06** | *** | 92% |
| STAT1 | STAT3/JAK | +1.86 *** | +1.70 | *** | 83% |
| JAK3 | STAT3/JAK | +1.54 *** | +1.12 | * | 95% |
| SOCS3 | STAT3/JAK | +2.16 *** | +0.86 | **n.s. (flip)** | 98% |
| IL1B | inflammasome | +1.63 n.s. | +1.09 | n.s. | 94% |
| CASP1 | inflammasome | +0.52 *** | +0.59 | *** | 72% |
| AIM2 | inflammasome | +3.35 *** | +2.84 | *** | 82% |
| IL18 | inflammasome | −0.37 *** | −0.05 | **n.s. (flip)** | 95% |
| IL36A | IL-36 | +10.06 *** | +9.03 | *** | 73% |
| S100A9 | output | +7.94 *** | +6.08 | *** | 98% |
| CXCL2 | output | +2.29 * | +1.19 | **n.s. (flip)** | 98% |

Two systematic effects, both expected and healthy:

1. **Large fold-changes shrank toward zero** (S100A9 7.9→6.1, SOCS3 2.2→0.86). Adding an independent lab/protocol regresses lab-specific magnitudes toward a more honest pooled value. The **direction of every gene is preserved**; no headline conclusion reverses.
2. **Three genes flipped to non-significant** — *SOCS3*, *IL18*, *CXCL2* — each now with I² of 95–98%. This is not a loss of signal but an **exposure of disagreement**: with only 3 studies the between-study conflict was hidden; the 4th made it visible, and the random-effects model correctly widens the CI. *SOCS3* is the notable one — a direct STAT3 target whose magnitude is unstable across cohorts, itself a finding.

## Interpretation for the IL-1β/STAT3 hypothesis

Both modules are elevated in psoriatic lesions, but they are not equally trustworthy targets. **STAT3/JAK-STAT is the tight, reproducible, druggable arm** (and is already being drugged clinically through JAK inhibition). **IL-1β at the ligand-transcript level is the least reliable signal** in the panel — its pathway activity in psoriasis is carried by the assembled inflammasome and its downstream output, not by a clean *IL1B* mRNA increase, and IL-36 is the dominant IL-1-family axis. A target hypothesis built on IL-1β should therefore rest on the inflammasome/output module or on IL-36, not on *IL1B* expression itself.

## Limitations

- **Healthy-anchored expansion.** SRP065812 adds healthy controls, so the pooled "normal" remains dominated by healthy donors. The healthy-vs-uninvolved control-type moderator we designed rests almost entirely on SRP165679 (the one study carrying both control types); it is estimable but under-powered and should be reported as exploratory.
- **Modest k for several genes.** *AIM2*, *CASP5*, *IL36A* rest on k=3 and *IL17A*, *RORC* on k=2 (these genes fail `filterByExpr` in the shallower studies); their estimates are less reliable than the k=4 majority and are flagged on the figure.
- **High heterogeneity is pervasive** (many genes I²>75%), reflecting genuine cross-cohort differences in lesion severity, platform, and depth — appropriately handled by the random-effects model but a caveat on any single pooled magnitude.

## Artifacts

- `fig_il1b_stat3_NNvsPP_4study.png` — the both-modules forest plot at k=4 (this section's figure).
- `meta_de_PPvsNN_4study.csv` — full genome-wide 4-study meta-analysis table.
- `meta_PPvsNN_k3_vs_k4_panel.csv` — side-by-side k=3 vs k=4 for every panel gene (Δlog₂FC, I², significance flips).
- `study_eligibility_PPvsNN.csv` — the 23-study recount3 screen with per-study decisions.
- `fig_il1b_stat3_NNvsPP.png` — the original k=3 version (peri-lesional excluded), for comparison.
