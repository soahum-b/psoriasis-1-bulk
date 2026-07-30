### 11.5a Robustness of the STAT3 gene-level estimate

The STAT3 forest pools three studies with considerable heterogeneity (I-squared = 91%). That is
exactly the regime - few studies, large I-squared - where the DerSimonian-Laird (DL) confidence
interval is known to be too narrow: DL builds its interval from a normal quantile while treating
tau-squared as if it were known exactly, but with only three studies tau-squared is itself
estimated with wide uncertainty, so the nominal 95% interval under-covers [44]. The
recommended remedy is the **Hartung-Knapp-Sidik-Jonkman (HKSJ)** modification, which replaces the
normal quantile with a t-distribution on k-1 degrees of freedom and re-weights, propagating the
tau-squared uncertainty into the interval [44, 45]. We therefore report HKSJ as the primary method
for this forest and retain DL as a sensitivity analysis, following the recommendation that analysts
using the modification be ready to justify it and show the conventional result alongside [45].

The point estimate is stable across methods; the interval is not (Figure 31A):

| Method | Estimate | 95% CI | p |
|---|---|---|---|
| DL + z (original) | +1.25 | [+0.81, +1.69] | 3e-08 |
| DL + Knapp-Hartung | +1.25 | [+0.02, +2.48] | 0.049 |
| **REML + HKSJ (primary)** | **+1.21** | **[-0.11, +2.53]** | **0.059** |

Under HKSJ the STAT3 gene-level interval widens roughly threefold and its lower bound crosses
zero (p = 0.059). The original p = 3e-08 was an artifact of the z-based interval at k = 3, not a
robust magnitude claim. The width is driven almost entirely by the smallest cohort (SRP126422,
n = 8, SE approximately 0.48), which carries negligible weight but, through its influence on
tau-squared, lets the t-correction open the interval; the two well-powered cohorts each place
STAT3 clearly up with tight intervals.

Crucially, this tempers the *magnitude* claim without touching the *direction* claim, which is
what the meta-analysis was built to test. Leave-one-out re-pooling (Figure 31B) keeps STAT3
positive in every case - +0.95, +1.39, and, with the large anchor cohort (SRP035988) removed,
**+1.12** - so the up-regulation is genuinely cross-cohort and not an artifact of the single
dominant study. The defensible statement is therefore stronger than a lone p-value: **STAT3's
up-regulation is reproducible across all cohorts (robust to leave-one-out, including
anchor removal); its pooled magnitude is approximately +1.2 but imprecisely estimated at k = 3,
which few-study-appropriate inference makes explicit.**

![Figure 31. STAT3 forest robustness]({{artifact:art_c793f88c-9806-46f9-95d2-0b570a3a8e50}})

*Figure 31. (A) Pooled STAT3 effect under three methods: the estimate is stable near +1.2 but the
interval triples in width and crosses zero once few-study uncertainty is honoured (HKSJ). (B)
Leave-one-out re-pooling (DerSimonian-Laird): every re-pool stays positive (+0.95 to +1.39), and
dropping the large anchor cohort keeps STAT3 up at +1.12 - the direction is not anchor-driven.*


---
**References for this section**

[44] IntHout J, Ioannidis JPA, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. *BMC Med Res Methodol.* 2014;14:25. doi:10.1186/1471-2288-14-25.

[45] Jackson D, Law M, Rücker G, Schwarzer G. The Hartung-Knapp modification for random-effects meta-analysis: A useful refinement but are there any residual concerns? *Stat Med.* 2017;36(25):3923-3934. doi:10.1002/sim.7411.
