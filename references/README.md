# references/

Reference library for the psoriasis project: annotated bibliography + the
source PDFs and machine-readable full text.

```
references/
├── README.md          ← you are here (file index)
├── references.md      ← the annotated library: per-paper claim + how it is used
├── pdfs/              ← source PDFs (citation-style filenames)
└── fulltext/          ← extracted text for grep/search
```

## PDFs on file

| File | Citation | Role in the project |
|---|---|---|
| `Heidenreich_2009_IJEP_angiogenesis-drives-psoriasis.pdf` | Heidenreich R, Röcken M, Ghoreschi K. *Int J Exp Pathol* 2009;90(3):232–248. DOI 10.1111/j.1365-2613.2009.00669.x | **Supports the vascular thesis.** Vessel formation starts with early psoriatic change and resolves before the epidermis; psoriatic capillaries convert to a venous phenotype — matches our top gradient genes ACKR1/PLVAP (post-capillary venule markers). |
| `Garshick_2020_ATVB_platelets-induce-endothelial-inflammation-COX1.pdf` | Garshick MS, Tawil M, Barrett TJ, et al. *Arterioscler Thromb Vasc Biol* 2020;40(5):1340–1351. DOI 10.1161/ATVBAHA.119.314008. PMID 32131611 | **Supplies the causal/interventional arm.** Platelets from psoriasis patients adhere to endothelium and induce IL8/IL1β/COX-2; randomised low-dose aspirin cut endothelial inflammatory transcripts >70% (NCT03228017). |
| `Cooper_1990_JID_IL1-activity-psoriatic-skin.pdf` | Cooper KD, et al. *J Invest Dermatol* 1990. | **Resolves the IL-1 direction paradox.** Historical precedent that IL-1 bioactivity is *reduced* in lesional psoriatic skin with a shift toward inhibitory activity — anticipates our endothelial IL1R1↓ / IL1RN↑ finding. |

## Full text on file (searchable)

| File | Note |
|---|---|
| `Garshick_2020_ATVB_fulltext.txt` | fetched via PMC |
| `Cai_2019_JID_IL1b-IL1R-pathway_fulltext.txt` | Cai Y, et al. *J Invest Dermatol* — DOI 10.1016/j.jid.2018.07.025, PMC6392027, PMID 30120937. Read carefully: it measures IL-1R **pathway activity** and γδ-T-cell/keratinocyte requirement, **not** IL1R1 transcript in endothelium — so it does not conflict with our result. |

## Cited but not held as PDF

Full citations and the specific claim taken from each are in `references.md`:
Sun et al. 2022 (Scissor), Wilks et al. 2021 (recount3), DerSimonian & Laird
1986, Ma et al. 2023 (GSE173706), Tsoi et al. 2019 (SRP165679), the psoriatic-
march reviews (Front Med 2022;9:864185; PMC9744099), IL-36/early-IL-23
(PMC7190273), S1PR3–STAT3 (s41419-025-07358-w), and the mLOY papers
(Nat Rev Cardiol 2023; Science 2022).

## Provenance note

Cooper 1990 is a scanned PDF with no text layer; its content was read visually,
so verify exact values against the page before quoting numbers. The Heidenreich
entry cites Braverman & Sibley 1982 second-hand through the review.
