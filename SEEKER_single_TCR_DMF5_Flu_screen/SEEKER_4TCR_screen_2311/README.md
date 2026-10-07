# SEEKER 4-TCR screen (2311 library)

Four single-TCR SEEKER screens run in parallel against the same full 9-mer NNK library
(SCD format, three enrichment rounds R1–R3): the **DMF5** melanoma TCR, the **1G4**
NY-ESO-1 TCR, an **EBV** (BMLF1 GLCTLVAML-reactive) TCR, and a **Flu** (M1
GILGFVFTL-reactive) TCR.

## Repository layout

```
SEEKER_4TCR_screen_2311/
├── data/
│   ├── top10000_combine_9mer_2311-4TCR-{DMF5|1G4|EBV|Flu}-SCD-round123.csv
│   │        one top-10,000 table per TCR (sorted by count_round3), columns:
│   │        nt / aa / length, count_round1-3, count_pct_round1-3,
│   │        r2_r1_ratio, r3_r2_ratio
│   └── Selected_Peptides.csv     selected-peptide list (306 peptides; aa, tcr_name)
├── code/
│   ├── make_selected_peptides.py        # builds Selected_Peptides.csv (306 peptides)
│   ├── Sanky_aa_usage_all_TCR.py        # lined heatmaps + frequency matrices
│   │                                    #   (peptide set read from Selected_Peptides.csv)
│   ├── Sanky_aa_usage.py                # archived Flu script (kept as-is)
│   └── {DMF5,1G4,EBV,Flu}_plot_seqlogo.R  # per-TCR seqlogo scripts (ggseqlogo)
└── results/
    ├── heatmaps/    lined_{TCR}_4TCR_heatmap.pdf (4) + frequency_{TCR}_4TCR_heatmap.csv (4)
    └── seqlogos/    bold_top200_9mer_{TCR}_mutation.pdf (4)
```

## Selected-peptide criteria

| TCR | Criterion | Peptides |
|---|---|---|
| DMF5 | top200, r3_r2_ratio ≥ 4, deduplicated | 147 |
| 1G4  | top150, r3_r2_ratio ≥ 4 | 105 |
| EBV  | top50, r3_r2_ratio ≥ 10 | 48 |
| Flu  | top25, r3_r2_ratio ≥ 40 | 6 |
| **total** | (no overlap between TCRs) | **306** |

The lined heatmaps and the frequency matrices are generated directly from
`Selected_Peptides.csv`, so all four TCR figures share the same peptide set.

## Reproduce

```bash
# 1. Selected peptides
python code/make_selected_peptides.py                     # writes data/Selected_Peptides.csv

# 2. Lined heatmaps + frequency matrices (pandas, numpy, matplotlib)
python code/Sanky_aa_usage_all_TCR.py                     # writes results/heatmaps/

# 3. Seqlogos (R: ggseqlogo, tidyverse)
Rscript code/DMF5_plot_seqlogo.R                          # repeat for 1G4 / EBV / Flu
```
