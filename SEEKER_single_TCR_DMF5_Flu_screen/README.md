# SEEKER Single TCR peptide screening

Analysis code and data for SEEKER 9-mer DMF5 and Flu-mut peptide screens:

| Screen | Library | Sorting rounds | Input table |
|--------|---------|----------------|-------------|
| **DMF5** | 9-mer NNK | R1–R3 | `data/top10000_combine_9mer_2305_DMF5_SCD_R123.csv` |
| **Flu mut** | 9-mer mutagenesis | R1–R3 | `data/top10000_combine_9mer_2308-Flu-mut-SCD-round123.csv` |
| **4-TCR parallel** (DMF5 / 1G4 / EBV / Flu) | 9-mer NNK | R1–R3 | [`SEEKER_4TCR_screen_2311/`](SEEKER_4TCR_screen_2311/README.md) |

Each input table contains the top 10,000 DNA barcodes ranked by combined abundance across rounds, with:

- `nt` / `aa` / `length` — nucleotide sequence, encoded 9-mer peptide, and length
- `count_roundN`, `count_pct_roundN` — raw and percentage counts per sorting round
- `r2_r1_ratio`, `r3_r2_ratio` — round-to-round enrichment ratios

## Repository layout

```
.
├── code/
│   ├── flu_mutation_screen/   # Flu screen: heatmap, cumulative curves, interactive network
│   ├── dmf5_screen/           # DMF5 screen: heatmap, interactive network
│   └── shared/                # Network pipeline & cross-screen utilities
├── data/                      # Top-10,000 combined count tables (one per screen)
├── results/                   # All outputs: tables, network figures
│   ├── flu_mutation_screen/   #   heatmaps, cumulative curves, selected peptides
│   ├── flu_network/           #   nodes/links/layout CSVs, network figures
│   ├── dmf5_screen/           #   heatmaps, seqlogo, frequency matrix
│   ├── dmf5_network/          #   nodes/links/layout CSVs, network figures
│   └── hamming_distance_distribution_Flu_vs_DMF5.*  # pairwise Hamming distance comparison
└── figures/                   # Publication figures (Fig 2D–2G, Fig 3C–3D)
```

The sub-folder [`SEEKER_4TCR_screen_2311/`](SEEKER_4TCR_screen_2311/README.md) is a
self-contained four-TCR parallel screen package (DMF5 / 1G4 / EBV / Flu, 2311 library) with
its own data, code, results and README.

## Publication figures

| Figure | Content | Source script |
|--------|---------|---------------|
| Fig 2D | Flu cumulative peptide frequency curves | `peptide_per_round_frequency.py` |
| Fig 2F | Flu amino-acid usage heatmap (+ seqlogo) | `sanky_aa_usage.py` / `plot_seqlogo.py` |
| Fig 2G | Flu peptide network, all edges (Hamming ≤ 3), 125 peptides | `plot_network_vector.py` |
| Fig 3C | DMF5 amino-acid usage heatmap (+ seqlogo) | `dmf5_plots.py` / `plot_seqlogo.py` |
| Fig 3D | DMF5 peptide network, Hamming ≤ 3 edges, 295 peptides | `plot_network_dmf5_light.py` |

Network figures mark the top enriched peptide (Flu: WT `GILGFVFTL`, gold node; DMF5: `TLLGFGSCV`, red star).

## Analysis overview

```
top-10k combined tables (data/)
        │
        ├─ per-screen figures: heatmaps, cumulative curves (Flu)
        │
        └─ peptide similarity networks (Hamming distance)
             1. enrichment filtering + peptide-level deduplication
             2. interactive force-directed network (pyecharts)  → *_interactive_network.html
             3. node/link extraction + hierarchical clustering  → nodes.csv / links.csv
             4. publication vector rendering (matplotlib)       → Fig 2G / Fig 3D
```

Alternative edge sets (`dist_le2`, `dist_eq1`, DMF5 `all_edges`) can be rendered from the same tables by editing `VERSIONS` in the plotting scripts.

## Reproduction

All scripts resolve paths relative to the repository root — no configuration needed. Run from the repo root.

```bash
# 1. Per-screen descriptive figures
python code/flu_mutation_screen/sanky_aa_usage.py                # Fig 2F heatmap
python code/flu_mutation_screen/peptide_per_round_frequency.py  # Fig 2D cumulative curves
python code/dmf5_screen/dmf5_plots.py                           # Fig 3C heatmap

# 2. Interactive networks (pyecharts)
python code/flu_mutation_screen/pyechart_flu.py                 # Flu network HTML (125 peptides)
python code/dmf5_screen/pyechart_dmf5.py                        # DMF5 network HTML (295 peptides)

# 3. Network tables (nodes.csv / links.csv / layout)
python code/shared/extract_network_data.py                      # Flu: HTML → nodes.csv / links.csv
python code/shared/patch_root_nodes.py                          # Flu: anchor flag + singleton colors
python code/shared/pipeline_all_networks.py                     # DMF5: clustering + tables + layout

# 4. Publication network figures (PDF / SVG / EPS / PNG)
python code/shared/plot_network_vector.py                       # Flu: Fig 2G (all_edges light)
python code/shared/plot_network_dmf5_light.py                   # DMF5: Fig 3D (dist_le3 light)
```

Optional analyses:

```bash
python code/shared/plot_seqlogo.py Flu                          # Fig 2F sequence logo
python code/shared/plot_seqlogo.py DMF5                         # Fig 3C sequence logo
python code/shared/plot_hamming_distribution.py                 # pairwise Hamming distance histograms
```

**Note on the network pipeline.** The committed `nodes.csv` / `links.csv` are the source of
truth for the publication figures — step 4 alone reproduces Fig 2G / Fig 3D exactly from these
tables. The interactive HTML (step 2) is an intermediate product and is not committed: it uses
a browser-side force-directed layout, so it never matches the fixed vector-figure layout.
Steps 2–3 re-derive the tables from the raw count data and yield the same nodes, edges, and
cluster colors; only the Flu label subsample differs run-to-run (isolated nodes always
labelled; connected nodes sampled by degree, see `DISPLAY_RATIO` in `pyechart_flu.py`).

**Note on label placement in the vector figures.** Node positions are fully deterministic
(fixed seeds in `spring_layout`), but `adjustText` places peptide labels with unseeded
randomness, so re-running the plotting scripts reproduces the same network with slightly
shifted label positions (~1 px jitter).

## Dependencies

**Python** ≥ 3.8:

```bash
pip install -r requirements.txt
```
