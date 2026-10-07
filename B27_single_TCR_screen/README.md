# B27\_single\_TCR\_screen


End-to-end processing pipelines for **HLA-B\*27:05 single-TCR antigen screens** performed with **SEEKER** coupled to deep sequencing: from raw sequencing data to the final primer-ready candidate peptide list.

This repository hosts one workflow per analysis stage / dataset, each in its own folder with
its own README. Currently available:

| Workflow | Data | Modality | Final output |
|---|---|---|---|
| `peptide_screen` (TCR59 / TCR63) | [scripts/peptide_screen](scripts/peptide_screen/README.md) | Illumina (short-read) | `<sample>_R1R2_annotated_with_R0_matched_filtered_with_primer.csv` |
| `virus_mimicry_analysis` (TCR59 / TCR63) | [virus_mimicry_analysis](virus_mimicry_analysis/README.md) | flow-validation follow-up | Fig11 top mimicry pairs + Fig14 virus-family × self-gene heat map (437 virus↔self pairs) |
| `mimicry_network` (TCR59 / TCR63) | [scripts/mimicry_network](scripts/mimicry_network/README.md) | figure assembly | published mimicry-network figures (PDF/SVG/EPS/PNG) + interactive label editor |

> More datasets / workflows will be added as they are processed.

The workflow reproduces, step by step, the analysis described in the manuscript:

> Single-TCR screening against a combined human proteome / virome peptide library (HLA-B\*27:05),
> followed by two-round enrichment analysis, reference annotation, and oligonucleotide primer
> design for candidate peptide synthesis.

***

## Repository layout

```
SEEKER-2026/
├── configs/
│   └── peptide_screen.yml.example        # peptide_screen config (TCR59/TCR63 example values)
├── scripts/
│   ├── peptide_screen/                   # Illumina peptide-only screen (see its README)
│   │   ├── README.md                     #   workflow documentation
│   │   ├── 01_extract_peptide_counts.py  #   FASTQ -> fullresult CSV (nt, count, aa, length)
│   │   ├── 02_combine_rounds.py          #   merge two enrichment rounds (+ %, r2/r1 ratio)
│   │   ├── 03_annotate_origin.py         #   match peptide DNA to human/virus peptide library
│   │   ├── 04_add_round0.py              #   merge pre-selection (R0) counts
│   │   ├── 05_filter_matched.py          #   keep peptides with a known origin (human/virus)
│   │   ├── 06_add_primers.py             #   sort by R2 count, number peptides, design primers
│   │   └── run_peptide_screen.py         #   config-driven driver (progress/ETA + logs + report)
│   └── mimicry_network/                  # mimicry network figures (see its README)
│       ├── README.md                     #   workflow documentation
│       ├── plot_mimicry_red.py           #   the two published figures
│       └── gen_label_editor.py           #   interactive label editor HTML
├── data/                                 # reference libraries, R0 counts, processed results
│   ├── reference/*.zip                   #   peptide-DNA + UniProt proteome libraries (unzip first)
│   ├── round0/*.csv                      #   pre-selection (R0) read counts
│   ├── results/*.csv                     #   processed TCR59/TCR63 results incl. flow validation
│   ├── mimicry_network/                  #   curated node layouts, render params,
│   │                                     #   hand-adjusted label overrides,
│   │                                     #   UniProt current entry names
│   └── sequencing/                       #   raw clean FASTQ (not committed)
├── figures/
│   └── mimicry_network/                  #   published figures + label editor HTML
├── virus_mimicry_analysis/               # molecular-mimicry analysis (see its README)
│   ├── README.md
│   ├── mimicry_analysis.py               #   derives the 437 virus↔self mimicry pairs
│   ├── gen_fig11_fig14_svg.py            #   Fig11 (top pairs) + Fig14 (family × gene heat map)
│   ├── mimicry_pairs.json                #   mimicry pairs (output)
│   ├── virus_peptides_annotated_v4.csv   #   annotated virus peptides + activation values
│   └── data/                             #   per-TCR validation summaries
├── environment.yml
└── README.md
```

***

## Installation

```bash
git clone https://github.com/ldi2009/SEEKER-2026.git
cd SEEKER-2026

# conda environment (python + pandas + biopython)
conda env create -f environment.yml
conda activate b27-tcr-screen
```

***

## Data

The peptide-DNA reference libraries, the pre-selection (R0) count tables and the processed
TCR59/TCR63 results (including the flow-validation `*_add_flow.csv` tables) are included under
[data/](data/README.md), together with the curated network layouts and label overrides under
[data/mimicry_network/](data/mimicry_network/). Only the raw clean FASTQ files stay local —
point the config's `samples:` section at them to re-run the pipeline.

***

## Usage

See the workflow READMEs: [scripts/peptide_screen/README.md](scripts/peptide_screen/README.md)
(processing) and [scripts/mimicry_network/README.md](scripts/mimicry_network/README.md)
(publication figures).

In short:

```bash
cp configs/peptide_screen.yml.example configs/my_screen.yml   # edit paths + sample section(s)
python scripts/peptide_screen/run_peptide_screen.py --config configs/my_screen.yml

python scripts/mimicry_network/plot_mimicry_red.py              # the two published figures
python scripts/mimicry_network/gen_label_editor.py              # interactive label editor
```

***

## Citation

If you use this code, please cite the manuscript (citation to be added upon publication).

## Contact

Questions: open an issue in this repository.
