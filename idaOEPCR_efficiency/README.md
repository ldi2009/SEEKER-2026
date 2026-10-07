# idaOEPCR efficiency test

**idaOEPCR** = *in droplet dual-asymmetric overlap extension PCR*.

Analysis code supporting the idaOEPCR validation panels of **Figure 4** in
*SEEKER: a genome-scale library-on-library screening platform for deciphering
T cell recognition of antigen*.

The scripts take the BLAST-annotated idaOEPCR sequencing tables as input, call
unique peptides, and regenerate the four Figure 4 idaOEPCR panels together with
the processed result tables. Raw sequencing reads are **not** included here
(see [Data availability](#9-data-availability)).

---

## 1. Repository layout

```
.
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── scripts/
│   ├── 01_idaOEPCR_unique_peptides.py
│   ├── 02_idaOEPCR_group_binning_cumulative.py
│   └── 03_figure4_idaOEPCR_ranking_distribution.py
├── data/                                 # 8 BLAST-annotated sample tables
│   ├── idaOEPCR-optimize-2_Not_1_top_10000_blast_sequences.csv
│   ├── idaOEPCR-optimize-2_Not_2_top_10000_blast_sequences.csv
│   ├── idaOEPCR-optimize-2_Not_3_top_10000_blast_sequences.csv
│   ├── idaOEPCR-optimize-2_Not_4_top_10000_blast_sequences.csv
│   ├── idaOEPCR-optimize-2_standard_1_top_10000_blast_sequences.csv
│   ├── idaOEPCR-optimize-2_standard_2_top_10000_blast_sequences.csv
│   ├── idaOEPCR-optimize-2_standard_3_top_10000_blast_sequences.csv
│   └── idaOEPCR-optimize-2_standard_4_top_10000_blast_sequences.csv
└── results/
    ├── Not_ranking_distribution.pdf
    ├── standard_ranking_distribution.pdf
    ├── Not_new_combined_distribution.png
    ├── standard_new_combined_distribution.png
    ├── hamming_summary_results.csv       # per-sample unique-peptide counts
    └── per_sample/                        # per-sample unique-peptide tables
        ├── idaOEPCR-optimize-2_Not_1_unique_peptides.csv
        ├── idaOEPCR-optimize-2_Not_2_unique_peptides.csv
        ├── idaOEPCR-optimize-2_Not_3_unique_peptides.csv
        ├── idaOEPCR-optimize-2_Not_4_unique_peptides.csv
        ├── idaOEPCR-optimize-2_standard_1_unique_peptides.csv
        ├── idaOEPCR-optimize-2_standard_2_unique_peptides.csv
        ├── idaOEPCR-optimize-2_standard_3_unique_peptides.csv
        └── idaOEPCR-optimize-2_standard_4_unique_peptides.csv
```

---

## 2. Scripts

| Script | Purpose | Input | Output |
| --- | --- | --- | --- |
| `01_idaOEPCR_unique_peptides.py` | Unique-peptide calling and per-sample tables | `data/*.csv` | `results/hamming_summary_results.csv`, `results/per_sample/*` |
| `02_idaOEPCR_group_binning_cumulative.py` | Group-level log10 binning and cumulative curves of unique peptides | `data/*.csv` | group bin/cumulative plots (not shipped in `results/`) |
| `03_figure4_idaOEPCR_ranking_distribution.py` | Final Figure 4 ranking and combined-distribution panels | `data/*.csv` | the four panels in `results/` |

All scripts take explicit `--input-dir` / `--output-dir` paths and can be run
from any directory.

---

## 3. Input data

### 3.1 BLAST-annotated sample tables (`*_top_10000_blast_sequences.csv`)

One table per idaOEPCR sample, containing the 10,000 most abundant sequences.
The eight tables used to generate the Figure 4 panels — four `Not` replicates
and four `standard` replicates — are provided in [`data/`](data).

| Column | Description |
| --- | --- |
| `peptideDNA` | Nucleotide sequence of the peptide-encoding region |
| `seq-TRAV-head` | Paired TRAV head sequence |
| `Translated_Peptide` | Translated peptide sequence |
| `Count` | Read count |
| `best_alpha_hit` | Best-matching alpha chain |
| `best_alpha_pident` | Percent identity of the best alpha-chain hit |
| `best_alpha_length` | Length of the best alpha-chain hit |

---

## 4. Core definition — unique peptide

The criterion is evaluated on the nucleotide sequence `peptideDNA`, using the
**Hamming distance** between equally long sequences, i.e. the number of
positions at which two sequences differ.

| Term | Definition |
| --- | --- |
| **Unique peptide count** | Number of peptides retained as *unique* in a sample. A peptide is unique when **no more abundant peptide (higher `Count`) lies within a Hamming distance of 2** — that is, it has no more abundant neighbour within two nucleotides. |

Candidate pool: the top 1,000 peptides per sample (by `Count`, descending).

---

## 5. Algorithm — Hamming distance and unique-peptide calling (script `01`)

This section documents the exact logic of script
`01_idaOEPCR_unique_peptides.py`. All comparisons are made on the nucleotide
sequence in the `peptideDNA` column, so a distance of 1 always means a
single-nucleotide difference.

### 5.1 Hamming distance

Two sequences of equal length are compared position by position; sequences of
different length are declared non-comparable and returned as `inf`, so they
can never be neighbours:

```python
def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return float("inf")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
```

The return value is therefore the number of positions at which the two
sequences differ (`0` for identical sequences, `1` for a single-nucleotide
difference, and so on).

### 5.2 Unique-peptide calling

The top 1,000 peptides (by `Count`) are scanned in descending abundance order.
A peptide is retained as **unique** when it lies at a Hamming distance `> 2`
from every peptide already retained — i.e. it has no more abundant neighbour
within two nucleotides:

```python
def unique_peptide_rows(df):
    top = df.nlargest(1000, "Count").reset_index(drop=True)
    kept_seqs = []
    kept_rows = []
    for _, row in top.iterrows():
        seq = str(row["peptideDNA"])
        if all(hamming_distance(seq, kept) > 2 for kept in kept_seqs):
            kept_seqs.append(seq)
            kept_rows.append(row)
    return pd.DataFrame(kept_rows)
```

Points that define the exact behaviour:

* **Candidate pool.** Only the 1,000 most abundant peptides of the sample are
  considered.
* **Descending abundance.** Peptides are visited from highest to lowest
  `Count`, so a peptide is judged against every more abundant peptide that was
  retained before it.
* **Distance strictly greater than 2.** The test is `> 2`, so a peptide is
  rejected only when a retained (more abundant) peptide differs by **one or
  two** nucleotides — equivalently, a peptide must differ from every more
  abundant neighbour by at least three positions.
* **Aggregation.** The number of retained peptides is the sample's **unique
  peptide count**, written to `hamming_summary_results.csv`. The retained
  peptides themselves are written to `per_sample/<sample>_unique_peptides.csv`.

---

## 6. Usage

```bash
# 01 — unique peptides + per-sample tables
python scripts/01_idaOEPCR_unique_peptides.py \
    --input-dir data \
    --output-dir results

# 02 — group-level binning and cumulative curves (optional)
python scripts/02_idaOEPCR_group_binning_cumulative.py \
    --input-dir data \
    --output-dir results/bin_cumulative

# 03 — Figure 4 ranking and combined-distribution panels
python scripts/03_figure4_idaOEPCR_ranking_distribution.py \
    --input-dir data \
    --output-dir results \
    --groups Not standard
```

---

## 7. Outputs

| Script | Output | Content |
| --- | --- | --- |
| `01` | `results/hamming_summary_results.csv` | Per sample: `unique_peptide_count` |
| `01` | `results/per_sample/<sample>_unique_peptides.csv` | Retained unique peptides per sample: `sample`, `peptideDNA`, `best_alpha_hit`, `Count` |
| `02` | `*_bin_cumulative_combined.png`, `group_unique_data.csv` | Per-group and combined bin/cumulative plots (not shipped) |
| `03` | `{group}_ranking_distribution.pdf` | Four stacked replicate panels, Rank vs log2(Count); unique peptides blue, remaining peptides orange |
| `03` | `{group}_new_combined_distribution.png` | The four replicates overlaid as a log10(read count) histogram of unique-peptide counts |

The four Figure 4 panels shipped in `results/` are:

| File | Content |
| --- | --- |
| `Not_ranking_distribution.pdf` | `Not` group, four replicate rank–abundance panels |
| `standard_ranking_distribution.pdf` | `standard` group, four replicate rank–abundance panels |
| `Not_new_combined_distribution.png` | `Not` group, overlaid unique-peptide count distribution |
| `standard_new_combined_distribution.png` | `standard` group, overlaid unique-peptide count distribution |

---

## 8. Requirements

Python 3.8+ with:

```
pandas
numpy
matplotlib
seaborn
```

Install with:

```bash
pip install -r requirements.txt
```

---

## 9. Data availability

Raw sequencing reads (Illumina PE150 for the SCD libraries, PacBio Revio for
the linked TCR–SCD libraries, and 10x Genomics scRNA-seq / V(D)J data) are
deposited in public repositories; accession codes are provided in the
manuscript. Human-derived single-cell data are deposited under controlled
access.

Only the eight BLAST-annotated idaOEPCR sample tables needed to regenerate the
Figure 4 panels are distributed in this repository, under [`data/`](data),
together with the processed result tables and figures. All other raw tables
are available from the corresponding author on request.

---

## 10. Notes

* `01` uses an explicit per-row Hamming-distance scan over the top-1,000
  peptides; it writes the summary table and the per-sample tables only — the
  figures are produced by `03`.
* `03` uses an exact but faster equivalent of the per-row unique-peptide call:
  for every row only the Hamming-distance-2 neighbourhood of its sequence is
  intersected with the set of sequences already seen at higher counts. The
  calling rule is the same as in `01` (no more abundant peptide within a
  Hamming distance of 2).
* `02` is retained as supporting code for the group-level binning / cumulative
  curves; its plots are not part of the shipped `results/`.
* Plot labels and titles in the scripts are in English and the source code
  contains no Chinese characters. The underlying analysis (thresholds,
  binning, Hamming-distance criteria) is unchanged from the working scripts.

---

## 11. License

Released under the MIT License — see [LICENSE](LICENSE).