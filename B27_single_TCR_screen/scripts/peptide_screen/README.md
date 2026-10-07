# Workflow 1 — `peptide_screen` (Illumina peptide-only screen, TCR59 / TCR63)

End-to-end processing of a **SEEKER**-style peptide-only single-TCR antigen screen: from
adapter-trimmed Illumina FASTQ to the final primer-ready candidate peptide list
(`<sample>_R1R2_annotated_with_R0_matched_filtered_with_primer.csv`).

## Pipeline

```
R2 clean FASTQ ──► 01 extract ──► 02 combine (R1+R2) ──► 03 annotate origin
                                                            │
                   final with_primer.csv ◄── 06 primers ◄── 05 filter ◄── 04 add R0
```

| Step | Script | Input -> Output |
|---|---|---|
| 1 | `01_extract_peptide_counts.py` | clean FASTQ -> `fullresult_*.csv` (nt, count, aa, length) |
| 2 | `02_combine_rounds.py` | R1 + R2 fullresult -> merged CSV (+ percentages, `r2_r1_ratio`) |
| 3 | `03_annotate_origin.py` | merged CSV -> + `fasta_name`, `origin`, `Gene` |
| 4 | `04_add_round0.py` | -> + `count_round0`, `count_pct_round0`, `r1_r0_ratio` |
| 5 | `05_filter_matched.py` | keep `origin` in {human, virus} |
| 6 | `06_add_primers.py` | sort by `count_round2` desc, renumber 1..N, add primer-F/R |

## Run

Copy the example config and edit paths + sample sections:

```bash
cp configs/peptide_screen.yml.example configs/my_screen.yml
python scripts/peptide_screen/run_peptide_screen.py --config configs/my_screen.yml
```

The driver runs all six steps per sample, streams every command into `logs/`, prints step-wise
progress with estimated remaining time, and writes a Markdown processing report next to the
results. Multiple samples can be listed under `samples:` in the config and are processed in
sequence.

Each step can also be run standalone, e.g.:

```bash
python scripts/peptide_screen/01_extract_peptide_counts.py R2_1.clean.fq TCR59-R2_plasmid_library.csv results/
python scripts/peptide_screen/02_combine_rounds.py fullresult_R1.csv fullresult_R2.csv TCR59_R1R2_combined.csv TCR59
python scripts/peptide_screen/06_add_primers.py TCR59_R1R2_annotated_with_R0_matched_filtered.csv TCR59_with_primer.csv
```

## Final output columns

| Column | Description |
|---|---|
| `nt`, `aa`, `length` | peptide DNA / amino-acid sequence / length |
| `count_pct_round0/1/2`, `count_round0/1/2` | counts and percentages in the pre-selection (R0) and enrichment rounds (R1, R2) |
| `r1_r0_ratio`, `r2_r1_ratio` | enrichment ratios between consecutive rounds |
| `fasta_name`, `origin`, `Gene` | matching peptide-library entry, origin (human/virus), gene name |
| `peptide_number` | rank after sorting by `count_round2` (descending) |
| `primer_F_name`, `primer-F` | forward primer name / sequence (`AGGCT` + peptide DNA) |
| `primer_R_name`, `primer-R` | reverse primer name / sequence (`ATCC` + reverse complement + `A`) |

## Processing details

- **Peptide extraction** scans reads for `CTGGAGGCT[peptide]GGATGC`, counts exact k-mers and
  translates them (9mer/10mer subsets are emitted separately). Files > 200M reads are processed
  in batches.
- **Annotation** matches the peptide DNA against the human and virome peptide-DNA libraries via a
  pre-built index (exact match) and maps UniProt IDs to gene names from the UniProt proteome
  FASTAs.
- **Round-0 merge** attaches pre-selection counts and computes `r1_r0_ratio`.
- **Filtering** keeps only peptides whose origin was resolved (human/virus).
- **Primer design** sorts by `count_round2` (descending, stable), renumbers peptides 1..N and
  derives synthesis oligos.
