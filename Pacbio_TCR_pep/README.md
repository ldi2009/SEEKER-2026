# Workflow — PacBio TCR+peptide screen (`Pacbio_TCR_pep`)

End-to-end processing of a **long-read (PacBio HiFi) TCR+peptide tandem-construct screen**: from
per-round Revio BAM files to the final primer-ready TCR-peptide candidate list
(`analysis/filtered_data_clean_with_primer.csv`). Example dataset:  **2704**, V2 sample,
enrichment rounds **R2 + R3**.

Because TCR and peptide are captured on the same molecule, every candidate is a **TCR-peptide
pair** (unlike the Illumina `peptide_screen` workflow, which is peptide-only).

## Pipeline

```
per-round BAM ──► samtools bam2fq ──► 01_extract_tcr_peptide.py (regex anchors + minimap2)
                                          │  per-read: TCRA/TCRB/TCRfull/peptide
                                          ▼
                                     02_final_qc.py (per round)
                                          │  full-M CIGAR + peptide QC,
                                          │  mismatch/missing detection, 1/10 threshold
                                          ▼
              final_QC.csv (R2) ──► 03_solo_analysis.py ──► final_QC.csv (R3)
                                          │  two-round merge -> origin lookup ->
                                          │  percent + fold change -> per-peptide dedup (1/20)
                                          │  -> 7-step clean filter -> primer design
                                          ▼
                       analysis/filtered_data_clean_with_primer.csv
```

| Step | Script                      | Input -> Output                                                                                 |
| ---- | --------------------------- | ----------------------------------------------------------------------------------------------- |
| 1    | `01_extract_tcr_peptide.py` | HiFi FASTQ -> `<prefix>_result.csv` / `_result_qc.csv` (+ combination counts)                   |
| 2    | `02_final_qc.py`            | per-round result.csv -> `final_QC.csv`, `questionable.csv`, mismatch/missing tables             |
| 3    | `03_solo_analysis.py`       | two `final_QC.csv` -> merged/origin/fold-change tables -> `filtered_data_clean_with_primer.csv` |

`data/reference/101TCRinfo_processed_with_counts.csv` is the processed
101-TCR reference panel table: per-TCR TRA/TRB V-J genes and CDR3 (aa + nt,
reference-encoded chain set), patient labels (`#1`-`#4`), and per-clone
cytotoxic-state cell counts (`total_cells`, `GZMBneg_GZMKneg`, `GZMBpos`,
`GZMKpos`, `GZMBpos_GZMKpos`).

## Run

```bash
cp configs/pacbio_screen.yml.example configs/my_pacbio.yml   # edit paths + rounds
python scripts/run_pacbio_screen.py --config configs/my_pacbio.yml
```

The driver converts BAM -> FASTQ (`samtools bam2fq`), runs steps 1-2 per round, then step 3,
streaming every command into `logs/` with step-wise progress and estimated remaining time, and
finishes with a Markdown processing report.

## Final output columns

| Column                                            | Description                                                |
| ------------------------------------------------- | ---------------------------------------------------------- |
| `peptide`, `peptideDNA`                           | peptide amino-acid / DNA sequence                          |
| `TCRA_align`, `TCRB_align`, `TCRfull_align`       | best TCR reference hit per insert (e.g. `TCR_038`)         |
| `count_round2`, `count_round3`                    | merged read counts in the two rounds (R2, R3)              |
| `percent_round2`, `percent_round3`, `fold_change` | per-round percentages and enrichment                       |
| `fasta_name`, `origin`, `Gene`                    | peptide-library entry, origin (human/virus), gene name     |
| `forward_primer` / `reverse_primer`               | synthesis oligos (`AGGCT`+peptideDNA / `ATCC`+revcomp+`A`) |
| `primer_number`, `forward_name`, `reverse_name`   | rank (1..N) and oligo names                                |

## The 7-step clean filter

1. remove control peptides containing `DMF5`
2. remove control peptides containing `TCR_074 (optional for removing high tonic signal TCR)`
3. remove rows with an empty TCRA/TCRB chain
4. remove peptides without a resolved origin (`INFO=MISSING` / `INFO=UNKNOWN`)
5. keep peptides of length >= `min_peptide_len` (default 8; an optional `special_peptides`
   whitelist can bypass this length rule but ships empty)
6. keep `fold_change` within a configurable window (2704-V2 example: 0.1 - 60)
7. keep only pairs whose TCR is in the allowed set (2704: TCR\_001..TCR\_049 + TCR\_101)

All thresholds are configurable via the `analysis:` block of the config.
