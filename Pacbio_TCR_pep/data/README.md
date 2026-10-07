# Data layout (`Pacbio_TCR_pep`)

This project reads its inputs from `data/` (paths are configured in the YAML config). Large
files stay local — nothing under `data/sequencing/` is committed. Create the layout below or
point the config anywhere you like:

```
data/
├── processed_data/                      # shipped processed results (patient 2704 / 2705)
│   ├── 2704_R23_filtered_data_clean_with_primer.csv
│   └── 2705_R23_filtered_data_clean_with_primer.csv
├── reference/
│   ├── 101TCR_reference.fasta        # TCR reference, included in this repo
│   │                                 #   (TCR_001..TCR_101 + DMF5-TCR_000 control;
│   │                                 #    minimap2 -ax map-hifi target)
│   ├── 101TCRinfo_processed_with_counts.csv
│   │                                 #   one row per TCR: V-J genes + CDR3 (aa/nt) of the
│   │                                 #   reference-encoded chain set, patient label (#1-#4),
│   │                                 #   per-clone cytotoxic-state cell counts
│   ├── human_proteome_B2705.fasta    # human peptide-DNA library (unzip from the repo root
│   │                                 #   B27_single_TCR_screen/data/reference/ or re-download)
│   └── virome_proteome_B2705.fasta   # virome peptide-DNA library
└── sequencing/
    ├── 2704_V2_good_R2/Revio/*/*.bam # per-round HiFi BAM (any glob you configure)
    └── 2704_V2_good_R3/Revio/*/*.bam
```

## Processed results (`processed_data/`)

Final TCR-peptide candidates from the three-source merged R23 analysis
(old + single + new FASTQs per round; FC window 0.1-60, 7-step clean filter), one file per
patient, sorted by `count_round3`:

| File | Rows | Description |
|---|---|---|
| `2704_R23_filtered_data_clean_with_primer.csv` | 59,812 | patient 2704, TCR_001..TCR_049 + TCR_101 |
| `2705_R23_filtered_data_clean_with_primer.csv` | 29,653 | patient 2705, TCR_050..TCR_100 |

Columns: `peptide`, `TCRA_align`/`TCRB_align`/`TCRfull_align`, `peptideDNA`,
`count_round2`/`count_round3` (merged read counts in the two enrichment rounds),
`fasta_name`/`origin`/`Gene`, `percent_round2`/`percent_round3`, `fold_change`,
`forward_primer`/`reverse_primer` (`AGGCT`+peptideDNA / `ATCC`+revcomp+`A`),
`primer_number`, `forward_name`, `reverse_name`.

## Notes

- **`101TCR_reference.fasta` is included** (`data/reference/`): 101 patient TCRs
  (`TCR_001`...`TCR_101`) plus the `DMF5-TCR_000` control entry. HiFi reads carrying
  TCRA/TCRB/TCRfull inserts are aligned with `minimap2 -ax map-hifi`; the best hit names the
  TCR. Reads assigned to the DMF5 control are removed later by the clean filter
  (`remove_peptide_keywords` / empty-chain / mismatch QC).

- The peptide-DNA libraries and UniProt proteomes used for gene-name lookup are the same files
  used by `B27_single_TCR_screen`; zipped copies ship in that project under
  [`B27_single_TCR_screen/data/reference/`](../B27_single_TCR_screen/data/reference/) — unzip
  them and point the config at the FASTAs, or reuse an existing local copy.
- TCR reference: one FASTA entry per TCR (`TCR_001`...). HiFi reads carrying TCRA/TCRB/TCRfull
  inserts are aligned with `minimap2 -ax map-hifi`; the best hit names the TCR.
- Pipeline outputs land in `analysis/` (final tables) and `output/` (per-round artifacts + logs +
  processing report) — both are git-ignored when produced inside the repository.
