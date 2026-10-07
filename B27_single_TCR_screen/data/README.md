# Data

Reference files, pre-selection counts and the processed screen results are shipped with this
repository; only the raw sequencing data stays local. The reference FASTA libraries are zipped
to keep the repository light — unzip them into the same folder before running the pipeline:

```bash
unzip 'data/reference/*.zip' -d data/reference/
```

```
data/
├── reference/                           # zipped FASTA libraries; unzip before use
│   ├── human_proteome_B2705.zip         #   human peptide-DNA library (146 MB -> 29 MB)
│   ├── virome_proteome_B2705.zip        #   virome peptide-DNA library (55 MB -> 11 MB)
│   ├── human_reviewed_unreviewed.zip    #   UniProt human proteome (gene-name lookup)
│   └── 10239_viruses_reviewed.zip       #   UniProt viral proteome (gene-name lookup)
├── round0/
│   ├── V2_B27_human_pep_count.csv       # pre-selection (R0) counts, human library
│   └── V2_B27_virome_pep_count.csv      # pre-selection (R0) counts, virome library
├── results/                             # processed screen results (TCR59 / TCR63)
│   ├── TCR59_R1R2_annotated_with_R0_matched_filtered_with_primer.csv
│   ├── TCR63_R1R2_annotated_with_R0_matched_filtered_with_primer.csv
│   ├── TCR59_R1R2_annotated_with_R0_matched_filtered_with_primer_add_flow.csv
│   └── TCR63_R1R2_annotated_with_R0_matched_filtered_with_primer_add_flow.csv
└── sequencing/                          # raw clean FASTQ per round (not committed)
```

## File formats

### Peptide-DNA library FASTA (`reference/human_proteome_B2705.zip`, `reference/virome_proteome_B2705.zip`)

Zipped FASTA; unzip before use. One entry per designed peptide; the peptide DNA is embedded
between the two constant flanking sequences shared by extraction and annotation:

```
>human+sp_Q67878_DPOL_+GRAVELHNF
CTGTGCTCGCGCTACTCTCTCTTTCTGGCCTGGAGGCT<peptide DNA>GGATGCGGAGGGTCCGGCGGGGGA
```

- Header convention: `<origin>+<UniProt/entry id>_<...>+<peptide AA>`; the gene name is resolved
  from the UniProt proteome FASTAs via the UniProt ID.
- 5' flank: `CTGTGCTCGCGCTACTCTCTCTTTCTGGCCTGGAGGCT`
- 3' flank: `GGATGCGGAGGGTCCGGCGGGGGA`

### UniProt proteome FASTAs (`reference/human_reviewed_unreviewed.zip`, `reference/10239_viruses_reviewed.zip`)

Zipped UniProt proteomes, used by `03_annotate_origin.py --proteome_files` to map UniProt IDs
to gene names. Standard UniProt headers:

```
>sp|P13497|BMP1_HUMAN Bone morphogenetic protein 1 OS=Homo sapiens ... GN=BMP1
```

### Round-0 count tables (`round0/*.csv`)

Pre-selection (R0) sequencing counts per peptide:

| Peptide | Read Count |
|---|---|
| GRAVELHNF | 304 |
| ... | ... |

### Processed results (`results/*.csv`)

One row per candidate peptide, as produced by the pipeline
(see [scripts/peptide_screen/README.md](../scripts/peptide_screen/README.md) for the full column
dictionary). The `*_add_flow.csv` variants extend the final tables with flow-cytometry
validation:

| Column | Description |
|---|---|
| `add_flow_single_copy` | single-copy percentage measured in the flow validation assay for the corresponding TCR-peptide pair |

### Workflow sequencing input (`sequencing/`, not committed)

Adapter-trimmed clean FASTQ per enrichment round (`.clean.fq` or `.clean.fq.gz`). Each read must
contain `CTGGAGGCT<peptide DNA>GGATGC`; the 27 nt / 30 nt inserts are reported as 9mer / 10mer.
