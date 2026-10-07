# SEEKER-2026

This is the repository for **SEEKER**, a library-on-library screening platform for deciphering T cell recognition of antigen.

Each project lives in its own folder with its own README. Currently available:

| Project | Description | README |
|---|---|---|
| [`B27_single_TCR_screen`](B27_single_TCR_screen/) | HLA-B\*27:05 single-TCR antigen screens with SEEKER coupled to deep sequencing (Illumina peptide-only; TCR59/TCR63) incl. virus–self molecular-mimicry analysis | [B27_single_TCR_screen/README.md](B27_single_TCR_screen/README.md) |
| [`Pacbio_TCR_pep`](Pacbio_TCR_pep/) | PacBio HiFi TCR+peptide tandem-construct screen (patients 2704/2705 V2, rounds R2+R3) | [Pacbio_TCR_pep/README.md](Pacbio_TCR_pep/README.md) |
| [`SEEKER_single_TCR_DMF5_Flu_screen`](SEEKER_single_TCR_DMF5_Flu_screen/) | Single-TCR 9-mer peptide screens (DMF5 NNK and Flu-M1 mutagenesis libraries, Illumina rounds R1–R3) with peptide similarity networks | [SEEKER_single_TCR_DMF5_Flu_screen/README.md](SEEKER_single_TCR_DMF5_Flu_screen/README.md) |
| [`idaOEPCR_efficiency`](idaOEPCR_efficiency/) | idaOEPCR (in-droplet dual-asymmetric overlap-extension PCR) library-construction efficiency test (Figure 4 validation panels) | [idaOEPCR_efficiency/README.md](idaOEPCR_efficiency/README.md) |

> More projects will be added as they are processed.

## Data availability

All processed/derived data tables supporting the manuscript (screening result tables, R0
counts, mimicry-pair lists, network tables, idaOEPCR BLAST tables, the 101-TCR reference panel
with per-clone cytotoxic-state counts, and the TCR reference FASTA) are deposited on **Dryad**
(DOI: to be assigned on submission) — see the package `README.md` there for a file-by-file
description. Raw sequencing reads are deposited separately (accessions to be added upon
acceptance). Small copies of the key result tables are also kept in this repository for
convenience, as noted in each project's README.

## License

Code: MIT — see [LICENSE](LICENSE). Data on Dryad: CC0 1.0.
