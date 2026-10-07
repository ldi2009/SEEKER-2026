# Molecular mimicry figures (TCR59 / TCR63)

This folder reproduces the two molecular-mimicry figures of the TCR59/TCR63
virus–self peptide cross-reactivity analysis:

- **Fig11** (`Fig11_Mimicry_Top_Pairs_Selected`) — sequence-level view of the
  six highest-confidence virus–self mimicry pairs (2 × TCR59, 4 × TCR63).
- **Fig14** (`Fig14_Mimicry_Network_Large`) — heat map of mimicry pairs
  aggregated by viral family (rows) and self-antigen gene (columns).

## Contents

| File | Role |
|---|---|
| `mimicry_analysis.py` | Derives the virus–self mimicry pairs (outputs `mimicry_pairs.json`) |
| `gen_fig11_fig14_svg.py` | Plotting script (matplotlib) |
| `mimicry_pairs.json` | 437 virus–self peptide mimicry pairs (output of `mimicry_analysis.py`) |
| `virus_peptides_annotated_v4.csv` | Virus peptides with UniProt annotation and activation values |
| `data/TCR59_peptide_validation_summary.csv` | TCR59 peptide validation summary (self-antigen side) |
| `data/TCR63_peptide_validation_summary.csv` | TCR63 peptide validation summary (self-antigen side) |
| `Fig11_Mimicry_Top_Pairs_Selected.png / .svg` | Figure 11 (raster + vector) |
| `Fig14_Mimicry_Network_Large.png / .svg` | Figure 14 (raster + vector) |

## Reproduce

```bash
pip install matplotlib numpy
python mimicry_analysis.py        # regenerates mimicry_pairs.json
python gen_fig11_fig14_svg.py     # regenerates both figures
```

All paths resolve relative to the script location, so the folder runs as-is
wherever it is cloned. Both a 300-dpi PNG and a vector SVG are written next
to the script; SVG text is kept as editable text (`svg.fonttype: none`).

## Input data

**`mimicry_pairs.json`** — one record per virus–self peptide pair, with
fields `virus_seq`, `human_seq`, `virus_gene`, `human_gene`,
`virus_organism`, `virus_family`, `virus_tcr`, `human_tcr`, `similarity`,
`jaccard`, `combined`, `common_motif`, `shared_kmers`, and the activation
values `virus_single` / `human_single`.

A pair was called when both peptides are activation-positive in addflow
validation (> 5%), differ in length by ≤ 6 residues, and either the combined
score exceeds 0.35 — combined = 0.6 × Levenshtein similarity (1 − edit
distance / longer peptide length) + 0.4 × 3-mer Jaccard index — or the two
peptides share a contiguous motif of ≥ 4 amino acids. In total 437 pairs
qualify: 60 with TCR59 virus peptides and 377 with TCR63 virus peptides.
`mimicry_analysis.py` recomputes this from the annotated virus CSV and the
two validation summaries in `data/` (activation flags: TCR59 single- or
high-copy addflow > 5% / > 30%; TCR63 summary and virus-peptide annotation
single-copy addflow, as TCR63 was validated without high-copy data).

**`virus_peptides_annotated_v4.csv`** — per virus peptide: source TCR,
sequence, UniProt AC / entry, gene name, organism, viral family, host and
human-pathogen annotation, plus the activation values `single_pct`,
`positive`, `strong`. The annotation was fetched from the UniProt REST API;
the figure script also uses this file to look up the full virus names
displayed in Fig11 (gene name taken from `gene` when it carries the
`GENE_SPECIES` suffix, otherwise from `uniprot_entry_id`).

## How each figure is built

**Fig11 (selected pairs).** From the 437 pairs, the script keeps pairs of
common human-infecting viruses (Rotavirus A, HSV-1/2, HCMV, HPV-30, Torque
teno virus, Hepatitis delta virus, Echovirus 12, Yellow fever virus), with
one pair per virus species and both peptides originating from the same TCR
(cross-TCR pairs excluded). The six displayed pairs are the curated
selection defined in `SELECTED_KEYS`:

| TCR | Virus peptide | Virus (gene) | Self-antigen gene |
|---|---|---|---|
| TCR59 | FRNSMHML | VP2 (Rotavirus A) | MYO5A |
| TCR59 | ARRALEASVR | HEPA (Human herpesvirus 2) | PRODH |
| TCR63 | RKARKTIKK | SHDAG (Hepatitis delta virus) | ADNP |
| TCR63 | RRRRKAVRR | CAPSD (Torque teno virus) | CACNA1H |
| TCR63 | RRRRKHVPYFL | VL2 (Human papillomavirus 30) | MRPS25 |
| TCR63 | HRRRKHLAVQR | UL87 (Human cytomegalovirus) | CYP26C1 |

Each row block aligns the virus peptide above the self peptide. Residues
that match positionally or fall inside the shared motif are filled with the
TCR color (mint #66C2A5 = TCR59, rose #FC8D62 = TCR63, white text);
mismatches are light grey with dark text. The virus is labeled as
`GENE (full virus name)` in blue italics; the self-antigen gene name is
colored by TCR.

**Fig14 (network heat map).** For each TCR, the script selects the five
self-antigen genes with the most mimicry pairs (TCR59: ZNF566, PRODH,
KIF2B, PRADC1, TRMT2B; TCR63: BMP1, CYP26C1, NOP2, CSNK1G2, PCDH15) and
counts all pairs between these genes and the 15 viral families that have at
least one pair among them (families ordered by total pair count; families
with no pair are omitted). Cells show pair counts on a white→orange
gradient (counts ≥ 5 switch to white text). Column labels are colored by the
TCR in which the self antigen was identified, row labels by the TCR origin
of the virus peptide (green TCR59, orange TCR63, dark grey shared by both),
and a black divider separates the TCR59 from the TCR63 gene block. Origin
legends sit below the heat map.

## Style

Arial throughout; TCR59 = #66C2A5 and TCR63 = #FC8D62 (ColorBrewer Set2);
oversized fonts and thick white grid lines for legibility in print.

> Note: the two scripts here cover the mimicry step only. The upstream
> peptide validation summaries in `data/` and the UniProt annotation of the
> virus peptides are produced by the screening pipeline maintained elsewhere
> in this repository.
