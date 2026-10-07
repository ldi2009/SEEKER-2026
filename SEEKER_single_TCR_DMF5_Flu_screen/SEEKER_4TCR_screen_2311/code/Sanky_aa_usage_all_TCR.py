# -*- coding: utf-8 -*-
"""4-TCR lined heatmap + frequency-matrix generator (seqlogo criteria, 2026-10-07).

Peptide source: data/Selected_Peptides.csv (306 peptides, built by
make_selected_peptides.py), identical to the four results/seqlogos/ PDFs:
  DMF5 147 / 1G4 105 / EBV 48 / Flu 6 — 306 in total, no overlap between TCRs.

Usage:
  python Sanky_aa_usage_all_TCR.py            # output into results/ of this package
  python Sanky_aa_usage_all_TCR.py <dir>      # output to a given directory
Writes: lined_{TCR}_4TCR_heatmap.pdf + frequency_{TCR}_4TCR_heatmap.csv
"""
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE, "data")
SELECTED_CSV = os.path.join(DATA_DIR, "Selected_Peptides.csv")

AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
               'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def load_selected_peptides():
    sel = pd.read_csv(SELECTED_CSV)
    return {t: g["aa"].tolist() for t, g in sel.groupby("tcr_name")}


def draw_peptide_heatmap(peptides, with_lines=False, output_file="heatmap.pdf", matrix_file=None):
    matrix = np.zeros((len(AMINO_ACIDS), len(peptides[0])))
    for i, aa in enumerate(AMINO_ACIDS):
        for j in range(len(peptides[0])):
            matrix[i, j] = sum(p[j] == aa for p in peptides) / len(peptides)

    if matrix_file:
        pd.DataFrame(matrix, index=AMINO_ACIDS,
                     columns=[f"Pos_{i+1}" for i in range(len(peptides[0]))]).to_csv(matrix_file)

    fig, ax = plt.subplots()
    im = ax.imshow(matrix, cmap=plt.get_cmap("Blues"))
    cbar = ax.figure.colorbar(im, ax=ax, ticks=np.arange(0, 1.01, 0.1))
    cbar.ax.set_ylabel("Frequency", rotation=-90, va="bottom")
    ax.set_xticks(np.arange(len(peptides[0])))
    ax.set_yticks(np.arange(len(AMINO_ACIDS)))
    ax.set_xticklabels(np.arange(1, len(peptides[0]) + 1))
    ax.set_yticklabels(AMINO_ACIDS)

    if with_lines:
        for p in peptides:
            xs, ys = [], []
            for i, aa in enumerate(p):
                if aa in AMINO_ACIDS:
                    xs.append(i)
                    ys.append(AMINO_ACIDS.index(aa))
            plt.plot(xs, ys, color="black", alpha=0.01)

    plt.xlabel("Peptide position")
    plt.ylabel("Amino acid")
    plt.title("Peptide sequence distribution")
    plt.savefig(output_file, format="pdf")
    plt.close()


if __name__ == "__main__":
    out_root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, "results")
    heat_dir = os.path.join(out_root, "heatmaps")
    os.makedirs(heat_dir, exist_ok=True)
    tcr_peptides = load_selected_peptides()
    for tcr, peptides in tcr_peptides.items():
        print(f"[{tcr}] Selected_Peptides (seqlogo criteria) -> {len(peptides)} peptides")
        draw_peptide_heatmap(peptides, with_lines=True,
                             output_file=os.path.join(heat_dir, f"lined_{tcr}_4TCR_heatmap.pdf"),
                             matrix_file=os.path.join(heat_dir, f"frequency_{tcr}_4TCR_heatmap.csv"))
    print(f"Done. Output root: {out_root}")
