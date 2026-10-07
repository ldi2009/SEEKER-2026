# -*- coding: utf-8 -*-
"""
Sequence logo (logomaker frequency-matrix version, published Fig 2F / Fig 3C)
Styled after the published version (R ggseqlogo: method=probability, col_scheme=chemistry):
chemistry colors / Probability axis / x axis 1-N ticks / bottom legend

Usage (from code/shared/):
    python plot_seqlogo.py Flu      # Flu mutagenesis screen: top 200 rows of the combined table
    python plot_seqlogo.py DMF5     # DMF5 single screen: all peptides in results/dmf5_network/selected_peptides.csv

Custom input:
    python plot_seqlogo.py --csv <path> [--column aa] [--top 200] \
        [--out-dir <dir>] [--out-prefix seqlogo] [--title "Title"]

Dependencies: pandas, matplotlib, logomaker
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import logomaker
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

PRESETS = {
    "Flu": dict(
        csv=os.path.join(REPO, "data", "top10000_combine_9mer_2308-Flu-mut-SCD-round123.csv"),
        column="aa", top=200,
        out_dir=os.path.join(REPO, "results", "flu_mutation_screen"),
        prefix="bold_top200_9mer_Flu_mutation",
    ),
    "DMF5": dict(
        csv=os.path.join(REPO, "results", "dmf5_network", "selected_peptides.csv"),
        column="aa", top=None,
        out_dir=os.path.join(REPO, "results", "dmf5_screen"),
        prefix="DMF5_seqlogo",
    ),
}

# Legend for the logomaker 'chemistry' color scheme (matches published ggseqlogo chemistry)
CHEMISTRY_LEGEND = [
    ("#008000", "Polar (G S T Y C)"),
    ("#800080", "Neutral (Q N)"),
    ("#0000ff", "Basic (K R H)"),
    ("#ff0000", "Acidic (D E)"),
    ("#000000", "Hydrophobic (A V L I P W F M)"),
]


def main():
    ap = argparse.ArgumentParser(
        description="Sequence logo from peptide CSV (frequency matrix, logomaker)")
    ap.add_argument("screen", nargs="?", choices=sorted(PRESETS),
                    help="preset: Flu / DMF5")
    ap.add_argument("--csv", help="input CSV path (overrides the preset)")
    ap.add_argument("--column", default="aa", help="sequence column name (default aa)")
    ap.add_argument("--top", type=int, default=None,
                    help="use only the first N sequences (default all)")
    ap.add_argument("--out-dir", default=None, help="output directory (overrides the preset)")
    ap.add_argument("--out-prefix", default=None, help="output filename prefix (overrides the preset)")
    ap.add_argument("--title", default=None, help="figure title (overrides the preset)")
    args = ap.parse_args()

    if args.screen:
        cfg = dict(PRESETS[args.screen])
        cfg["column"] = args.column
        cfg["title"] = None                      # presets have no title, matching the published figure
    else:
        if not args.csv:
            ap.error("provide a preset name (Flu/DMF5) or a --csv path")
        cfg = dict(csv=args.csv, column=args.column, top=None,
                   out_dir=os.path.join(REPO, "results", "seqlogo"),
                   prefix="seqlogo", title="Sequence Logo")
    if args.csv:
        cfg["csv"] = args.csv
    if args.top is not None:
        cfg["top"] = args.top
    if args.out_dir:
        cfg["out_dir"] = args.out_dir
    if args.out_prefix:
        cfg["prefix"] = args.out_prefix
    if args.title:
        cfg["title"] = args.title

    if not os.path.exists(cfg["csv"]):
        raise FileNotFoundError(f"{cfg['csv']} not found")

    df = pd.read_csv(cfg["csv"])
    if cfg["column"] not in df.columns:
        raise KeyError(f"column '{cfg['column']}' not found in {cfg['csv']}")

    seqs = df[cfg["column"]].dropna().astype(str).tolist()
    seqs = [s for s in seqs if len(s) > 0]
    lengths = [len(s) for s in seqs]
    most_common_length = max(set(lengths), key=lengths.count)
    seqs = [s for s in seqs if len(s) == most_common_length]
    if cfg["top"]:
        seqs = seqs[:cfg["top"]]
    if not seqs:
        raise ValueError("no valid sequences found")
    print(f"sequences: {len(seqs)}, length: {most_common_length}")

    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    matrix = pd.DataFrame(0, index=list(amino_acids), columns=range(most_common_length))
    for seq in seqs:
        for i, aa in enumerate(seq):
            if aa in amino_acids:
                matrix.loc[aa, i] += 1
    matrix = matrix.div(matrix.sum(axis=0), axis=1).T

    logo = logomaker.Logo(matrix, shade_below=.5, fade_below=.5,
                          font_name="Arial", color_scheme="chemistry")
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.style_xticks(anchor=0, spacing=1, rotation=0)
    logo.ax.set_xticklabels([str(i + 1) for i in range(most_common_length)])
    logo.ax.set_ylabel("Probability")
    logo.ax.set_xlabel("Position")
    logo.ax.set_ylim(-0.01, 1.01)
    logo.ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
    if cfg["title"]:
        logo.ax.set_title(cfg["title"])

    handles = [Patch(facecolor=c, label=t) for c, t in CHEMISTRY_LEGEND]
    logo.ax.legend(handles=handles, loc="upper center",
                   bbox_to_anchor=(0.5, -0.10), ncol=3, frameon=False,
                   fontsize=8, handlelength=1.2, columnspacing=1.2)

    os.makedirs(cfg["out_dir"], exist_ok=True)
    plt.tight_layout()
    plt.gcf().set_size_inches(8, 6)
    pdf_path = os.path.join(cfg["out_dir"], cfg["prefix"] + ".pdf")
    png_path = os.path.join(cfg["out_dir"], cfg["prefix"] + ".png")
    plt.savefig(pdf_path, format="pdf", dpi=300)
    plt.savefig(png_path, format="png", dpi=300)
    plt.close()
    print(f"saved: {pdf_path}\n       {png_path}")


if __name__ == "__main__":
    main()
