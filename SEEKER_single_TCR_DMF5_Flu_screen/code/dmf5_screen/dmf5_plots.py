"""
DMF5 screen amino-acid usage heatmaps: plain and single-peptide-line versions,
plus the frequency matrix CSV (sequence logos are produced by code/shared/plot_seqlogo.py).

Usage:
- Run: `python dmf5_plots.py` (from `code/dmf5_screen/`).
- Outputs saved to `results/dmf5_screen/`.

Requirements:
- Python 3.x
- pandas, numpy, matplotlib
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Input file path
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
INPUT_FILE_PATH = os.path.join(REPO, "data", "top10000_combine_9mer_2305_DMF5_SCD_R123.csv")
# Output directory
OUTPUT_DIR = os.path.join(REPO, "results", "dmf5_screen")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Check that the input file exists
if not os.path.exists(INPUT_FILE_PATH):
    raise FileNotFoundError(f"File {INPUT_FILE_PATH} not found. Please check the path.")

# Amino acid alphabet
AMINO_ACIDS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def read_input_file():
    """
    Read the input CSV file.
    :return: DataFrame with the data
    """
    df = pd.read_csv(INPUT_FILE_PATH, sep=",", header=0)
    print("Columns in the file:", df.columns)
    df.columns = df.columns.str.strip()
    if "aa" not in df.columns:
        raise KeyError("The column 'aa' is missing in the file.")
    return df
def draw_peptide_heatmap(peptides, with_lines, output_file, matrix_file=None):
    """
    Generate and save the peptide heatmap.
    :param peptides: list of peptide sequences
    :param with_lines: whether to trace single-peptide lines
    :param output_file: output file name
    :param matrix_file: frequency matrix file name (optional)
    """
    matrix = np.zeros((len(AMINO_ACIDS), len(peptides[0])))
    for i, aa in enumerate(AMINO_ACIDS):
        for j in range(len(peptides[0])):
            freq = sum([p[j] == aa for p in peptides]) / len(peptides)
            matrix[i, j] = freq
    if matrix_file:
        df_matrix = pd.DataFrame(matrix, index=AMINO_ACIDS, columns=[f"Pos_{i+1}" for i in range(len(peptides[0]))])
        matrix_path = os.path.join(OUTPUT_DIR, matrix_file)
        df_matrix.to_csv(matrix_path, index=True)
        print(f"Frequency matrix saved to {matrix_path}")
    fig, ax = plt.subplots()
    cmap = plt.get_cmap('Blues')
    im = ax.imshow(matrix, cmap=cmap)
    cbar = ax.figure.colorbar(im, ax=ax, ticks=np.arange(0, 1.01, 0.1))
    cbar.ax.set_ylabel("Frequency", rotation=-90, va="bottom")
    ax.set_xticks(np.arange(len(peptides[0])))
    ax.set_yticks(np.arange(len(AMINO_ACIDS)))
    ax.set_xticklabels(np.arange(1, len(peptides[0]) + 1))
    ax.set_yticklabels(AMINO_ACIDS)
    if with_lines:
        for p in peptides:
            xs = []
            ys = []
            for i, aa in enumerate(p):
                if aa in AMINO_ACIDS:
                    x = i
                    y = AMINO_ACIDS.index(aa)
                    xs.append(x)
                    ys.append(y)
            plt.plot(xs, ys, color='black', alpha=0.01)
    plt.xlabel('Peptide position')
    plt.ylabel('Amino acid')
    plt.title('Peptide sequence distribution')
    plt.tight_layout()
    pdf_path = os.path.join(OUTPUT_DIR, output_file)
    png_path = os.path.join(OUTPUT_DIR, output_file.replace(".pdf", ".png"))
    plt.savefig(png_path, format="png")
    plt.savefig(pdf_path, format="pdf")
    plt.close()
    print(f"Heatmap saved to {pdf_path} and {png_path}")

def main():
    """
    Main entry point: draw the heatmaps.
    """
    df = read_input_file()
    # Take the top 300 peptide sequences
    peptide_sequences = df.iloc[:300, 1].str.upper().tolist()
    # Draw the heatmaps
    draw_peptide_heatmap(peptide_sequences, with_lines=False, output_file="DMF5_9mer_heatmap.pdf", matrix_file="frequency_DMF5_9mer_heatmap.csv")
    draw_peptide_heatmap(peptide_sequences, with_lines=True, output_file="lined_DMF5_9mer_heatmap.pdf")

if __name__ == "__main__":
    main()