import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Input file path
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
INPUT_FILE_PATH = os.path.join(REPO, "data", "top10000_combine_9mer_2308-Flu-mut-SCD-round123.csv")
# Output directory
OUTPUT_DIR = os.path.join(REPO, "results", "flu_mutation_screen")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Check that the input file exists
if not os.path.exists(INPUT_FILE_PATH):
    raise FileNotFoundError(f"File {INPUT_FILE_PATH} not found. Please check the path.")

df = pd.read_csv(INPUT_FILE_PATH, header=0)
# Extract peptide sequences
peptide_sequences = df.iloc[:200]  # top 200 rows
# Drop peptides containing a stop codon (*)
peptide_sequences = peptide_sequences[~peptide_sequences['aa'].str.contains(r'\*', na=False)]

# Optional enrichment filter on r3_r2_ratio (disabled, matching the published figure)
if 'r3_r2_ratio' not in peptide_sequences.columns:
    raise KeyError("The column 'r3_r2_ratio' is missing in the input file.")
#peptide_sequences = peptide_sequences[peptide_sequences['r3_r2_ratio'] >= 40]

# Save the selected rows as a new CSV file
selected_file_path = os.path.join(OUTPUT_DIR, "selected_Flu_mutation_screen_peptide.csv")
peptide_sequences.to_csv(selected_file_path, index=False)
print(f"Selected peptide rows saved to {selected_file_path}")

# Extract the 'aa' column as a list
peptide_sequences = peptide_sequences['aa'].dropna().astype(str).tolist()
# Amino acid alphabet
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


# Heatmap plotting function
def draw_peptide_heatmap(peptides, amino_acids, with_lines=False, output_file="heatmap.pdf", matrix_file=None):
    # Build the frequency matrix
    matrix = np.zeros((len(amino_acids), len(peptides[0])))
    for i, aa in enumerate(amino_acids):
        for j in range(len(peptides[0])):
            freq = sum([p[j] == aa for p in peptides]) / len(peptides)
            matrix[i, j] = freq

    # Save the matrix as CSV if a path is given
    if matrix_file:
        df_matrix = pd.DataFrame(matrix, index=amino_acids, columns=[f"Pos_{i+1}" for i in range(len(peptides[0]))])
        df_matrix.to_csv(matrix_file, index=True)

    # Draw the heatmap
    fig, ax = plt.subplots()
    cmap = plt.get_cmap('Blues')
    im = ax.imshow(matrix, cmap=cmap)

    # Colorbar
    cbar = ax.figure.colorbar(im, ax=ax, ticks=np.arange(0, 1.01, 0.1))
    cbar.ax.set_ylabel("Frequency", rotation=-90, va="bottom")

    # Tick labels
    ax.set_xticks(np.arange(len(peptides[0])) + 0)
    ax.set_yticks(np.arange(len(amino_acids)) + 0)
    ax.set_xticklabels(np.arange(1, len(peptides[0]) + 1))
    ax.set_yticklabels(amino_acids)

    # Optionally trace each single peptide through the grid
    if with_lines:
        for p in peptides:
            xs = []
            ys = []
            for i, aa in enumerate(p):
                if aa in amino_acids:  # skip non-standard residues
                    x = i
                    y = amino_acids.index(aa)
                    xs.append(x)
                    ys.append(y)
            plt.plot(xs, ys, color='black', alpha=0.01)

    # Labels, title and export
    plt.xlabel('Peptide position')
    plt.ylabel('Amino acid')
    plt.title('Peptide sequence distribution')
    plt.savefig(output_file, format='pdf')  # save as PDF
    plt.close()


# Heatmap without lines + frequency matrix
draw_peptide_heatmap(
    peptide_sequences,
    amino_acids,
    with_lines=False,
    output_file=os.path.join(OUTPUT_DIR, "Flu_mutation_heatmap.pdf"),
    matrix_file=os.path.join(OUTPUT_DIR, "frequency_Flu_mutation_heatmap.csv")
)

# Heatmap with single-peptide lines
draw_peptide_heatmap(
    peptide_sequences,
    amino_acids,
    with_lines=True,
    output_file=os.path.join(OUTPUT_DIR, "lined_Flu_mutation_heatmap.pdf")
)