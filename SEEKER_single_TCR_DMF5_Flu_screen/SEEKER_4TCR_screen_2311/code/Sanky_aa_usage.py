"""Archived Flu-screen script (kept as-is, translated).

Original single-TCR workflow for the Flu screen: plain + lined heatmaps,
frequency matrix and an amino-acid transition Sankey diagram.
Inputs come from ../data/, outputs go to ../results/flu_archive/.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go  # used for the Sankey diagram

# input file (Flu top-10,000 table)
INPUT_FILE_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "data",
    "top10000_combine_9mer_2311-4TCR-Flu-SCD-round123.csv")
# output directory
OUTPUT_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "flu_archive")
# output file names
OUTPUT_HEATMAP_FILE = "Flu_4TCR_heatmap.pdf"
OUTPUT_LINED_HEATMAP_FILE = "lined_Flu_4TCR_heatmap.pdf"
OUTPUT_MATRIX_FILE = "frequency_Flu_4TCR_heatmap.csv"
OUTPUT_SANKY_FILE = "sankey_Flu_4TCR.html"

if not os.path.exists(INPUT_FILE_PATH):
    raise FileNotFoundError(f"Input file not found: {INPUT_FILE_PATH}")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# read the table (top rows are ranked by count_round3)
df = pd.read_csv(INPUT_FILE_PATH, header=0)
# take the top 50 rows, drop starred peptides
peptide_sequences = df.iloc[:50]
peptide_sequences = peptide_sequences[~peptide_sequences['aa'].str.contains(r'\*', na=False)]

# keep r3_r2_ratio >= 40
if 'r3_r2_ratio' not in peptide_sequences.columns:
    raise KeyError("The column 'r3_r2_ratio' is missing in the input file.")
peptide_sequences = peptide_sequences[peptide_sequences['r3_r2_ratio'] >= 40]

# extract the 'aa' column as a list
peptide_sequences = peptide_sequences['aa'].dropna().astype(str).tolist()

# amino-acid alphabet
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def draw_peptide_heatmap(peptides, amino_acids, with_lines=False, output_file="heatmap.pdf", matrix_file=None):
    # build the frequency matrix
    matrix = np.zeros((len(amino_acids), len(peptides[0])))
    for i, aa in enumerate(amino_acids):
        for j in range(len(peptides[0])):
            freq = sum([p[j] == aa for p in peptides]) / len(peptides)
            matrix[i, j] = freq

    # save the matrix as CSV if requested
    if matrix_file:
        df_matrix = pd.DataFrame(matrix, index=amino_acids, columns=[f"Pos_{i+1}" for i in range(len(peptides[0]))])
        df_matrix.to_csv(matrix_file, index=True)

    # heatmap
    fig, ax = plt.subplots()
    cmap = plt.get_cmap('Blues')
    im = ax.imshow(matrix, cmap=cmap)

    # color bar
    cbar = ax.figure.colorbar(im, ax=ax, ticks=np.arange(0, 1.01, 0.1))
    cbar.ax.set_ylabel("Frequency", rotation=-90, va="bottom")

    # tick labels
    ax.set_xticks(np.arange(len(peptides[0])) + 0)
    ax.set_yticks(np.arange(len(amino_acids)) + 0)
    ax.set_xticklabels(np.arange(1, len(peptides[0]) + 1))
    ax.set_yticklabels(amino_acids)

    # optional per-peptide polylines through the cell centers
    if with_lines:
        for p in peptides:
            xs = []
            ys = []
            for i, aa in enumerate(p):
                if aa in amino_acids:  # skip non-standard residues
                    xs.append(i)
                    ys.append(amino_acids.index(aa))
            plt.plot(xs, ys, color='black', alpha=0.01)

    plt.xlabel('Peptide position')
    plt.ylabel('Amino acid')
    plt.title('Peptide sequence distribution')
    plt.savefig(output_file, format='pdf')
    plt.close()


def draw_sankey(peptides, amino_acids, output_file="sankey.html"):
    # build the Sankey transition data
    links = {}
    for p in peptides:
        for i in range(len(p) - 1):
            source = f"{p[i]}_{i+1}"
            target = f"{p[i+1]}_{i+2}"
            if (source, target) not in links:
                links[(source, target)] = 0
            links[(source, target)] += 1

    # nodes and links
    nodes = sorted(set([key[0] for key in links.keys()] + [key[1] for key in links.keys()]))
    node_indices = {node: i for i, node in enumerate(nodes)}

    sankey_links = {
        "source": [node_indices[source] for source, target in links.keys()],
        "target": [node_indices[target] for source, target in links.keys()],
        "value": list(links.values())
    }

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes
        ),
        link=dict(
            source=sankey_links["source"],
            target=sankey_links["target"],
            value=sankey_links["value"]
        )
    )])

    fig.update_layout(title_text="Sankey Diagram of Amino Acid Transitions", font_size=10)
    fig.write_html(output_file)
    print(f"Sankey diagram saved to {output_file}")


# plain heatmap + frequency matrix
draw_peptide_heatmap(
    peptide_sequences,
    amino_acids,
    with_lines=False,
    output_file=os.path.join(OUTPUT_DIR, OUTPUT_HEATMAP_FILE),
    matrix_file=os.path.join(OUTPUT_DIR, OUTPUT_MATRIX_FILE)
)

# lined heatmap
draw_peptide_heatmap(
    peptide_sequences,
    amino_acids,
    with_lines=True,
    output_file=os.path.join(OUTPUT_DIR, OUTPUT_LINED_HEATMAP_FILE)
)

# Sankey diagram
draw_sankey(
    peptide_sequences,
    amino_acids,
    output_file=os.path.join(OUTPUT_DIR, OUTPUT_SANKY_FILE)
)
