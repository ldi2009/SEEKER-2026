import os
import math
import itertools
from typing import List, Dict

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
from pyecharts.charts import Graph
from pyecharts.options import GraphNode, GraphLink, LabelOpts

# Input file path
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
INPUT_FILE_PATH = os.path.join(REPO, "data", "top10000_combine_9mer_2305_DMF5_SCD_R123.csv")
# Output directory
OUTPUT_DIR = os.path.join(REPO, "results", "dmf5_network")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output file names
OUTPUT_HTML_FILE = "dmf5_interactive_network.html"
OUTPUT_DENDROGRAM_FILE = "DMF5_dendrogram.png"  # clustering dendrogram output

# Parameters
DISTANCE_ALPHA_LUT = {
    1: 1.0,
    2: 0.8,
    3: 0.6,
    4: 0.1
}
REPULSION = 10  # repulsion between nodes
NODE_SYMBOL_SIZE = 5  # node size
EDGE_LINE_WIDTH = 1.0  # edge width
DISPLAY_RATIO = 1  # label display ratio

# Fixed color palette
COLOR_PALETTE = [
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
    "#AEC7E8", "#FFBB78"
]

# Read the CSV file
df = pd.read_csv(INPUT_FILE_PATH, header=0)

# Check that the 'aa' column exists
if 'aa' not in df.columns:
    raise KeyError("The 'aa' column is not found in the CSV file. Please check the column names.")

# Extract peptide sequences
peptide_sequences = df['aa'].dropna().astype(str).tolist()

# Keep only sequences of the most common length
sequence_lengths = [len(seq) for seq in peptide_sequences]
most_common_length = max(set(sequence_lengths), key=sequence_lengths.count)
valid_sequences = [seq for seq in peptide_sequences if len(seq) == most_common_length]

# Take the first 300 sequences
valid_sequences = valid_sequences[:300]
valid_sequences = sorted(set(valid_sequences))  # stable ordering
print(valid_sequences)

# Export the selected peptides as a CSV file
selected_peptides_df = pd.DataFrame(valid_sequences, columns=['aa'])
selected_peptides_csv_path = os.path.join(OUTPUT_DIR, "selected_peptides.csv")
selected_peptides_df.to_csv(selected_peptides_csv_path, index=False)
print(f"Selected peptides saved to {selected_peptides_csv_path}")

# Check that enough sequences remain
if len(valid_sequences) == 0:
    raise ValueError("No valid sequences found with the most common length.")

# Hamming distance matrix
def hamming_distance_matrix(sequences: List[str]) -> np.ndarray:
    n = len(sequences)
    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            distance = sum(c1 != c2 for c1, c2 in zip(sequences[i], sequences[j]))
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance
    return distance_matrix

# Compute the Hamming distance matrix
hamming_matrix = hamming_distance_matrix(valid_sequences)
print(hamming_matrix)

# Hierarchical clustering
n_clusters = 12  # number of clusters
# sklearn >= 1.2: the 'affinity' parameter was renamed to 'metric'
clustering = AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage='average')
clusters = clustering.fit_predict(hamming_matrix)

# Assign colors from the fixed palette
cluster_colors = [COLOR_PALETTE[i % len(COLOR_PALETTE)] for i in range(n_clusters)]
print(cluster_colors)

from scipy.spatial.distance import squareform

# Draw the dendrogram
linkage_matrix = linkage(squareform(hamming_matrix), method='average')

# Color map for leaves
leaf_colors = {valid_sequences[i]: cluster_colors[clusters[i]] for i in range(len(valid_sequences))}

# Custom dendrogram color function
def dendrogram_color_func(link_id):
    # Leaf node: return its cluster color
    if link_id < len(valid_sequences):
        return leaf_colors.get(valid_sequences[link_id], "#000000")
    # Branch node: return the color if both children match, else black
    else:
        left_child = int(linkage_matrix[link_id - len(valid_sequences), 0])
        right_child = int(linkage_matrix[link_id - len(valid_sequences), 1])
        left_color = dendrogram_color_func(left_child)
        right_color = dendrogram_color_func(right_child)
        return left_color if left_color == right_color else "#000000"

plt.figure(figsize=(10, 7))
dendrogram(
    linkage_matrix,
    labels=valid_sequences,
    leaf_rotation=90,
    leaf_font_size=8,
    link_color_func=lambda x: dendrogram_color_func(x)  # custom color function
)
plt.title("Dendrogram of Peptide Clustering")
plt.xlabel("Peptides")
plt.ylabel("Hamming Distance")
plt.tight_layout()
dendrogram_path = os.path.join(OUTPUT_DIR, OUTPUT_DENDROGRAM_FILE)
plt.savefig(dendrogram_path)
plt.close()
print(f"Dendrogram saved to {dendrogram_path}")

# Edges
raw_links = [
    (valid_sequences[i], valid_sequences[j])
    for i, j in itertools.combinations(range(len(valid_sequences)), 2)
    if hamming_matrix[i, j] in DISTANCE_ALPHA_LUT
]

# Node connection counts
node_connections = {node: 0 for node in valid_sequences}
for source, target in raw_links:
    node_connections[source] += 1
    node_connections[target] += 1

# GraphLink objects
links = [
    GraphLink(
        source=source,
        target=target,
        linestyle_opts={
            "color": "#000000",
            "width": EDGE_LINE_WIDTH if hamming_matrix[valid_sequences.index(source), valid_sequences.index(target)] in DISTANCE_ALPHA_LUT else 0.5,
            "opacity": DISTANCE_ALPHA_LUT.get(hamming_matrix[valid_sequences.index(source), valid_sequences.index(target)], 0.3)
        }
    )
    for source, target in raw_links
]

# Nodes
nodes = [
    GraphNode(
        name=valid_sequences[i],
        symbol="circle",
        symbol_size=NODE_SYMBOL_SIZE,
        itemstyle_opts={"color": cluster_colors[clusters[i]]},
        label_opts={
            # boost label display probability
            "show": node_connections[valid_sequences[i]] == 0 or np.random.rand() < (DISPLAY_RATIO * 1.5) / (1 + node_connections[valid_sequences[i]]),
            "position": "top" if node_connections[valid_sequences[i]] == 0 else "right",
            "font_size": 12
        }
    )
    for i in range(len(valid_sequences))
]

# Render with pyecharts
graph = Graph(init_opts={"width": "1000px", "height": "1000px"})
graph.add(
    series_name="Peptide Network",
    nodes=nodes,
    links=links,
    repulsion=REPULSION,
    label_opts=LabelOpts(is_show=True, position="top", font_size=8)
)

# Save the chart as HTML
output_html_path = os.path.join(OUTPUT_DIR, OUTPUT_HTML_FILE)
graph.render(output_html_path)

print(f"Graph saved to {output_html_path}")