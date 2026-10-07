"""
Flu mutagenesis screen: interactive peptide network (pyecharts force-directed)
-> flu_interactive_network.html

Nodes: top-200 barcodes, stop-codon-free -> r1_r2_fold_change <= 500 -> P2=I & P9=L
       -> 9-mer -> first 200, deduplicated = 125 peptides (matches the published figure)
Clustering: Hamming distance + hierarchical clustering (average linkage, 12 clusters)
Labels: isolated nodes always shown; connected nodes sampled at
        rand < DISPLAY_RATIO*1.5/(1+degree)
        (the published figure shows 46/125 labels; re-running samples a different
         label set - to reproduce the published figure exactly, use the committed
         results/flu_network/nodes.csv + plot_network_vector.py)
Next step: code/shared/extract_network_data.py extracts nodes.csv / links.csv from this HTML

Usage:
- Run: `python pyechart_flu.py` (from `code/flu_mutation_screen/`).
- Output: `results/flu_network/flu_interactive_network.html`.

Requirements:
- Python 3.x, pandas, numpy, scikit-learn, pyecharts
"""
import os
import itertools
from typing import List

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from pyecharts.charts import Graph
from pyecharts.options import GraphNode, GraphLink, LabelOpts

# Input file path
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
INPUT_FILE_PATH = os.path.join(REPO, "data", "top10000_combine_9mer_2308-Flu-mut-SCD-round123.csv")
# Output directory
NETWORK_DIR = os.path.join(REPO, "results", "flu_network")
os.makedirs(NETWORK_DIR, exist_ok=True)

# Output file name
OUTPUT_HTML_FILE = "flu_interactive_network.html"

# Parameters
DISTANCE_ALPHA_LUT = {
    1: 1.0,  # closest distance, highest opacity
    2: 0.8,  # intermediate distance, higher opacity
    3: 0.6,  # far distance, medium opacity
}
REPULSION = 10       # repulsion between nodes
NODE_SYMBOL_SIZE = 5  # node size
EDGE_LINE_WIDTH = 1.0  # edge width
DISPLAY_RATIO = 3    # label display ratio (rand < DISPLAY_RATIO*1.5/(1+degree))
SEED = 42            # label random-sampling seed (not fixed in the published HTML,
                     # so re-runs sample a different label set)

# Published cluster palette (by sklearn cluster label 0-11, identical to the
# published HTML / nodes.csv):
# label 2 = main cluster, 94 peptides (purple); 0 = 18 peptides (red);
# 1 = 4 peptides (pink); 3-11 = singleton clusters
# (singleton colors are for interactive browsing only; patch_root_nodes.py
#  forces them to gray after extraction)
cluster_colors = [
    "#F44336",  # label 0 -> 18 peptides
    "#E91E63",  # label 1 -> 4 peptides
    "#673AB7",  # label 2 -> 94 peptides (main cluster)
    "#9C27B0",  # label 3-11 -> singleton clusters
    "#3F51B5", "#2196F3", "#03A9F4", "#00BCD4",
    "#009688", "#4CAF50", "#8BC34A", "#CDDC39",
]

# Read the CSV file
df = pd.read_csv(INPUT_FILE_PATH, sep=",")

# Check that the 'aa' column exists
if 'aa' not in df.columns:
    raise KeyError("The column 'aa' is not found in the CSV file. Please check the column names.")

# Extract peptide sequences
peptide_sequences = df.iloc[:200]  # top 200 rows
# Drop peptides containing a stop codon (*)
peptide_sequences = peptide_sequences[~peptide_sequences['aa'].str.contains(r'\*', na=False)]

# Check that the required columns exist
required_columns = ['count_pct_round1', 'count_pct_round2', 'count_pct_round3']
for col in required_columns:
    if col not in peptide_sequences.columns:
        raise KeyError(f"The column '{col}' is missing in the input file.")

# Compute r1-r2 fold change for filtering
peptide_sequences = peptide_sequences.copy()
peptide_sequences['r1_r2_fold_change'] = peptide_sequences['count_pct_round2'] / peptide_sequences['count_pct_round1']

# Drop rows with r1_r2_fold_change > 500
peptide_sequences = peptide_sequences[peptide_sequences['r1_r2_fold_change'] <= 500]

# Extract the 'aa' column as a list
peptide_sequences = peptide_sequences['aa'].dropna().astype(str).tolist()

# Keep only sequences with I at position 2 and L at position 9
valid_sequences = [seq for seq in peptide_sequences if len(seq) >= 9 and seq[1] == 'I' and seq[8] == 'L']

# Keep only sequences of the most common length
sequence_lengths = [len(seq) for seq in valid_sequences]
most_common_length = max(set(sequence_lengths), key=sequence_lengths.count)
valid_sequences = [seq for seq in valid_sequences if len(seq) == most_common_length]

# Take the first 200 sequences
valid_sequences = valid_sequences[:200]
valid_sequences = sorted(set(valid_sequences))
print(f"{len(valid_sequences)} peptides selected")

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

# Hierarchical clustering
n_clusters = 12  # number of clusters
# sklearn >= 1.2: the 'affinity' parameter was renamed to 'metric'
clustering = AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage='average')
clusters = clustering.fit_predict(hamming_matrix)
print("cluster sizes:", dict(zip(*np.unique(clusters, return_counts=True))))

# Edges (as tuples)
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

seq_index = {s: i for i, s in enumerate(valid_sequences)}

# GraphLink objects
links = [
    GraphLink(
        source=source,
        target=target,
        linestyle_opts={
            "color": "#000000",
            "width": EDGE_LINE_WIDTH,
            # edge opacity by Hamming distance: the closer, the more opaque
            "opacity": DISTANCE_ALPHA_LUT.get(hamming_matrix[seq_index[source], seq_index[target]], 0.5),
        }
    )
    for source, target in raw_links
]

# Nodes (label decision computed separately for reporting)
np.random.seed(SEED)
label_show = [
    node_connections[s] == 0
    or np.random.rand() < (DISPLAY_RATIO * 1.5) / (1 + node_connections[s])
    for s in valid_sequences
]
print(f"edges: {len(raw_links)}, labels shown: {sum(label_show)} / {len(valid_sequences)}")

nodes = [
    GraphNode(
        name=valid_sequences[i],
        symbol="circle",
        symbol_size=NODE_SYMBOL_SIZE,
        itemstyle_opts={"color": cluster_colors[clusters[i]]},
        label_opts={
            "show": label_show[i],
            "position": "top" if node_connections[valid_sequences[i]] == 0 else "right",
            "font_size": 12,
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
output_html_path = os.path.join(NETWORK_DIR, OUTPUT_HTML_FILE)
graph.render(output_html_path)

print(f"Graph saved to {output_html_path}")
