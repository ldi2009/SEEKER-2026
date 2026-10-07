# -*- coding: utf-8 -*-
"""Extract peptide network data from flu_interactive_network.html (pyecharts force
graph) -> nodes.csv / links.csv"""
import re
import csv
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(REPO, "results", "flu_network", "flu_interactive_network.html")
OUT_DIR = os.path.join(REPO, "results", "flu_network")

with open(SRC, "r", encoding="utf-8") as f:
    html = f.read()

data_start = html.index('"data":')
links_start = html.index('"links":')
data_sec = html[data_start:links_start]
links_sec = html[links_start:]

# ---- Nodes ----
node_blocks = re.findall(
    r'\{\s*"name":\s*"([^"]+)",[\s\S]*?"symbolSize":\s*([\d.]+),[\s\S]*?"color":\s*"(#[0-9A-Fa-f]{6})"[\s\S]*?"label":\s*\{([^}]*)\}',
    data_sec,
)
nodes = []
for name, size, color, label_body in node_blocks:
    label_show = '"show":true' in re.sub(r"\s", "", label_body)
    nodes.append({"name": name, "symbol_size": float(size), "color": color, "label_show": label_show})
print(f"nodes parsed: {len(nodes)}")

# ---- Edges ----
edge_blocks = re.findall(
    r'\{\s*"source":\s*"([^"]+)",\s*"target":\s*"([^"]+)",\s*"lineStyle":\s*\{[\s\S]*?"width":\s*([\d.]+),\s*"opacity":\s*([\d.]+)',
    links_sec,
)
OPACITY_TO_DIST = {1.0: 1, 0.8: 2, 0.6: 3}
edges, mismatch = [], 0
for source, target, width, opacity in edge_blocks:
    dist = OPACITY_TO_DIST.get(float(opacity))
    actual = sum(a != b for a, b in zip(source, target))
    if dist is None or actual != dist:
        mismatch += 1
    edges.append({"source": source, "target": target, "hamming": actual, "opacity": float(opacity)})
print(f"edges parsed: {len(edges)}, opacity<->hamming mismatch: {mismatch}")

# ---- Statistics ----
degree = {n["name"]: 0 for n in nodes}
deg1 = {n["name"]: 0 for n in nodes}
for e in edges:
    degree[e["source"]] += 1
    degree[e["target"]] += 1
    if e["hamming"] == 1:
        deg1[e["source"]] += 1
        deg1[e["target"]] += 1

color_counts = Counter(n["color"] for n in nodes)
print("cluster sizes by color:", dict(color_counts))
labeled = [n["name"] for n in nodes if n["label_show"]]
print(f"labeled: {len(labeled)}, labeled degree range: {min(degree[x] for x in labeled)}-{max(degree[x] for x in labeled)}")
unlabeled = [n["name"] for n in nodes if not n["label_show"]]
print(f"unlabeled degree range: {min(degree[x] for x in unlabeled)}-{max(degree[x] for x in unlabeled)}")
print("special nodes:", [(n["name"], n["symbol_size"]) for n in nodes if n["symbol_size"] != 5])

# color -> cluster id (by cluster size, descending)
order = [c for c, _ in color_counts.most_common()]
color_to_cluster = {c: i for i, c in enumerate(order)}

# ---- Export ----
with open(os.path.join(OUT_DIR, "nodes.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["name", "cluster", "color", "label_show", "degree", "degree_dist1"])
    for n in nodes:
        w.writerow([n["name"], color_to_cluster[n["color"]], n["color"], int(n["label_show"]),
                    degree[n["name"]], deg1[n["name"]]])

with open(os.path.join(OUT_DIR, "links.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["source", "target", "hamming"])
    for e in edges:
        w.writerow([e["source"], e["target"], e["hamming"]])

print("Saved nodes.csv / links.csv ->", OUT_DIR)
