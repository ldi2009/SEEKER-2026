# -*- coding: utf-8 -*-
"""
Peptide Hamming-distance network - publication vector figure (PDF / SVG / EPS / PNG)
Data: results/flu_network/nodes.csv / links.csv (extracted from flu_interactive_network.html
by extract_network_data.py)
Outputs 3 edge versions: all edges (1-3) / distance<=2 / distance=1, sharing one layout
"""
import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
from adjustText import adjust_text

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
NET_DIR = os.path.join(REPO, "results", "flu_network")
WT = "GILGFVFTL"
SEED = 42

# ---------------- Tunable parameters ----------------
FIGSIZE_IN = (7.09, 4.92)          # 180 x 125 mm wide two-column figure
NODE_SIZE = 26                     # node area (pt^2)
WT_NODE_SIZE = 60                  # WT epitope node area
LABEL_FONTSIZE = 5.5
WT_LABEL_FONTSIZE = 7
LAYOUT_K = 0.32                    # spring layout node spacing
LAYOUT_ITER = 1000
RING_GAP = 0.008                   # gap between the isolate ring and the main component
                                   # (direction-aware, hugs the cluster outline)

# Light-color version: node/edge colors blended toward white, outputs get a _light
# suffix (dark version kept), and WT is annotated with an arrow
LIGHT = True
LIGHT_KEEP = 0.55                  # fraction of cluster color kept (rest blended to white)
LIGHT_GRAY_KEEP = 0.75             # fraction of singleton gray kept
WT_ARROW_COLOR = "#B26A00"         # WT arrow and text (amber)
ANN_DX, ANN_DY = -0.018, -0.028    # WT arrow text offset relative to the WT node

EDGE_STYLE = {                     # per Hamming distance: (color, width, alpha)
    1: ("#1a1a1a", 0.9, 0.90),
    2: ("#555555", 0.5, 0.45),
    3: ("#9e9e9e", 0.28, 0.18),
}
LIGHT_EDGE_STYLE = {               # light-version edge colors: one step softer
    1: ("#333333", 0.9, 0.90),
    2: ("#777777", 0.5, 0.45),
    3: ("#b3b3b3", 0.28, 0.18),
}
VERSIONS = {                       # filename suffix -> distance bins kept
    "all_edges": [1, 2, 3],
    "dist_le2": [1, 2],
    "dist_eq1": [1],
}

plt.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,             # embedded TrueType, journals require editable text
    "ps.fonttype": 42,
    "svg.fonttype": "none",         # keep text in SVG (not paths), editable in AI/Inkscape
})

# ---------------- Data ----------------
nodes_df = pd.read_csv(os.path.join(NET_DIR, "nodes.csv"))
links_df = pd.read_csv(os.path.join(NET_DIR, "links.csv"))

G = nx.Graph()
for _, r in nodes_df.iterrows():
    G.add_node(r["name"], cluster=int(r["cluster"]), color=r["color"],
               label_show=bool(r["label_show"]), degree=int(r["degree"]))
for _, r in links_df.iterrows():
    d = int(r["hamming"])
    G.add_edge(r["source"], r["target"], hamming=d, spring_w=1.0 / d ** 2)

# ---------------- Layout (fixed seed, shared by all three versions) ----------------
# Strategy: d<=2 edges dominate the layout (w=1/d^3); nodes connected only by d=3
# edges attach weakly (w=1/27); fully isolated nodes (degree=0) go to an outer ring
# so they don't stretch the canvas
H12 = nx.Graph()
H12.add_nodes_from(G.nodes)
for u, v, d in G.edges(data="hamming"):
    if d <= 2:
        H12.add_edge(u, v)
main_comp = max(nx.connected_components(H12), key=len)
isolates = [n for n in G.nodes if n not in main_comp and G.degree(n) == 0]

L = nx.Graph()
L.add_nodes_from(G.nodes)
for u, v, d in G.edges(data="hamming"):
    if d <= 2:
        L.add_edge(u, v, spring_w=1.0 / d ** 3)
    elif u not in main_comp or v not in main_comp:
        # nodes linked only by d=3 edges (independent sequence families) attach
        # weakly; they don't pull on the main component's interior
        L.add_edge(u, v, spring_w=1.0 / 27)

pos = nx.spring_layout(L, k=0.55, iterations=2000, seed=SEED, weight="spring_w")
main_xy = np.array([pos[n] for n in main_comp])
_c = main_xy.mean(axis=0)
_r_max = np.linalg.norm(main_xy - _c, axis=1).max()
main_r = np.linalg.norm(main_xy - _c, axis=1)
main_ang = np.arctan2(main_xy[:, 1] - _c[1], main_xy[:, 0] - _c[0])

# Anisotropic gap: the wide canvas scales x and y differently, so the y gap is
# reduced by the render ratio to keep the ring-cluster visual gap isotropic.
# Note: only nodes with edges are counted - degree-0 isolates fly to +-1 in the
# spring layout and would badly skew the aspect-ratio estimate
_active_xy = np.array([pos[n] for n in G.nodes if G.degree(n) > 0])
_ar = (FIGSIZE_IN[0] / FIGSIZE_IN[1]) * (np.ptp(_active_xy[:, 1]) / np.ptp(_active_xy[:, 0]))
GAP_X, GAP_Y = RING_GAP, RING_GAP * max(_ar, 0.25)

def ring_pos(ang):
    # direction-aware: farthest main-component node radius within +-25 deg of this
    # direction, plus the gap - the ring hugs the cluster outline instead of a circle
    d = np.abs((main_ang - ang + np.pi) % (2 * np.pi) - np.pi)
    sel = d <= np.deg2rad(25)
    r_dir = main_r[sel].max() if sel.any() else 0.0
    r_dir = max(r_dir, _r_max * 0.5)
    ux, uy = np.cos(ang), np.sin(ang)
    return (_c[0] + (r_dir + GAP_X) * ux, _c[1] + (r_dir + GAP_Y) * uy)

for i, n in enumerate(sorted(isolates)):
    pos[n] = ring_pos(np.pi / 2 + 2 * np.pi * i / len(isolates))

# Singleton-cluster nodes that still have weak d=3 edges get flung far away by the
# spring layout (RIKEARFTL/EINVWVRTL); pull them back to free angles on the isolate
# ring so the main structure fills the canvas
cl_size = nodes_df["cluster"].value_counts()
weak_singletons = [n for n in G.nodes
                   if n not in main_comp and G.degree(n) > 0
                   and int(cl_size.loc[int(G.nodes[n]["cluster"])]) == 1]
if weak_singletons and isolates:
    N = len(isolates)
    mids = [np.pi / 2 + 2 * np.pi * (i + 0.5) / N for i in range(N)]  # gaps between isolates
    used = [False] * N

    def cur_ang(n):
        return np.arctan2(pos[n][1] - _c[1], pos[n][0] - _c[0]) % (2 * np.pi)

    for n in sorted(weak_singletons, key=cur_ang):
        a = cur_ang(n)
        best_i, best_d = None, None
        for i, m in enumerate(mids):
            if used[i]:
                continue
            d = abs((m - a + np.pi) % (2 * np.pi) - np.pi)
            if best_d is None or d < best_d:
                best_d, best_i = d, i
        used[best_i] = True
        pos[n] = ring_pos(mids[best_i])
elif weak_singletons:
    # without isolates, fall back to an even spread around the main component
    for i, n in enumerate(sorted(weak_singletons)):
        pos[n] = ring_pos(np.pi / 2 + 2 * np.pi * i / len(weak_singletons))

pd.DataFrame(
    [{"name": n, "x": pos[n][0], "y": pos[n][1]} for n in G.nodes]
).to_csv(os.path.join(NET_DIR, "layout_positions.csv"), index=False)

node_names = list(G.nodes)
xy = np.array([pos[n] for n in node_names])
name_to_idx = {n: i for i, n in enumerate(node_names)}

def blend(color, alpha, bg="#ffffff"):
    """EPS does not support transparency: blend alpha into the white background
    to get the equivalent solid color"""
    c = np.array(matplotlib.colors.to_rgb(color))
    b = np.array(matplotlib.colors.to_rgb(bg))
    return matplotlib.colors.to_hex((1 - alpha) * b + alpha * c)

def lighten(color, keep=LIGHT_KEEP, bg="#ffffff"):
    """Light color: blend the color with white at the given keep fraction"""
    return blend(color, keep, bg)

edge_style = LIGHT_EDGE_STYLE if LIGHT else EDGE_STYLE
suffix_out = "_light" if LIGHT else ""

# ---------------- Plotting ----------------
for suffix, keep in VERSIONS.items():
    fig, ax = plt.subplots(figsize=FIGSIZE_IN)
    ax.set_facecolor("white")

    # --- Edges ---
    segs, cols, lws = [], [], []
    for u, v, d in G.edges(data="hamming"):
        if d not in keep:
            continue
        c, w, a = edge_style[d]
        segs.append([pos[u], pos[v]])
        cols.append(blend(c, a))
        lws.append(w)
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws,
                                     capstyle="round", zorder=1))

    # --- Nodes ---
    # Singleton clusters (cluster id >= 3, one node each) are forced to gray so the
    # legend matches; set to False for the original colors
    USE_GRAY_SINGLETONS = True
    singleton_color = lighten("#9E9E9E", LIGHT_GRAY_KEEP) if LIGHT else "#9E9E9E"
    wt_i = node_names.index(WT)
    reg = [i for i in range(len(node_names)) if i != wt_i]
    face_reg = []
    for i in reg:
        n = node_names[i]
        if USE_GRAY_SINGLETONS and G.nodes[n]["cluster"] >= 3:
            face_reg.append(singleton_color)
        else:
            face_reg.append(lighten(G.nodes[n]["color"]) if LIGHT else G.nodes[n]["color"])
    ax.scatter(xy[reg, 0], xy[reg, 1], s=NODE_SIZE, c=face_reg,
               edgecolors="white", linewidths=0.4, zorder=3)
    # WT node drawn above ordinary labels so the gold highlight is never covered
    ax.scatter([xy[wt_i, 0]], [xy[wt_i, 1]], s=WT_NODE_SIZE, c="#FFD54F",
               edgecolors="#000000", linewidths=1.1, zorder=4.6)

    # --- Labels (the 46 from the original figure + WT) ---
    texts = []
    for n in node_names:
        if n == WT:
            continue
        if G.nodes[n]["label_show"]:
            texts.append(ax.text(pos[n][0], pos[n][1], n, fontsize=LABEL_FONTSIZE,
                                 color="#212121", zorder=4,
                                 path_effects=[pe.withStroke(linewidth=1.6, foreground="white")]))
    # adjustText>=1.4: avoidance points must be passed as separate x= / y= arrays,
    # otherwise they are silently ignored. Extra guard ring around the WT node keeps
    # labels off the enlarged gold node
    avoid_xy = xy.copy()
    wt_xy = np.array(pos[WT])
    for r in (0.02, 0.05, 0.10):
        ang = np.linspace(0, 2 * np.pi, 16, endpoint=False)
        avoid_xy = np.vstack([avoid_xy, wt_xy + r * np.column_stack([np.cos(ang), np.sin(ang)])])
    # WT arrow annotation (light version): the text box and arrow path also join the
    # avoidance set so peptide labels don't cover the amber text
    ann_x = pos[WT][0] + ANN_DX
    ann_y = pos[WT][1] + ANN_DY
    if LIGHT:
        for fx in np.linspace(-0.040, 0.003, 7):
            for fy in (-0.004, 0.0, 0.004):
                avoid_xy = np.vstack([avoid_xy, [ann_x + fx, ann_y + fy]])
        for t in np.linspace(0.2, 0.8, 4):
            avoid_xy = np.vstack([avoid_xy,
                                  [pos[WT][0] + t * (ann_x - pos[WT][0]),
                                   pos[WT][1] + t * (ann_y - pos[WT][1])]])
    if texts:
        adjust_text(texts, x=avoid_xy[:, 0], y=avoid_xy[:, 1], ax=ax,
                    expand=(1.15, 1.35),
                    force_text=(0.35, 0.55), force_static=(0.25, 0.45),
                    avoid_self=True, ensure_inside_axes=True)
    if not LIGHT:
        ax.text(pos[WT][0], pos[WT][1] + 0.035, WT, fontsize=WT_LABEL_FONTSIZE,
                fontweight="bold", ha="center", color="#000000", zorder=5,
                path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])

    # --- WT arrow annotation (light version): amber arrow from the lower left
    # pointing at the gold WT node, labelled with the full WT sequence ---
    if LIGHT:
        ann = ax.annotate(
            f"WT {WT}", xy=(pos[WT][0], pos[WT][1]),
            xytext=(ann_x, ann_y),
            fontsize=7.5, fontweight="bold", color=WT_ARROW_COLOR,
            ha="right", va="center", zorder=5.6,
            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
            arrowprops=dict(arrowstyle="-|>", color=WT_ARROW_COLOR,
                            lw=1.1, mutation_scale=13, shrinkA=2, shrinkB=7),
            annotation_clip=False)

    # --- Legend ---
    edge_handles = [
        Line2D([0], [0], color=blend(edge_style[d][0], edge_style[d][2]),
               lw=edge_style[d][1], label=f"Hamming distance = {d}")
        for d in keep
    ]
    c1, c2, c3 = ("#673AB7", "#F44336", "#E91E63")
    if LIGHT:
        c1, c2, c3 = lighten(c1), lighten(c2), lighten(c3)
    cluster_handles = [
        Patch(facecolor=c1, label="Cluster 1 (n = 94)"),
        Patch(facecolor=c2, label="Cluster 2 (n = 18)"),
        Patch(facecolor=c3, label="Cluster 3 (n = 4)"),
        Patch(facecolor=singleton_color, label="Singleton clusters (n = 9)"),
        Patch(facecolor="#FFD54F", edgecolor="#000000", label=f"WT {WT}"),
    ]
    # Wide canvas: the left side is occupied by the main cluster, so both legends
    # stack vertically in the free upper-right area (above the red/pink clusters)
    leg1 = ax.legend(handles=edge_handles, title="Edges", loc="upper right",
                     frameon=True, fontsize=6, title_fontsize=6.5,
                     borderpad=0.6, handlelength=1.8, framealpha=1.0)
    leg1.get_frame().set_edgecolor("#cccccc")
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=cluster_handles, title="Nodes", loc="upper right",
                     bbox_to_anchor=(1.0, 0.84),
                     frameon=True, fontsize=6, title_fontsize=6.5,
                     borderpad=0.6, handlelength=1.2, framealpha=1.0)
    leg2.get_frame().set_edgecolor("#cccccc")

    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.margins(0.03)

    for ext in ("pdf", "svg", "eps", "png"):
        out = os.path.join(NET_DIR, f"peptide_network_{suffix}{suffix_out}.{ext}")
        fig.savefig(out, format=ext, dpi=600 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.02)
        print("saved", out)
    plt.close(fig)

print("done")
