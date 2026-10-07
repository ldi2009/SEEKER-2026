# -*- coding: utf-8 -*-
"""Mimicry network, red-edge emphasis edition (plot_mimicry_red.py).
- Same source network as the main figure: all activated peptides, main-figure
  coordinates (nodes_layout_lev_*_final.csv), cluster colors cut=6.0 tab20,
  node shape = origin
- Visual hierarchy (weaken ordinary edges, highlight mimicry in red):
  red edge = mimicry edge (cross-origin virus<->self with d<=3), #d32f2f,
  line width decays with d; light gray edges = all remaining d<=4 edges
  (same-origin + d=4 cross-origin), heavily weakened
- Participant node outlines colored by origin (virus #d32f2f / human self
  #1565c0, 1.0 pt), labels in the same color family (virus dark red / self
  blue); all other nodes get thin outlines. Node sizes reduced: 34->26 /
  singleton 18->13
- Labels: every mimicry participant is annotated (cross-origin edge with
  d<=3, no gene-name deduplication). Gene edition is two-line: sequence on
  top, gene name below; virus gene names use the current UniProt entry
  (looked up by fasta accession in uniprot_virus_current.json, batch-checked
  2026-09, e.g. Q9ENL4 old VP2_C -> current CAPSD_CTFVL), self peptides use
  the Gene field
- Outward bias: labels preferentially placed on the side of the node facing
  away from the network centroid, reducing leader-line crossings
- Every label carries a dark-red dashed leader line (#a85757, dash 3-2)
  pointing to its node (user-approved final style)
- Exports {tcr}_labelpos_d3.json: exact node/edge/label positions for the
  proofreading/editor HTML (both gene and nogene label sets); if
  labelpos_override_d3.json exists (exported from the drag editor), those
  user-edited positions win
- The nogene edition is computed for the editor JSON only; figures are
  saved for the gene edition
Outputs: tcr{59,63}_mimicry_red_d3_gene.{pdf,svg,eps,png}
"""
import json
import os
import re
import time

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))          # -> B27_single_TCR_screen
DATA = os.path.join(REPO, "data", "mimicry_network")   # curated layouts + JSONs
RAW = os.path.join(REPO, "data", "results")            # peptide_screen outputs
OUT = os.path.join(REPO, "figures", "mimicry_network")  # generated figures
os.makedirs(OUT, exist_ok=True)
JSON_PATH = os.path.join(DATA, "tcr5963_layout_TCR59.json")
UNIPROT_CUR = json.load(open(os.path.join(DATA, "uniprot_virus_current.json"),
                             encoding="utf-8"))

FIGSIZE_IN = (6.5, 5.2)
NODE_SIZE = 26          # reduced (was 34)
TRI_FACTOR = 1.45
SINGLE_SIZE = 13        # reduced (was 18)
LABEL_FONTSIZE = 6.0
TITLE_FONTSIZE = 11
LIGHT_KEEP = 0.55
GRAY_SINGLE = "#c9c9c9"
MIM_K = 3

RED = "#d32f2f"
RED_LW = {1: 1.6, 2: 1.3, 3: 1.1}          # mimicry edge width decays with d
WEAK_EDGE = {1: ("#b0b0b0", 0.8), 2: ("#bdbdbd", 0.55),
             3: ("#cccccc", 0.30), 4: ("#d8d8d8", 0.20)}  # weakened ordinary edges
LAB_VIR = "#b71c1c"    # virus label dark red
LAB_SELF = "#1565c0"   # human self label blue
LINE_VIR = "#a85757"   # virus label leader line
LINE_SELF = "#7d94b5"  # self label leader line

CSV = {
    "TCR59": os.path.join(RAW, "TCR59_R1R2_annotated_with_R0_matched_filtered_with_primer_add_flow.csv"),
    "TCR63": os.path.join(RAW, "TCR63_R1R2_annotated_with_R0_matched_filtered_with_primer_add_flow.csv"),
}

plt.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})


def blend(color, keep, bg="#ffffff"):
    c = np.array(matplotlib.colors.to_rgb(color))
    b = np.array(matplotlib.colors.to_rgb(bg))
    return matplotlib.colors.to_hex((1 - keep) * b + keep * c)


def save_retry(fig, path, **kw):
    for attempt in range(4):
        try:
            fig.savefig(path, **kw)
            return
        except PermissionError:
            if attempt == 3:
                raise
            time.sleep(1.2)


def lev(p, q):
    if p == q:
        return 0
    lp, lq = len(p), len(q)
    prev = list(range(lq + 1))
    for i in range(1, lp + 1):
        cur = [i] + [0] * lq
        pi = p[i - 1]
        for j in range(1, lq + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1,
                         prev[j - 1] + (pi != q[j - 1]))
        prev = cur
    return prev[lq]


def load(f):
    df = pd.read_csv(f)
    flow = df[["aa", "origin", "Gene", "fasta_name",
               "add_flow_single_copy"]].dropna(
        subset=["add_flow_single_copy"]).copy()
    flow = flow.groupby("aa").agg(
        origin=("origin", "first"), gene=("Gene", "first"),
        fasta=("fasta_name", "first"),
        flow=("add_flow_single_copy", "max")).reset_index()
    return flow[flow.flow > 5].reset_index(drop=True)


def sane_gene(g, origin, fasta=""):
    """virus: look up the current UniProt entry by accession (gene_species,
    checked 2026-09); every entry inside the fasta is outdated (truncated or
    renamed species codes), e.g. Q9ENL4 old VP2_C -> current CAPSD_CTFVL;
    self: the Gene field; placeholder name when neither has information"""
    if origin == "virus":
        parts = str(fasta).split("+")
        if len(parts) >= 2 and parts[1].startswith(("sp_", "tr")):
            acc = parts[1].split("_")[1]
            hit = UNIPROT_CUR.get(acc, {}).get("entry", "")
            if hit and hit != "?":
                return hit
        m = re.match(r"^[^+]*\+sp_[A-Za-z0-9]+_(.+?)(\+.*)?$", str(fasta))
        if m:
            e = m.group(1).strip("_")
            if e:
                return e
    g0 = str(g).strip()
    if len(g0) <= 2 or g0.replace(".", "").replace("-", "").isdigit():
        return "virus" if origin == "virus" else "self"
    return g0


def render(tcr, sub, D, G, pos, th, cl_color, pal3, n_sv, n_ss, out_base,
           label_gene, dpi_png=600, label_override=None, save=True):
    n = len(sub)
    deg = dict(G.degree())
    sub = sub.copy()
    sub["deg"] = [deg.get(i, 0) for i in range(n)]

    fig, ax = plt.subplots(figsize=FIGSIZE_IN)
    ax.set_position([0, 0, 1, 1])

    # edges split in two groups: red mimicry (zorder 2.5) + weakened ordinary (zorder 2)
    red_segs, red_lws, gray_segs, gray_cols, gray_lws = [], [], [], [], []
    for i in range(n):
        for j in range(i + 1, n):
            d = int(D[i, j])
            if d > th:
                continue
            cross = sub.loc[i, "origin"] != sub.loc[j, "origin"]
            if cross and d <= MIM_K:
                red_segs.append([pos[i], pos[j]])
                red_lws.append(RED_LW[d])
            else:
                c, lw = WEAK_EDGE[d]
                gray_segs.append([pos[i], pos[j]])
                gray_cols.append(c)
                gray_lws.append(lw)
    ax.add_collection(LineCollection(gray_segs, colors=gray_cols,
                                     linewidths=gray_lws, zorder=2))
    ax.add_collection(LineCollection(red_segs, colors=RED,
                                     linewidths=red_lws, zorder=2.5,
                                     capstyle="round"))

    for i in range(n):
        r = sub.loc[i]
        marker = "o" if r.origin == "virus" else "^"
        size = NODE_SIZE if r.origin == "virus" else NODE_SIZE * TRI_FACTOR
        color = cl_color.get(r.cluster)
        mim = r.cross_deg > 0
        if color is None:
            ax.scatter(*pos[i], s=SINGLE_SIZE, marker=marker, c=GRAY_SINGLE,
                       edgecolors="#8a8a8a", linewidths=0.4, zorder=3.0)
        elif mim:
            stroke_m = RED if r.origin == "virus" else LAB_SELF
            ax.scatter(*pos[i], s=size, marker=marker, c=blend(color, LIGHT_KEEP),
                       edgecolors=stroke_m, linewidths=1.0, zorder=3.3)
        else:
            ax.scatter(*pos[i], s=size, marker=marker, c=blend(color, LIGHT_KEEP),
                       edgecolors=blend(color, 0.9), linewidths=0.6, zorder=3.0)

    all_xy = np.array([pos[i] for i in range(n)])
    top_y = float(all_xy[:, 1].max())
    ax.text(0, top_y + 0.05, tcr, fontsize=TITLE_FONTSIZE, fontweight="bold",
            color="#555555", ha="center", va="bottom", zorder=4)

    mx = float(np.abs(all_xy[:, 0]).max()) + 0.10
    my = float(np.abs(all_xy[:, 1]).max()) + 0.20
    S_pt = min(468.0 / (2 * mx), 374.4 / (2 * my))
    ax.set_xlim(-mx, mx)
    ax.set_ylim(-my, my)

    ann_boxes = []
    tw = len(tcr) * TITLE_FONTSIZE * 0.62
    th_pt = 13.5
    ann_boxes.append((-tw / (2 * S_pt), top_y + 0.05,
                      tw / (2 * S_pt), top_y + 0.05 + th_pt / S_pt))

    def box_dist(b1, b2):
        dx = max(b2[0] - b1[2], b1[0] - b2[2], 0)
        dy = max(b2[1] - b1[3], b1[1] - b2[3], 0)
        return np.hypot(dx, dy) * S_pt

    cand = []
    order = sorted(sub.index, key=lambda k: (-sub.loc[k, "cross_deg"],
                                             -sub.loc[k, "deg"]))
    for k in order:
        if sub.loc[k, "cross_deg"] > 0:
            cand.append(k)          # label every participant, no gene dedup
    cand.sort(key=lambda k: -sub.loc[k, "deg"])

    n_skip = n_free = n_line = 0
    label_anns = []
    label_info = []
    labels_pos = []

    FREE_RADII = (13, 16, 19, 23)
    LINE_RADII = (26, 36, 48, 64, 86, 110, 128)
    centroid = all_xy.mean(axis=0)

    for i in cand:
        r = sub.loc[i]
        if label_gene:
            lines = [r.aa, r.gene_sane]
        else:
            lines = [r.aa]
        label = "\n".join(lines)
        nx_, ny_ = pos[i]
        w_pt = max(len(t) for t in lines) * LABEL_FONTSIZE * 0.58
        h_pt = LABEL_FONTSIZE * 1.12 * len(lines) + 0.3
        others_xy = np.array([pos[j] for j in range(n) if j != i])

        vx, vy = nx_ - centroid[0], ny_ - centroid[1]
        vnorm = float(np.hypot(vx, vy))
        out_ang = float(np.arctan2(vy, vx)) if vnorm > 0.02 else None

        ov = (label_override or {}).get(int(i))
        if ov is not None:
            best = (1.0, float(ov["cx"]), float(ov["cy"]), bool(ov["ha"]),
                    99.0, 0.0)
            line_free = True
        else:
            def search(radii, min_cl_other, free):
                best = None
                for radius_pt in radii:
                    for k in range(16):
                        ang = 2 * np.pi * k / 16 + 0.11
                        cx = nx_ + np.cos(ang) * radius_pt / S_pt
                        cy = ny_ + np.sin(ang) * radius_pt / S_pt
                        x0 = cx if np.cos(ang) >= 0 else cx - w_pt / S_pt
                        y0 = cy - h_pt / (2 * S_pt)
                        box = (x0, y0, x0 + w_pt / S_pt, y0 + h_pt / S_pt)
                        if box[0] < -mx + 0.02 or box[2] > mx - 0.02 \
                                or box[1] < -my + 0.02 or box[3] > my - 0.02:
                            continue
                        dxa = max(box[0] - nx_, nx_ - box[2], 0)
                        dya = max(box[1] - ny_, ny_ - box[3], 0)
                        cl_anchor = float(np.hypot(dxa, dya)) * S_pt
                        if cl_anchor < 6:
                            continue
                        dx = np.maximum(np.maximum(others_xy[:, 0] - box[0],
                                                   box[2] - others_xy[:, 0]), 0)
                        dy = np.maximum(np.maximum(others_xy[:, 1] - box[1],
                                                   box[3] - others_xy[:, 1]), 0)
                        cl_other = float(np.min(np.hypot(dx, dy))) * S_pt
                        if cl_other < min_cl_other:
                            continue
                        if free and cl_anchor >= cl_other:
                            continue
                        cl_lab = min((box_dist(box, b) for b in ann_boxes), default=999)
                        if cl_lab < 5.0:
                            continue
                        score = min(cl_other, 42.0) + (8.0 if free else 0.0) \
                            - 0.08 * radius_pt
                        if out_ang is not None:
                            score += 6.0 * np.cos(ang - out_ang)
                        if best is None or score > best[0]:
                            best = (score, cx, cy, np.cos(ang) >= 0, cl_other,
                                    radius_pt)
                return best

            best = search(FREE_RADII, 7.5, True)
            line_free = best is not None
            if not line_free:
                best = search(LINE_RADII, 0.0, False)
        if best is None:
            n_skip += 1
            continue
        _, cx, cy, ha, cl_out, radius_pt = best
        is_vir = sub.loc[i, "origin"] == "virus"
        aps = dict(arrowstyle="-", color=LINE_VIR if is_vir else LINE_SELF,
                   lw=0.7, linestyle=(0, (3, 2)), shrinkA=1.5, shrinkB=4)
        ann = ax.annotate(label, xy=(nx_, ny_), xytext=(cx, cy),
                          fontsize=LABEL_FONTSIZE,
                          color=LAB_VIR if is_vir else LAB_SELF,
                          ha="left" if ha else "right", va="center", zorder=5,
                          path_effects=[pe.withStroke(linewidth=1.6, foreground="white")],
                          arrowprops=aps,
                          annotation_clip=False)
        label_anns.append(ann)
        label_info.append((label.replace("\n", " / "), cl_out))
        labels_pos.append(dict(i=int(i), lines=lines,
                               cx=round(float(cx), 5), cy=round(float(cy), 5),
                               ha=int(ha), free=int(line_free)))
        if line_free:
            n_free += 1
        else:
            n_line += 1
        w_d = w_pt / S_pt
        h_d = h_pt / S_pt
        x0 = cx if ha else cx - w_d
        ann_boxes.append((x0, cy - h_d / 2, x0 + w_d, cy + h_d / 2))

    edge_handles = [
        Line2D([0], [0], color=RED, lw=1.4,
               label="Mimicry edge (virus\u2013self, Levenshtein d \u2264 3)"),
        Line2D([0], [0], color="#bdbdbd", lw=0.6,
               label="Other edges (Levenshtein d = 1\u20134)"),
    ]
    leg1 = ax.legend(handles=edge_handles, title="Edges", loc="center left",
                     bbox_to_anchor=(1.01, 0.70), frameon=True, fontsize=6,
                     title_fontsize=6.5, borderpad=0.6, handlelength=1.8,
                     framealpha=1.0)
    leg1.get_frame().set_edgecolor("#cccccc")
    ax.add_artist(leg1)

    node_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.5,
               markerfacecolor="#ffffff", markeredgecolor="#8a8a8a",
               markeredgewidth=0.7),
        Line2D([0], [0], marker="^", linestyle="none", markersize=6.5,
               markerfacecolor="#ffffff", markeredgecolor="#8a8a8a",
               markeredgewidth=0.7),
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.5,
               markerfacecolor="#d2d2d2", markeredgecolor=RED,
               markeredgewidth=1.0),
        Line2D([0], [0], marker="^", linestyle="none", markersize=6.5,
               markerfacecolor="#d2d2d2", markeredgecolor=LAB_SELF,
               markeredgewidth=1.0),
    ]
    node_labels = ["Virus antigen", "Self antigen",
                   "Mimicry participant (virus)",
                   "Mimicry participant (self)"]
    if n_sv > 0:
        node_handles.append(Line2D([0], [0], marker="o", linestyle="none",
                                   markersize=4.0, markerfacecolor=GRAY_SINGLE,
                                   markeredgecolor="#8a8a8a", markeredgewidth=0.4))
        node_labels.append(f"Virus singleton (n = {n_sv})")
    if n_ss > 0:
        node_handles.append(Line2D([0], [0], marker="^", linestyle="none",
                                   markersize=4.0, markerfacecolor=GRAY_SINGLE,
                                   markeredgecolor="#8a8a8a", markeredgewidth=0.4))
        node_labels.append(f"Self singleton (n = {n_ss})")
    leg2 = ax.legend(handles=node_handles, labels=node_labels, title="Nodes",
                     loc="center left", bbox_to_anchor=(1.01, 0.24), frameon=True,
                     fontsize=6, title_fontsize=6.5, borderpad=0.6, handlelength=1.5,
                     handletextpad=0.7, labelspacing=0.8, framealpha=1.0)
    leg2.get_frame().set_edgecolor("#cccccc")

    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_aspect("equal", adjustable="box")

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # leg1 added via add_artist is not part of the default tight bbox ->
    # legend text would be clipped on the right, so pass it explicitly
    extra = (leg1, leg2)
    tb = fig.get_tightbbox(renderer, bbox_extra_artists=extra)
    pad_in = 0.02
    boxes = []
    for a in label_anns:
        e = a.get_window_extent(renderer)
        boxes.append(dict(
            text=a.get_text(),
            x0=round((e.x0 / fig.dpi - (tb.x0 - pad_in)) * 600, 1),
            x1=round((e.x1 / fig.dpi - (tb.x0 - pad_in)) * 600, 1),
            ytop=round(((tb.y1 + pad_in) - e.y1 / fig.dpi) * 600, 1),
            ybot=round(((tb.y1 + pad_in) - e.y0 / fig.dpi) * 600, 1)))
    if save:
        with open(out_base + "_labelboxes.json", "w", encoding="utf-8") as f:
            json.dump(dict(boxes=boxes, label_info=[[t, round(c, 1)] for t, c in label_info],
                           n_red_edges=len(red_segs)),
                      f, ensure_ascii=False, indent=1)

        for ext in ("pdf", "svg", "eps", "png"):
            save_retry(fig, out_base + "." + ext,
                       dpi=dpi_png if ext == "png" else None,
                       bbox_inches="tight", pad_inches=0.02,
                       bbox_extra_artists=extra)
    plt.close(fig)

    # payload for the proofreading HTML: node styles / edges / exact label positions
    nodes_js, edges_js = [], []
    for i in range(n):
        r = sub.loc[i]
        color = cl_color.get(r.cluster)
        mim = r.cross_deg > 0
        if color is None:
            fill, stroke, slw = GRAY_SINGLE, "#8a8a8a", 0.4
            r_pt = (SINGLE_SIZE ** 0.5) / 2
        else:
            fill = blend(color, LIGHT_KEEP)
            if mim:
                stroke = RED if r.origin == "virus" else LAB_SELF
                slw = 1.0
            else:
                stroke, slw = blend(color, 0.9), 0.6
            sz = NODE_SIZE if r.origin == "virus" else NODE_SIZE * TRI_FACTOR
            r_pt = (sz ** 0.5) / 2
        nodes_js.append(dict(aa=r.aa, o=r.origin, g=r.gene_sane,
                             x=round(pos[i][0], 5), y=round(pos[i][1], 5),
                             tri=int(r.origin != "virus"), fill=fill,
                             stroke=stroke, slw=slw, r=round(r_pt, 3),
                             mim=int(mim), x3=int(r.cross_deg),
                             cl=int(r.cluster)))
    for i in range(n):
        for j in range(i + 1, n):
            d = int(D[i, j])
            if d > th:
                continue
            cross = sub.loc[i, "origin"] != sub.loc[j, "origin"]
            edges_js.append([i, j, d, int(cross and d <= MIM_K)])
    payload = dict(nodes=nodes_js, edges=edges_js, labels=labels_pos,
                   geo=dict(S=round(S_pt, 2), mx=round(mx, 4), my=round(my, 4)))
    return len(red_segs), n_free, n_line, n_skip, payload


cfg = json.load(open(JSON_PATH, encoding="utf-8"))
th = int(cfg["params"]["dispTh"])
cut = float(cfg["params"]["cut"])
print(f"JSON: dispTh=d<={th} cut={cut} | mimicry: cross-origin d<={MIM_K}")

for tcr, csv_path in CSV.items():
    sub = load(csv_path)
    seqs = sub.aa.tolist()
    org = (sub.origin == "virus").tolist()
    n = len(seqs)
    D = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            d = lev(seqs[i], seqs[j])
            D[i, j] = D[j, i] = d

    Z = linkage(squareform(D.astype(float)), method="average")
    sub["cluster"] = fcluster(Z, t=cut, criterion="distance")
    csizes = sub.groupby("cluster").size()
    multi = sorted([c for c in csizes.index if csizes[c] >= 2],
                   key=lambda c: (-csizes[c], c))
    palette = list(plt.cm.tab20.colors)
    cl_color = {c: palette[k % len(palette)] for k, c in enumerate(multi)}
    sub["gene_sane"] = [sane_gene(g, o, fa) for g, o, fa
                        in zip(sub.gene, sub.origin, sub.fasta)]

    cross_deg = np.zeros(n, dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] <= MIM_K and org[i] != org[j]:
                cross_deg[i] += 1
                cross_deg[j] += 1
    sub["cross_deg"] = cross_deg

    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] <= th:
                G.add_edge(i, j)

    co = pd.read_csv(os.path.join(DATA, f"nodes_layout_lev_{tcr}_final.csv"))
    assert co.aa.tolist() == seqs, "layout CSV does not align with the data"
    xy = co[["x", "y"]].to_numpy(dtype=float)
    pos = {i: (xy[k, 0], xy[k, 1]) for k, i in enumerate(range(n))}

    n_mim = int((cross_deg > 0).sum())
    sing = sub.cluster.map(lambda c: csizes[c] == 1)
    n_sv = int((sing & (sub.origin == "virus")).sum())
    n_ss = int((sing & (sub.origin != "virus")).sum())
    print(f"===== {tcr}: n={n}, participants {n_mim}, red (cross-origin "
          f"d<={MIM_K}) edges {int(cross_deg.sum() // 2)} =====")

    ov_path = os.path.join(DATA, "labelpos_override_d3.json")
    ov_all = {}
    if os.path.exists(ov_path):
        ov_all = json.load(open(ov_path, encoding="utf-8"))
        print(f"  [override] using user-edited positions {ov_path}")
    payloads = {}
    for label_gene, tag in ((True, "gene"), (False, "nogene")):
        ovr = {int(e["i"]): e for e in
               (ov_all.get(tcr, {}).get(f"labels_{tag}", []))}
        out_base = os.path.join(OUT, f"{tcr.lower()}_mimicry_red_d3_{tag}")
        n_red, n_free, n_line, n_sk, pl = render(
            tcr, sub, D, G, pos, th, cl_color, palette[:3],
            n_sv, n_ss, out_base, label_gene, label_override=ovr or None,
            save=label_gene)
        payloads[tag] = pl
        if label_gene:
            print(f"  saved ..._d3_{tag}.* (red edges {n_red}, labels "
                  f"{n_free + n_line + n_sk} all with leader lines, "
                  f"skipped {n_sk})")
        else:
            print("  nogene label positions computed for the editor JSON "
                  "(figures not saved)")
    with open(os.path.join(OUT, f"{tcr.lower()}_labelpos_d3.json"), "w",
              encoding="utf-8") as f:
        json.dump(dict(nodes=payloads["gene"]["nodes"],
                       edges=payloads["gene"]["edges"],
                       geo=payloads["gene"]["geo"],
                       labels_gene=payloads["gene"]["labels"],
                       labels_nogene=payloads["nogene"]["labels"]),
                  f, ensure_ascii=False)

print("done")
