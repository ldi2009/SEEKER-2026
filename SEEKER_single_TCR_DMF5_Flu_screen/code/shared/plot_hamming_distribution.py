# -*- coding: utf-8 -*-
"""
Flu vs DMF5 peptide library pairwise Hamming distance distribution comparison (publication vector figure)
(a) Percentage of peptide pairs per distance bin (grouped bar chart)
(b) Cumulative percentage curves + dashed lines marking the network granularity cut-offs
Data: results/flu_network/nodes.csv (Flu, n=125) + results/dmf5_network/nodes.csv (n=295)
Output: hamming_distance_distribution_Flu_vs_DMF5.{pdf,svg,eps,png}
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

plt.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

AA = "ACDEFGHIKLMNPQRSTVWY"

def pair_dist_hist(nodes_csv):
    names = pd.read_csv(nodes_csv)["name"].tolist()
    arr = np.array([[AA.index(c) for c in s] for s in names])
    D = pdist(arr, metric="hamming") * len(names[0])
    counts = np.array([(D == d).sum() for d in range(1, 10)])
    return counts, counts.sum()

def blend(color, alpha, bg="#ffffff"):
    c = np.array(matplotlib.colors.to_rgb(color))
    b = np.array(matplotlib.colors.to_rgb(bg))
    return matplotlib.colors.to_hex((1 - alpha) * b + alpha * c)

flu_cnt, flu_tot = pair_dist_hist(os.path.join(REPO, "results", "flu_network", "nodes.csv"))
dmf5_cnt, dmf5_tot = pair_dist_hist(os.path.join(REPO, "results", "dmf5_network", "nodes.csv"))

flu_pct = flu_cnt / flu_tot * 100
dmf5_pct = dmf5_cnt / dmf5_tot * 100
flu_cum = np.cumsum(flu_pct)
dmf5_cum = np.cumsum(dmf5_pct)

print(f"Flu  n_pairs={flu_tot}: cum={np.round(flu_cum, 1)}")
print(f"DMF5 n_pairs={dmf5_tot}: cum={np.round(dmf5_cum, 1)}")

C_FLU = "#4C72B0"
C_DMF5 = "#DD8452"
DS = np.arange(1, 10)

fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.1))
ax1, ax2 = axes

# ---- (a) Distance distribution bar chart ----
w = 0.38
ax1.bar(DS - w / 2, flu_pct, width=w, color=blend(C_FLU, 0.85),
        edgecolor=C_FLU, linewidth=0.4, label="Flu-M1 mutant library (n = 125)")
ax1.bar(DS + w / 2, dmf5_pct, width=w, color=blend(C_DMF5, 0.85),
        edgecolor=C_DMF5, linewidth=0.4, label="DMF5 library (n = 295)")
ax1.set_xlabel("Hamming distance", fontsize=8)
ax1.set_ylabel("Peptide pairs (%)", fontsize=8)
ax1.set_xticks(DS)
ax1.set_xlim(0.4, 9.6)
ax1.set_ylim(0, 30)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.tick_params(labelsize=7)
ax1.text(0.02, 0.97, "a", transform=ax1.transAxes, fontsize=10, fontweight="bold",
         va="top", ha="left")

# ---- (b) Cumulative curves ----
# Cut-off dashed lines use pre-blended colors instead of alpha (EPS has no transparency)
C_FLU_D = blend(C_FLU, 0.55)
C_DMF5_D = blend(C_DMF5, 0.55)
ax2.axvline(2, color=C_FLU_D, ls="--", lw=1.1, zorder=1)
ax2.axvline(3, color=C_DMF5_D, ls="--", lw=1.1, zorder=1)
ax2.plot(DS, flu_cum, "-o", color=C_FLU, ms=3.5, lw=1.4, label="Flu-M1 mutant library")
ax2.plot(DS, dmf5_cum, "-s", color=C_DMF5, ms=3.5, lw=1.4, label="DMF5 library")
# Key-point values: Flu offsets point down-right (curve rises to the upper right,
# lower right is empty), DMF5 offsets point up-left (curve stays low, upper left is empty)
for d, off in {1: (10, -3), 2: (8, -8), 3: (8, -8), 4: (8, -8)}.items():
    ax2.annotate(f"{flu_cum[d-1]:.1f}", (d, flu_cum[d-1]), textcoords="offset points",
                 xytext=off, fontsize=5.5, color=C_FLU, ha="left", va="center")
for d in (2, 3, 4):
    ax2.annotate(f"{dmf5_cum[d-1]:.1f}", (d, dmf5_cum[d-1]), textcoords="offset points",
                 xytext=(-4, 2), fontsize=5.5, color=C_DMF5, ha="right", va="center")
ax2.set_xlabel("Hamming distance", fontsize=8)
ax2.set_ylabel("Cumulative pairs (%)", fontsize=8)
ax2.set_xticks(DS)
ax2.set_xlim(0.5, 9.5)
ax2.set_ylim(0, 108)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.tick_params(labelsize=7)
ax2.text(0.02, 0.97, "b", transform=ax2.transAxes, fontsize=10, fontweight="bold",
         va="top", ha="left")

# ---- Figure-level legend (shared by both panels + cut-off line explanation) ----
handles = [
    Line2D([0], [0], color=C_FLU, lw=1.4, marker="o", ms=3.5,
           label=f"Flu-M1 mutant library (n = 125, {flu_tot:,} pairs)"),
    Line2D([0], [0], color=C_DMF5, lw=1.4, marker="s", ms=3.5,
           label=f"DMF5 library (n = 295, {dmf5_tot:,} pairs)"),
    Line2D([0], [0], color=C_FLU_D, lw=1.2, ls="--",
           label="Flu network cut-off (d \u2264 2)"),
    Line2D([0], [0], color=C_DMF5_D, lw=1.2, ls="--",
           label="DMF5 network cut-off (d \u2264 3)"),
]
fig.legend(handles=handles, loc="upper center", ncol=2, fontsize=6.5,
           frameon=False, bbox_to_anchor=(0.5, 1.0), columnspacing=1.6,
           handlelength=2.2)

fig.tight_layout(rect=[0, 0, 1, 0.88])
for ext in ("pdf", "svg", "eps", "png"):
    out = os.path.join(REPO, "results", f"hamming_distance_distribution_Flu_vs_DMF5.{ext}")
    fig.savefig(out, format=ext, dpi=600 if ext == "png" else None)
plt.close(fig)
print("saved hamming_distance_distribution_Flu_vs_DMF5.*")
