#!/usr/bin/env python3
"""Reproduce the Figure 4 idaOEPCR ranking and distribution panels.

For each requested group (default: Not, standard) two figures are produced:

* ``{group}_ranking_distribution.pdf``
  Four stacked panels, one per replicate, showing the rank-abundance
  distribution (Rank vs log2(Count)). Unique peptides are drawn in blue and
  all remaining peptides in orange.
* ``{group}_new_combined_distribution.png``
  The four replicates overlaid as a log10(read count) histogram of unique
  peptide counts.

Unique-peptide calling (a peptide is unique when no peptide with a strictly
higher Count lies within a Hamming distance of 2) is an exact but faster
equivalent of the original per-row loop: for every row only the
Hamming-distance-2 neighbourhood of its sequence is intersected with the set
of sequences already seen at higher counts.

The peptide nucleotide sequence is read from the ``peptideDNA`` column.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NUCLEOTIDES = "ACGT"
NEIGHBOR_CACHE = {}


def neighbors(seq):
    """All sequences at Hamming distance <= 2 from seq, including seq itself."""
    cached = NEIGHBOR_CACHE.get(seq)
    if cached is not None:
        return cached
    chars = list(seq)
    n = len(chars)
    out = {seq}
    for i in range(n):
        orig_i = chars[i]
        for a in NUCLEOTIDES:
            if a == orig_i:
                continue
            chars[i] = a
            out.add("".join(chars))
            for j in range(i + 1, n):
                orig_j = chars[j]
                for b in NUCLEOTIDES:
                    if b == orig_j:
                        continue
                    chars[j] = b
                    out.add("".join(chars))
                chars[j] = orig_j
        chars[i] = orig_i
    NEIGHBOR_CACHE[seq] = out
    return out


def unique_flags(df):
    """Mark unique peptides; df must be sorted by Count descending."""
    seqs = df["peptideDNA"].astype(str).values
    counts = df["Count"].values
    n = len(df)
    flags = np.empty(n, dtype=bool)
    seen = set()
    i = 0
    while i < n:
        j = i
        while j < n and counts[j] == counts[i]:
            j += 1
        for k in range(i, j):
            flags[k] = not any(x in seen for x in neighbors(seqs[k]))
        for k in range(i, j):
            seen.add(seqs[k])
        i = j
    return flags


def filter_data(df):
    """Apply the abundance cut-off and annotate ranking / uniqueness."""
    threshold = df["Count"].max() / 300
    filtered = df[df["Count"] >= threshold].copy()
    if filtered.empty:
        return filtered
    filtered["Ranking"] = np.arange(1, len(filtered) + 1)
    filtered["log2Count"] = np.log2(filtered["Count"])
    filtered["unique_peptide"] = unique_flags(filtered)
    return filtered


def group_of(filename):
    for keyword in ("NE", "Not", "standard", "heat"):
        if keyword.lower() in filename.lower():
            return keyword
    return "unknown"


def calculate_global_breaks(directory, csv_files):
    """Shared log10 bin edges across all samples."""
    all_log = []
    for filename in csv_files:
        df = pd.read_csv(os.path.join(directory, filename))
        unique_df = filter_data(df)
        if unique_df.empty:
            continue
        count_table = unique_df[unique_df["unique_peptide"]]["Count"].value_counts()
        valid = count_table[count_table.index > 0].index.dropna()
        if not valid.empty:
            all_log.extend(np.log10(valid))
    if not all_log:
        return np.linspace(1.5, 4.0, num=30)
    return np.linspace(min(all_log) - 0.1, max(all_log) + 0.1, num=30)


def bin_sums_for(directory, filename, breaks):
    """Per-bin sum of unique peptides for one sample."""
    df = pd.read_csv(os.path.join(directory, filename))
    unique_df = filter_data(df)
    if unique_df.empty:
        return None
    count_table = unique_df[unique_df["unique_peptide"]]["Count"].value_counts()
    valid = count_table[count_table.index > 0].index.dropna()
    if valid.empty:
        return None
    binned = pd.cut(np.log10(valid), bins=breaks, include_lowest=True)
    valid_count = count_table[count_table.index > 0].dropna()
    return valid_count.groupby(binned, observed=False).sum().fillna(0).values


def plot_ranking(directory, group, filenames, output_dir):
    fig, axes = plt.subplots(4, 1, figsize=(12, 24))
    for i, filename in enumerate(filenames[:4]):
        df = pd.read_csv(os.path.join(directory, filename)).copy()
        df["Ranking"] = np.arange(1, len(df) + 1)
        df["log2Count"] = np.log2(df["Count"])
        df["unique_peptide"] = unique_flags(df)

        short_label = filename.split("_top_10000")[0]
        unique_df = df[df["unique_peptide"]]
        other_df = df[~df["unique_peptide"]]

        axes[i].bar(unique_df["Ranking"], unique_df["log2Count"], alpha=0.7,
                    label=f"{short_label} - Unique", color="blue")
        axes[i].bar(other_df["Ranking"], other_df["log2Count"], alpha=0.7,
                    label=f"{short_label} - Other", color="orange")
        axes[i].set_title(f"File: {short_label} - Sequence Count Distribution",
                          fontsize=14)
        axes[i].set_xlabel("Ranking", fontsize=12)
        axes[i].set_ylabel("log2(Count)", fontsize=12)
        axes[i].set_ylim(5, 15)
        axes[i].set_xlim(0, 1200)
        axes[i].grid(False)
        axes[i].legend(fontsize=8, loc="upper right")

    plt.tight_layout()
    path = os.path.join(output_dir, f"{group}_ranking_distribution.pdf")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {path}")


def plot_combined(group, sums, breaks, labels, output_dir):
    plt.figure(figsize=(8, 6))
    bin_centers = (breaks[:-1] + breaks[1:]) / 2
    colors = ["red", "green", "blue", "orange"]
    for i, (label, values) in enumerate(zip(labels, sums)):
        plt.bar(bin_centers, values, width=(breaks[1] - breaks[0]),
                edgecolor="darkblue", color=colors[i % len(colors)],
                alpha=0.5, label=label)
    plt.xlabel("log10(Read Count)")
    plt.ylabel("Unique Peptide Count Sum")
    plt.title(f"{group} group: Unique Peptide Count Distribution")
    plt.ylim(0, 50)
    plt.legend()
    path = os.path.join(output_dir, f"{group}_new_combined_distribution.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True,
                        help="directory with *_top_10000_blast_sequences.csv")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--groups", nargs="+", default=["Not", "standard"])
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    csv_files = sorted(f for f in os.listdir(args.input_dir)
                       if "blast" in f and f.endswith(".csv"))
    if not csv_files:
        raise SystemExit("No *_top_10000_blast_sequences.csv found.")

    breaks = calculate_global_breaks(args.input_dir, csv_files)

    for group in args.groups:
        group_files = [f for f in csv_files if group_of(f) == group]
        if not group_files:
            print(f"No files found for group {group}, skipping.")
            continue

        plot_ranking(args.input_dir, group, group_files, args.output_dir)

        labels, sums = [], []
        for filename in group_files:
            values = bin_sums_for(args.input_dir, filename, breaks)
            if values is None:
                continue
            labels.append(filename.split("_top_10000")[0])
            sums.append(values)
        if sums:
            plot_combined(group, sums, breaks, labels, args.output_dir)


if __name__ == "__main__":
    main()