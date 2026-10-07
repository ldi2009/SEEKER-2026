#!/usr/bin/env python3
"""Group-level unique-peptide binning and cumulative curves.

For every BLAST-annotated sample table the 1000 most abundant peptides are
kept and scanned in descending abundance order; a peptide is retained as
"unique" when it lies at a Hamming distance > 2 from every peptide already
retained. Retained peptides are grouped (heat / NE / Not / standard), binned on
a log10 read-count axis, and displayed as a bin histogram together with a
cumulative curve (per group and across all groups).

The peptide nucleotide sequence is read from the ``peptideDNA`` column.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return float("inf")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def extract_group(filename):
    for group in ["heat", "NE", "Not", "standard"]:
        if group in filename:
            return group
    return "other"


def create_log_bins(data, num_bins=10):
    min_val = data.min()
    max_val = data.max()
    if min_val <= 0:
        min_val = 1e-1
    return np.logspace(np.log10(min_val), np.log10(max_val), num_bins + 1)


def process_data(directory):
    all_samples = [
        f for f in os.listdir(directory) if "blast" in f and f.endswith(".csv")
    ]
    all_unique = []
    all_samples_unique = {}

    for filename in all_samples:
        sample_name = filename.split("_top_10000")[0]
        df = pd.read_csv(os.path.join(directory, filename)).nlargest(1000, "Count")

        unique_peptides = []
        for _, row in df.iterrows():
            seq = row["peptideDNA"]
            if all(hamming_distance(seq, p) > 2 for p in unique_peptides):
                unique_peptides.append(seq)
        all_samples_unique[sample_name] = set(unique_peptides)

    for filename in all_samples:
        sample_name = filename.split("_top_10000")[0]
        df = pd.read_csv(os.path.join(directory, filename)).nlargest(1000, "Count")

        unique_peptides = []
        for _, row in df.iterrows():
            seq = row["peptideDNA"]
            if all(hamming_distance(seq, p) > 2 for p in unique_peptides):
                unique_peptides.append(seq)
                all_unique.append({
                    "sample": sample_name,
                    "group": extract_group(sample_name),
                    "peptideDNA": seq,
                    "Count": row["Count"],
                    "log10_count": np.log10(row["Count"]) if row["Count"] > 0 else 0,
                    "is_repeated": any(
                        seq in seqs
                        for s, seqs in all_samples_unique.items()
                        if s != sample_name
                    ),
                })

    return pd.DataFrame(all_unique)


def plot_combined_bin_cumulative(data, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    groups = sorted(data["group"].unique())
    num_bins = 10

    for group in groups:
        group_data = data[data["group"] == group].copy()
        if group_data.empty:
            print(f"No data for group {group}, skipping")
            continue

        log_counts = group_data["log10_count"]
        bins = create_log_bins(log_counts, num_bins)
        group_data["bin"] = pd.cut(log_counts, bins=bins)

        bin_stats = (
            group_data.groupby("bin", observed=False).size().reset_index(name="count")
        )
        bin_centers = [(b.left + b.right) / 2 for b in bin_stats["bin"]]
        bin_labels = [f"{b.left:.1f}-{b.right:.1f}" for b in bin_stats["bin"]]

        sorted_log = np.sort(group_data["log10_count"])
        cumulative = np.arange(1, len(sorted_log) + 1)

        fig, ax1 = plt.subplots(figsize=(12, 7))
        ax1.bar(bin_centers, bin_stats["count"],
                width=(bins[1] - bins[0]) * 0.8,
                color="skyblue", edgecolor="black", alpha=0.7,
                label=f"{group} bin counts")
        ax1.set_xlabel("log10(read count)")
        ax1.set_ylabel("Unique peptides in bin", color="blue")
        ax1.tick_params(axis="y", labelcolor="blue")
        ax1.set_xticks(bin_centers)
        ax1.set_xticklabels(bin_labels, rotation=45, ha="right")
        ax1.grid(axis="y", alpha=0.3)

        ax2 = ax1.twinx()
        ax2.plot(sorted_log, cumulative, color="red", linewidth=2.5,
                 label=f"{group} cumulative")
        ax2.set_ylabel("Cumulative unique peptides", color="red")
        ax2.tick_params(axis="y", labelcolor="red")
        ax2.set_ylim(0, 400)
        ax2.set_yticks(np.arange(0, 401, 50))

        plt.title(f"{group}: unique peptide binning and cumulative curve")
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

        plt.tight_layout()
        save_path = os.path.join(output_dir, f"{group}_bin_cumulative_combined.png")
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"Saved: {save_path}")

    fig, ax1 = plt.subplots(figsize=(12, 7))
    colors = sns.color_palette("Set1", n_colors=len(groups))

    all_bin_counts = []
    for group in groups:
        group_data = data[data["group"] == group]
        bins = create_log_bins(group_data["log10_count"], num_bins)
        group_data = group_data.copy()
        group_data["bin"] = pd.cut(group_data["log10_count"], bins=bins)
        all_bin_counts.append(group_data.groupby("bin", observed=False).size())

    bin_means = (
        pd.concat(all_bin_counts, axis=1).mean(axis=1).reset_index(name="mean_count")
    )
    bin_centers = [(b.left + b.right) / 2 for b in bin_means["bin"]]
    ax1.bar(bin_centers, bin_means["mean_count"],
            width=(bins[1] - bins[0]) * 0.8,
            color="gray", edgecolor="black", alpha=0.6,
            label="Mean bin count (all groups)")

    ax2 = ax1.twinx()
    for i, group in enumerate(groups):
        sorted_log = np.sort(data[data["group"] == group]["log10_count"])
        cumulative = np.arange(1, len(sorted_log) + 1)
        ax2.plot(sorted_log, cumulative, color=colors[i], linewidth=2, label=group)

    ax1.set_xlabel("log10(read count)")
    ax1.set_ylabel("Unique peptides in bin (mean)", color="black")
    ax2.set_ylabel("Cumulative unique peptides", color="red")
    ax2.tick_params(axis="y", labelcolor="red")
    ax2.set_ylim(0, 400)
    ax2.set_yticks(np.arange(0, 401, 50))

    plt.title("All groups: mean bin counts and cumulative curves")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    plt.tight_layout()
    save_path = os.path.join(output_dir, "all_groups_bin_cumulative_combined.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved: {save_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, help="directory with BLAST CSVs")
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    output_dir = args.output_dir or os.path.join(
        args.input_dir, "combined_bin_cumulative_plots"
    )
    os.makedirs(output_dir, exist_ok=True)

    unique_data = process_data(args.input_dir)
    if unique_data.empty:
        raise SystemExit("No valid data found; nothing to plot.")

    unique_data.to_csv(os.path.join(output_dir, "group_unique_data.csv"), index=False)
    plot_combined_bin_cumulative(unique_data, output_dir)
    print(f"All outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()