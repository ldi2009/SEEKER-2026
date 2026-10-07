#!/usr/bin/env python3
"""Unique-peptide calling and per-sample unique-peptide tables (idaOEPCR).

For each BLAST-annotated sample table the 1000 most abundant peptides are kept
and scanned in descending abundance order. A peptide is retained as **unique**
when it lies at a Hamming distance > 2 from every peptide already retained,
i.e. it has no more abundant neighbour within two nucleotides. All comparisons are
made on the nucleotide sequence in the ``peptideDNA`` column; sequences of
different length are never neighbours.

Outputs (written to ``--output-dir``):

* ``hamming_summary_results.csv`` — per sample: ``unique_peptide_count``;
* ``per_sample/<sample>_unique_peptides.csv`` — the retained unique peptides
  of each sample (``sample``, ``peptideDNA``, ``best_alpha_hit``, ``Count``).
"""

import argparse
import os

import pandas as pd

SEQ_COLUMN = "peptideDNA"
TOP_N = 1000


def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return float("inf")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def unique_peptide_rows(df):
    """Greedy unique-peptide scan over the top-``TOP_N`` peptides by Count."""
    top = df.nlargest(TOP_N, "Count").reset_index(drop=True)
    kept_seqs = []
    kept_rows = []
    for _, row in top.iterrows():
        seq = str(row[SEQ_COLUMN])
        if all(hamming_distance(seq, kept) > 2 for kept in kept_seqs):
            kept_seqs.append(seq)
            kept_rows.append(row)
    return pd.DataFrame(kept_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True,
                        help="directory with *_top_10000_blast_sequences.csv")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    per_sample_dir = os.path.join(args.output_dir, "per_sample")
    os.makedirs(per_sample_dir, exist_ok=True)

    csv_files = sorted(f for f in os.listdir(args.input_dir)
                       if "blast" in f and f.endswith(".csv"))
    if not csv_files:
        raise SystemExit("No matching CSV files found.")

    samples = {}
    for filename in csv_files:
        sample = filename.split("_top_10000")[0]
        df = pd.read_csv(os.path.join(args.input_dir, filename))
        unique = unique_peptide_rows(df)
        samples[sample] = {"df": df, "unique": unique}
        print(f"{sample}: {len(unique)} unique peptides")

    summary = []
    for sample, data in samples.items():
        unique = data["unique"]

        table = pd.DataFrame({
            "sample": sample,
            SEQ_COLUMN: unique[SEQ_COLUMN].astype(str).values,
            "best_alpha_hit": unique["best_alpha_hit"].values,
            "Count": unique["Count"].values,
        })
        table.to_csv(os.path.join(per_sample_dir,
                                  f"{sample}_unique_peptides.csv"), index=False)

        summary.append({
            "file": sample,
            "unique_peptide_count": len(unique),
        })

    summary_path = os.path.join(args.output_dir, "hamming_summary_results.csv")
    pd.DataFrame(summary).to_csv(summary_path, index=False)
    print(f"Saved {summary_path}")
    print(f"Saved per-sample tables to {per_sample_dir}")


if __name__ == "__main__":
    main()