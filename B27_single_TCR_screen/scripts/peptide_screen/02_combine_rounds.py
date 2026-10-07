#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2 (peptide_screen): merge two enrichment rounds.

Input : two fullresult CSVs (nt, count, aa, length; headerless), one per round.
Output: merged CSV with header, percentages per round and the R2/R1 enrichment
        ratio, sorted by count_round2 then count_round1 (descending).

Usage:
    python 02_combine_rounds.py <r1_fullresult.csv> <r2_fullresult.csv> <output.csv> <sample_name>
"""

import time
import argparse
from datetime import datetime

import pandas as pd


def combine_two_rounds(r1_file, r2_file, output_file, sample_name=""):
    print(f"\n{'='*60}")
    print(f"Sample: {sample_name}")
    print(f"{'='*60}")
    print(f"R1 file: {r1_file}")
    print(f"R2 file: {r2_file}")
    print(f"Output file: {output_file}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    start_time = time.time()

    seq_r1 = pd.read_csv(r1_file, header=None,
                         names=["nt", "count", "aa", "length"]).drop_duplicates()
    seq_r2 = pd.read_csv(r2_file, header=None,
                         names=["nt", "count", "aa", "length"]).drop_duplicates()

    print(f"\nR1 unique sequences: {len(seq_r1):,}")
    print(f"R2 unique sequences: {len(seq_r2):,}")

    total_r1 = seq_r1["count"].sum()
    total_r2 = seq_r2["count"].sum()
    print(f"R1 total reads: {total_r1:,}")
    print(f"R2 total reads: {total_r2:,}")

    seq_r1["count_pct_round1"] = seq_r1["count"] / total_r1 * 100
    seq_r2["count_pct_round2"] = seq_r2["count"] / total_r2 * 100
    seq_r1 = seq_r1.rename(columns={"count": "count_round1"})
    seq_r2 = seq_r2.rename(columns={"count": "count_round2"})

    print("\nMerging data...")
    merged = pd.merge(
        seq_r1[["nt", "aa", "length", "count_round1", "count_pct_round1"]],
        seq_r2[["nt", "aa", "length", "count_round2", "count_pct_round2"]],
        on=["nt", "aa", "length"],
        how="outer"
    ).fillna(0)

    merged["r2_r1_ratio"] = round(
        merged["count_pct_round2"] / merged["count_pct_round1"].replace(0, float('nan')), 3
    ).fillna(0)

    cols = ["nt", "aa", "length",
            "count_pct_round1", "count_pct_round2",
            "count_round1", "count_round2",
            "r2_r1_ratio"]
    merged = merged[cols].sort_values(
        by=["count_round2", "count_round1"], ascending=[False, False]
    ).drop_duplicates()

    merged.to_csv(output_file, index=False)

    elapsed = time.time() - start_time
    print(f"\nMerge complete!")
    print(f"Total sequences: {len(merged):,}")
    print(f"Shared by both rounds: {((merged['count_round1'] > 0) & (merged['count_round2'] > 0)).sum():,}")
    print(f"R1 only: {((merged['count_round1'] > 0) & (merged['count_round2'] == 0)).sum():,}")
    print(f"R2 only: {((merged['count_round1'] == 0) & (merged['count_round2'] > 0)).sum():,}")
    print(f"Elapsed: {elapsed:.2f}s")

    print(f"\nTop 10 sequences (by R2 count):")
    print(merged.head(10)[["aa", "length", "count_round1", "count_round2",
                           "r2_r1_ratio"]].to_string())
    return merged


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Merge two enrichment rounds (R1 + R2)")
    ap.add_argument('r1_file', help='R1 fullresult CSV (headerless: nt, count, aa, length)')
    ap.add_argument('r2_file', help='R2 fullresult CSV (headerless: nt, count, aa, length)')
    ap.add_argument('output_file', help='merged output CSV')
    ap.add_argument('sample_name', nargs='?', default='', help='sample name for display')
    args = ap.parse_args()
    combine_two_rounds(args.r1_file, args.r2_file, args.output_file, args.sample_name)
