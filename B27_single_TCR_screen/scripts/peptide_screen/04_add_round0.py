#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 4 (peptide_screen): merge pre-selection (R0) peptide counts into the
annotated result.

R0 inputs are CSVs with columns `Peptide` and `Read Count` (human + virome
tables). Adds `count_round0`, `count_pct_round0`, `r1_r0_ratio` and reorders
columns to the canonical layout.

Usage:
    python 04_add_round0.py <input.csv> <output.csv> --r0_virus CSV --r0_human CSV
"""

import time
import argparse
from datetime import datetime

import pandas as pd


def load_r0_data(r0_virus_file, r0_human_file):
    print("Loading R0 (pre-selection) data...")
    t0 = time.time()

    virus_df = pd.read_csv(r0_virus_file)
    human_df = pd.read_csv(r0_human_file)

    r0_df = pd.concat([virus_df, human_df], ignore_index=True)
    r0_df = r0_df.rename(columns={'Peptide': 'aa', 'Read Count': 'count_round0'})

    total_r0 = r0_df['count_round0'].sum()
    r0_df['count_pct_round0'] = r0_df['count_round0'] / total_r0 * 100

    # same peptide may appear in both libraries -> sum counts and percentages
    r0_df = r0_df.groupby('aa').agg(
        count_round0=('count_round0', 'sum'),
        count_pct_round0=('count_pct_round0', 'sum')
    ).reset_index()

    print(f"  Virus table: {len(virus_df):,} rows")
    print(f"  Human table: {len(human_df):,} rows")
    print(f"  After merge by peptide: {len(r0_df):,} rows")
    print(f"  R0 total reads: {total_r0:,}")
    print(f"  Load time: {time.time()-t0:.1f}s")
    return r0_df


def add_r0_to_output(input_csv, output_csv, r0_df, sample_name=""):
    print(f"\nProcessing {sample_name}: {input_csv}")
    t0 = time.time()

    df = pd.read_csv(input_csv)
    print(f"  Input records: {len(df):,}")

    merged = pd.merge(df, r0_df[['aa', 'count_round0', 'count_pct_round0']],
                      on='aa', how='left')
    merged['count_round0'] = merged['count_round0'].fillna(0).astype(int)
    merged['count_pct_round0'] = merged['count_pct_round0'].fillna(0)

    if 'origin' in merged.columns:
        merged['origin'] = merged['origin'].fillna('INFO=UNKNOWN')
    if 'Gene' in merged.columns:
        merged['Gene'] = merged['Gene'].fillna('INFO=UNKNOWN')

    if 'count_round1' in merged.columns and 'count_pct_round1' in merged.columns:
        merged['r1_r0_ratio'] = round(
            merged['count_pct_round1'] / merged['count_pct_round0'].replace(0, float('nan')), 3
        ).fillna(0)

    r0_match = int((merged['count_round0'] > 0).sum())

    # canonical column order
    cols = ['aa', 'count_pct_round0', 'count_round0']
    if 'r1_r0_ratio' in merged.columns:
        cols.append('r1_r0_ratio')
    if 'count_pct_round1' in merged.columns:
        cols.extend(['count_pct_round1', 'count_round1'])
    if 'count_pct_round2' in merged.columns:
        cols.extend(['count_pct_round2', 'count_round2'])
    if 'r2_r1_ratio' in merged.columns:
        cols.append('r2_r1_ratio')
    if 'nt' in merged.columns:
        cols.insert(0, 'nt')
    if 'length' in merged.columns:
        cols.insert(2, 'length')
    if 'fasta_name' in merged.columns:
        cols.extend(['fasta_name', 'origin', 'Gene'])

    merged = merged[[c for c in cols if c in merged.columns]]
    merged.to_csv(output_csv, index=False)

    print(f"  R0 matched:     {r0_match:,} ({r0_match/len(df)*100:.1f}%)")
    print(f"  R0 not matched: {len(df)-r0_match:,} ({(len(df)-r0_match)/len(df)*100:.1f}%)")
    print(f"  Output records: {len(merged):,}")
    print(f"  Output: {output_csv}")
    print(f"  Elapsed: {time.time()-t0:.2f}s")
    return merged


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Add R0 (pre-selection) counts to results")
    ap.add_argument('input_csv', help='input CSV (combined/annotated)')
    ap.add_argument('output_csv', help='output CSV with R0 columns')
    ap.add_argument('--r0_virus', required=True, help='R0 virus peptide counts CSV (Peptide, Read Count)')
    ap.add_argument('--r0_human', required=True, help='R0 human peptide counts CSV (Peptide, Read Count)')
    args = ap.parse_args()

    print(f"\n{'='*60}")
    print("Add R0 (pre-selection) counts")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}")

    start_time = time.time()
    r0_df = load_r0_data(args.r0_virus, args.r0_human)
    add_r0_to_output(args.input_csv, args.output_csv, r0_df)

    print(f"\n{'='*60}")
    print(f"All done! Total time: {time.time()-start_time:.2f}s")
    print(f"{'='*60}")
