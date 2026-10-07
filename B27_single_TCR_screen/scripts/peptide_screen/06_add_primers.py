#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 6 (peptide_screen): sort candidates by R2 count (descending, stable) and
design synthesis oligonucleotides.

Adds:
    peptide_number   rank after sorting (1 = highest count_round2)
    primer_F_name    "<n>_F"
    primer-F         AGGCT + peptide DNA
    primer_R_name    "<n>_R"
    primer-R         ATCC + reverse-complement(peptide DNA) + A

Usage:
    python 06_add_primers.py <input.csv> <output.csv> [--no_sort]
"""

import time
import argparse
from datetime import datetime

import pandas as pd

COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def reverse_complement(sequence):
    return ''.join(COMPLEMENT.get(base, base) for base in reversed(str(sequence)))


def add_primer_columns(input_file, output_file, sort_by_count=True, sample_name=""):
    print(f"\n{'='*60}")
    print(f"Add primer sequences - {sample_name}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}")
    t0 = time.time()

    df = pd.read_csv(input_file, low_memory=False)
    print(f"Input records: {len(df):,}")

    if sort_by_count and 'count_round2' in df.columns:
        print("Sorting by count_round2 (descending, stable)...")
        df['_sort_key'] = df['count_round2'].fillna(0)
        df = df.sort_values('_sort_key', ascending=False, kind='mergesort')
        df = df.drop(columns=['_sort_key']).reset_index(drop=True)

    # number peptides after sorting: 1 = highest R2 count
    df['peptide_number'] = range(1, len(df) + 1)

    df['primer_F_name'] = df['peptide_number'].astype(str) + '_F'
    df['primer-F'] = 'AGGCT' + df['nt'].astype(str)

    df['primer_R_name'] = df['peptide_number'].astype(str) + '_R'
    df['primer-R'] = 'ATCC' + df['nt'].apply(reverse_complement) + 'A'

    df.to_csv(output_file, index=False)

    print(f"Output records: {len(df):,}")
    print(f"Output: {output_file}")
    print(f"Elapsed: {time.time()-t0:.2f}s")

    print("\nTop 5 examples:")
    cols = [c for c in ['aa', 'nt', 'primer_F_name', 'primer-F',
                        'primer_R_name', 'primer-R'] if c in df.columns]
    print(df[cols].head().to_string())
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Sort candidates and add primer sequences")
    ap.add_argument('input_csv')
    ap.add_argument('output_csv')
    ap.add_argument('sample_name', nargs='?', default='', help='sample name for display')
    ap.add_argument('--no_sort', action='store_true',
                    help='do not re-sort by count_round2 (keep input order)')
    args = ap.parse_args()
    add_primer_columns(args.input_csv, args.output_csv,
                       sort_by_count=not args.no_sort, sample_name=args.sample_name)
