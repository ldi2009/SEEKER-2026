#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 5 (peptide_screen): keep only peptides whose origin was resolved
(origin in {human, virus}); drops INFO=MISSING / NA rows.

Usage:
    python 05_filter_matched.py <input.csv> <output.csv> [--origins human virus]
"""

import time
import argparse
from datetime import datetime

import pandas as pd


def filter_matched(input_csv, output_csv, origins=('human', 'virus')):
    print(f"\n{'='*60}")
    print(f"Filter matched peptides (keep origin in {', '.join(origins)})")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}")
    t0 = time.time()

    df = pd.read_csv(input_csv)
    before = len(df)
    print(f"Before filter: {before:,} rows")

    filtered = df[df['origin'].isin(origins)].copy()
    after = len(filtered)
    print(f"After filter:  {after:,} rows (removed {before-after:,} unmatched rows)")

    filtered.to_csv(output_csv, index=False)

    print(f"Output: {output_csv}")
    print(f"Elapsed: {time.time()-t0:.2f}s")
    return filtered


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Keep only rows with resolved origin")
    ap.add_argument('input_csv')
    ap.add_argument('output_csv')
    ap.add_argument('--origins', nargs='+', default=['human', 'virus'],
                    help='origin values to keep (default: human virus)')
    args = ap.parse_args()
    filter_matched(args.input_csv, args.output_csv, tuple(args.origins))
