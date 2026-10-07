#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2 (pacbio_tcr_peptide): per-round Final QC on a result.csv produced by
01_extract_tcr_peptide.py.

Steps:
  1. read result.csv
  2. raw unique-combination counts (TCRA/TCRB/TCRfull/peptideDNA/peptide)
  3. drop rows with empty peptide
  4. QC filter (TCRfull CIGAR all-M + peptide non-empty)
  5. QC combination counts
  6. TCR mismatch / missing-chain detection
  7. per-peptide max-count dictionary from mismatch+missing rows, then remove
     questionable rows whose count < max_count / 10

Outputs in --output_dir:
    raw_combined_counts.csv, raw_peptide_filtered.csv, result_qc.csv,
    result_qc_count_combined_counts.csv, tcr_mismatch.csv, tcr_missing.csv,
    peptide_max_count_dict.csv, final_QC.csv, questionable.csv

Usage:
    python 02_final_qc.py <result.csv> <output_dir> <sample_name>
"""

import os
import re
import time
import argparse
from datetime import datetime

import pandas as pd


def has_only_M(cigar):
    if not cigar or not isinstance(cigar, str):
        return False
    return bool(re.fullmatch(r'\d+M', cigar))


def fmt_time(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        return f"{seconds/60:.1f}min"
    return f"{seconds/3600:.1f}h"


def combination_counts(df, progress_label=""):
    counts = {}
    total = len(df)
    t0 = time.time()
    for i, (_, row) in enumerate(df.iterrows()):
        if i % 1000000 == 0 and i > 0:
            elapsed = time.time() - t0
            speed = i / elapsed if elapsed > 0 else 0
            eta = (total - i) / speed if speed > 0 else 0
            print(f"    {progress_label} progress: {i:,}/{total:,} "
                  f"({i/total*100:.1f}%) | {speed:,.0f} rows/s | ETA {fmt_time(eta)}")
        tcra = row.get('TCRA_align', '') if pd.notna(row.get('TCRA_align', '')) else ''
        tcrb = row.get('TCRB_align', '') if pd.notna(row.get('TCRB_align', '')) else ''
        tcrfull = row.get('TCRfull_align', '') if pd.notna(row.get('TCRfull_align', '')) else ''
        pdna = row.get('peptideDNA', '') if pd.notna(row.get('peptideDNA', '')) else ''
        pep = row.get('peptide', '') if pd.notna(row.get('peptide', '')) else ''
        key = (tcra, tcrb, tcrfull, pdna, pep)
        counts[key] = counts.get(key, 0) + 1

    out = pd.DataFrame([
        {'TCRA_align': k[0], 'TCRB_align': k[1], 'TCRfull_align': k[2],
         'peptideDNA': k[3], 'peptide': k[4], 'count': v}
        for k, v in counts.items()
    ]).sort_values('count', ascending=False)
    return out


def run_final_qc(result_file, output_dir, sample_name):
    print(f"\n{'='*80}")
    print(f"Final QC - {sample_name}")
    print(f"{'='*80}")
    sample_start = time.time()

    step_start = time.time()
    print(f"  [2.1] Reading result.csv...")
    df = pd.read_csv(result_file, low_memory=False)
    total_rows = len(df)
    print(f"    Records: {total_rows:,} | "
          f"File size: {os.path.getsize(result_file)/1024/1024/1024:.2f} GB | "
          f"{fmt_time(time.time()-step_start)}")

    step_start = time.time()
    print(f"  [2.2] Raw combination counts...")
    raw_counts_df = combination_counts(df, "raw")
    raw_counts_df.to_csv(os.path.join(output_dir, 'raw_combined_counts.csv'), index=False)
    print(f"    Combinations: {len(raw_counts_df):,}, total count: "
          f"{raw_counts_df['count'].sum():,} | {fmt_time(time.time()-step_start)}")

    step_start = time.time()
    print(f"  [2.3] Dropping empty peptides...")
    filtered_raw_df = raw_counts_df[(raw_counts_df['peptide'].notna()) &
                                    (raw_counts_df['peptide'] != '')]
    filtered_raw_df.to_csv(os.path.join(output_dir, 'raw_peptide_filtered.csv'), index=False)
    print(f"    {len(raw_counts_df):,} -> {len(filtered_raw_df):,} | "
          f"{fmt_time(time.time()-step_start)}")

    step_start = time.time()
    print(f"  [2.4] QC filter (TCRfull CIGAR all-M + peptide non-empty)...")
    qc_df = df[(df['TCRfull_cigar'].apply(has_only_M)) &
               (df['peptide'].notna()) & (df['peptide'] != '')]
    qc_df.to_csv(os.path.join(output_dir, 'result_qc.csv'), index=False)
    print(f"    {len(df):,} -> {len(qc_df):,} | {fmt_time(time.time()-step_start)}")
    del df

    step_start = time.time()
    print(f"  [2.5] QC combination counts...")
    qc_counts_df = combination_counts(qc_df, "qc")
    qc_counts_df.to_csv(os.path.join(output_dir, 'result_qc_count_combined_counts.csv'),
                        index=False)
    print(f"    QC combinations: {len(qc_counts_df):,}, total count: "
          f"{qc_counts_df['count'].sum():,} | {fmt_time(time.time()-step_start)}")
    del qc_df

    step_start = time.time()
    print(f"  [2.6] TCR mismatch / missing-chain detection...")
    all_present = ((filtered_raw_df['TCRA_align'].notna()) & (filtered_raw_df['TCRA_align'] != '') &
                   (filtered_raw_df['TCRB_align'].notna()) & (filtered_raw_df['TCRB_align'] != '') &
                   (filtered_raw_df['TCRfull_align'].notna()) & (filtered_raw_df['TCRfull_align'] != ''))
    mismatch_mask = all_present & (
        (filtered_raw_df['TCRA_align'] != filtered_raw_df['TCRB_align']) |
        (filtered_raw_df['TCRA_align'] != filtered_raw_df['TCRfull_align']) |
        (filtered_raw_df['TCRB_align'] != filtered_raw_df['TCRfull_align']))
    mismatch_df = filtered_raw_df[mismatch_mask].sort_values('count', ascending=False)
    mismatch_df.to_csv(os.path.join(output_dir, 'tcr_mismatch.csv'), index=False)
    missing_df = filtered_raw_df[~all_present].sort_values('count', ascending=False)
    missing_df.to_csv(os.path.join(output_dir, 'tcr_missing.csv'), index=False)
    print(f"    Mismatch: {len(mismatch_df):,}, missing: {len(missing_df):,} | "
          f"{fmt_time(time.time()-step_start)}")

    step_start = time.time()
    print(f"  [2.7] Final QC filter (1/10 of per-peptide max mismatch/missing count)...")
    combined_mm = pd.concat([mismatch_df, missing_df], ignore_index=True)
    if len(combined_mm) > 0:
        peptide_max = combined_mm.groupby('peptide')['count'].max().reset_index()
        peptide_max.columns = ['peptide', 'max_count']
    else:
        peptide_max = pd.DataFrame(columns=['peptide', 'max_count'])
    peptide_max.to_csv(os.path.join(output_dir, 'peptide_max_count_dict.csv'), index=False)

    max_dict = dict(zip(peptide_max['peptide'], peptide_max['max_count']))
    qc_counts_df['peptide_max_count'] = qc_counts_df['peptide'].map(max_dict)
    qc_counts_df['threshold'] = qc_counts_df['peptide_max_count'] * 0.1
    qc_counts_df['is_questionable'] = (qc_counts_df['count'] < qc_counts_df['threshold']).fillna(False)

    final_qc_df = qc_counts_df[~qc_counts_df['is_questionable']].drop(
        columns=['peptide_max_count', 'threshold', 'is_questionable'])
    questionable_df = qc_counts_df[qc_counts_df['is_questionable']].drop(
        columns=['peptide_max_count', 'threshold', 'is_questionable'])
    final_qc_df = final_qc_df.sort_values('count', ascending=False)
    questionable_df = questionable_df.sort_values('count', ascending=False)

    final_qc_file = os.path.join(output_dir, 'final_QC.csv')
    final_qc_df.to_csv(final_qc_file, index=False)
    questionable_df.to_csv(os.path.join(output_dir, 'questionable.csv'), index=False)
    print(f"    Final QC: {len(final_qc_df):,} rows "
          f"(total count {final_qc_df['count'].sum():,}), "
          f"questionable: {len(questionable_df):,} | {fmt_time(time.time()-step_start)}")

    print(f"\n  {sample_name} Final QC done! Total: {fmt_time(time.time()-sample_start)}")
    return final_qc_file


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Per-round Final QC for the PacBio workflow")
    ap.add_argument('result_csv', help='result.csv from 01_extract_tcr_peptide.py')
    ap.add_argument('output_dir', help='output directory for QC tables')
    ap.add_argument('sample_name', help='round/sample name')
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    run_final_qc(args.result_csv, args.output_dir, args.sample_name)
