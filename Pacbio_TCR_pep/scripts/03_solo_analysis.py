#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 3 (pacbio_tcr_peptide): single-sample "solo" analysis of two rounds.

Pipeline:
  1. read the two per-round final_QC.csv tables
  2. merge on (TCRA_align, TCRB_align, peptideDNA), no scaling
  3. peptide-origin lookup against the human/virus peptide-DNA libraries
  4. per-round percentages + fold change
  5. per-peptide dedup (keep highest count in the last round) + 1/20 threshold
  6. 7-step clean filter:
       remove control peptides (configurable keywords, e.g. DMF5/TCR_074),
       remove empty TCR chains, remove unannotated peptides,
       keep special peptides or length >= min_peptide_len,
       fold-change window, allowed TCR set
  7. primer design (AGGCT+peptideDNA / ATCC+revcomp+A), numbered 1..N

Outputs in --output_dir:
    <sample>_analysis_result_combined.csv
    <sample>_analysis_result_with_origin.csv
    <sample>_analysis_result_solo_filtered.csv
    filtered_data_clean.csv
    filtered_data_clean_with_primer.csv   <- final result

Usage:
    python 03_solo_analysis.py <roundA_final_QC.csv> <roundB_final_QC.csv> \
        --output_dir analysis --sample_name 2704-V2-R23 \
        --human_proteome PATH --virus_proteome PATH \
        --allowed_tcr_numbers "1-49,101" \
        [--fold_change_min 0.1 --fold_change_max 60] \
        [--min_count_ratio 0.05] [--min_peptide_len 8] \
        [--special_peptides P1 P2 ...] [--remove_keywords K1 K2 ...]
"""

import os
import re
import time
import argparse
from datetime import datetime

import numpy as np
import pandas as pd

FLANK_5PRIME = 'CTGTGCTCGCGCTACTCTCTCTTTCTGGCCTGGAGGCT'
FLANK_3PRIME = 'GGATGCGGAGGGTCCGGCGGGGGA'


def fmt_time(seconds):
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        return f"{seconds/60:.1f}min"
    return f"{seconds/3600:.1f}h"


def parse_fasta(fasta_path):
    sequences = {}
    current_id, current_seq = None, []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq).upper()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        sequences[current_id] = ''.join(current_seq).upper()
    return sequences


def build_peptidedna_map(fasta_file, origin_type):
    """peptideDNA -> {fasta_name, origin, Gene} from a flank-flanked peptide library."""
    print(f"    Reading FASTA: {fasta_file}")
    sequences = parse_fasta(fasta_file)
    print(f"      Sequences: {len(sequences):,}")
    peptide_map = {}
    matched = 0
    for seq_id, seq in sequences.items():
        if FLANK_5PRIME in seq and FLANK_3PRIME in seq:
            start = seq.find(FLANK_5PRIME)
            end = seq.find(FLANK_3PRIME, start + len(FLANK_5PRIME))
            if end > start + len(FLANK_5PRIME):
                peptideDNA = seq[start + len(FLANK_5PRIME):end]
                if peptideDNA not in peptide_map:
                    gene_name = ''
                    parts = seq_id.split('+')
                    if len(parts) >= 2:
                        id_part = parts[1]
                        id_parts = id_part.split('_')
                        if len(id_parts) >= 3:
                            gene_name = id_parts[2].rstrip('_')
                        elif len(id_parts) == 2:
                            gene_name = id_parts[1]
                        else:
                            gene_name = id_part
                    peptide_map[peptideDNA] = {'fasta_name': seq_id,
                                               'origin': origin_type,
                                               'Gene': gene_name}
                    matched += 1
    print(f"      peptideDNA extracted: {matched:,} (unique: {len(peptide_map):,})")
    return peptide_map


def combine_two_rounds(df_a, df_b, label_a, label_b):
    """Merge two rounds on (TCRA_align, TCRB_align, peptideDNA); no scaling."""
    print("  " + "=" * 70)
    print(f"  Two-round merge ({label_a} + {label_b}, no scaling)")
    print("  " + "=" * 70)
    start_time = time.time()
    print(f"    {label_a}: {len(df_a):,} rows")
    print(f"    {label_b}: {len(df_b):,} rows")

    df_a = df_a.copy()
    df_b = df_b.copy()
    df_a['key'] = df_a.apply(lambda x: (x['TCRA_align'], x['TCRB_align'], x['peptideDNA']), axis=1)
    df_b['key'] = df_b.apply(lambda x: (x['TCRA_align'], x['TCRB_align'], x['peptideDNA']), axis=1)

    counts_a = dict(zip(df_a['key'], df_a['count']))
    counts_b = dict(zip(df_b['key'], df_b['count']))

    row_a, row_b = {}, {}
    for _, row in df_a.iterrows():
        row_a.setdefault(row['key'], {'peptide': row['peptide'],
                                      'TCRfull_align': row['TCRfull_align']})
    for _, row in df_b.iterrows():
        row_b.setdefault(row['key'], {'peptide': row['peptide'],
                                      'TCRfull_align': row['TCRfull_align']})

    all_keys = set(counts_a) | set(counts_b)
    print(f"    Total combination keys: {len(all_keys):,}")

    combined = []
    for key in all_keys:
        tcra, tcrb, peptide_dna = key
        # prefer the later round's annotation, fall back to the earlier round
        info = row_b.get(key) or row_a.get(key) or {}
        combined.append({
            'TCRA_align': tcra, 'TCRB_align': tcrb,
            'TCRfull_align': info.get('TCRfull_align', ''),
            'peptideDNA': peptide_dna,
            'peptide': info.get('peptide', ''),
            f'count_{label_a}': counts_a.get(key, 0),
            f'count_{label_b}': counts_b.get(key, 0),
        })

    combined_df = pd.DataFrame(combined).sort_values(f'count_{label_b}', ascending=False)
    print(f"    Combined rows: {len(combined_df):,} | {fmt_time(time.time()-start_time)}")
    return combined_df


def find_peptide_origin(df, human_map, virus_map):
    print("  Peptide-origin lookup...")
    unique_pdna = df['peptideDNA'].unique()
    print(f"    Unique peptideDNA: {len(unique_pdna):,}")
    origin_info = {}
    found_human = found_virus = not_found = 0
    for ptdna in unique_pdna:
        if not ptdna or pd.isna(ptdna):
            origin_info[ptdna] = {'fasta_name': 'INFO=MISSING', 'origin': 'unknown', 'Gene': '-'}
            not_found += 1
            continue
        ptdna_upper = str(ptdna).upper()
        if ptdna_upper in human_map:
            origin_info[ptdna] = human_map[ptdna_upper]
            found_human += 1
        elif ptdna_upper in virus_map:
            origin_info[ptdna] = virus_map[ptdna_upper]
            found_virus += 1
        else:
            origin_info[ptdna] = {'fasta_name': 'INFO=MISSING', 'origin': 'unknown', 'Gene': '-'}
            not_found += 1
    print(f"    Human: {found_human:,} | Virus: {found_virus:,} | Not found: {not_found:,}")
    df['fasta_name'] = df['peptideDNA'].map(
        lambda x: origin_info.get(x, {}).get('fasta_name', 'INFO=MISSING'))
    df['origin'] = df['peptideDNA'].map(
        lambda x: origin_info.get(x, {}).get('origin', 'unknown'))
    df['Gene'] = df['peptideDNA'].map(
        lambda x: origin_info.get(x, {}).get('Gene', '-'))
    return df


def calculate_percent_fold_change(df, round_a_label, round_b_label):
    print("  Percent + fold change...")
    total_a = df[f'count_{round_a_label}'].sum()
    total_b = df[f'count_{round_b_label}'].sum()
    print(f"    {round_a_label} total count: {total_a:,}")
    print(f"    {round_b_label} total count: {total_b:,}")
    df[f'percent_{round_a_label}'] = df[f'count_{round_a_label}'] / total_a * 100
    df[f'percent_{round_b_label}'] = df[f'count_{round_b_label}'] / total_b * 100
    df['fold_change'] = np.where(
        df[f'percent_{round_a_label}'] > 0,
        df[f'percent_{round_b_label}'] / df[f'percent_{round_a_label}'], np.inf)
    return df


def solo_dedup_and_filter(df, last_round_label, min_count_ratio):
    print("  Per-peptide dedup + 1/20 threshold...")
    print(f"    Input rows: {len(df):,}")
    df_sorted = df.sort_values(['peptide', f'count_{last_round_label}'],
                               ascending=[True, False])
    solo_df = df_sorted.groupby('peptide').first().reset_index()
    print(f"    After dedup: {len(solo_df):,}")
    peptide_max = solo_df.groupby('peptide')[f'count_{last_round_label}'].transform('max')
    solo_df['max_count'] = peptide_max
    filtered = solo_df[solo_df[f'count_{last_round_label}'] >=
                       solo_df['max_count'] * min_count_ratio]
    filtered = filtered.drop(columns=['max_count']).sort_values(
        f'count_{last_round_label}', ascending=False)
    print(f"    After threshold: {len(filtered):,}")
    return filtered


def clean_filter(df, args, allowed_tcr_numbers, last_round_label):
    print("  7-step clean filter...")
    initial = len(df)
    print(f"    Input rows: {initial:,}")

    for kw in args.remove_keywords:
        df = df[~df['peptide'].str.contains(kw, na=False)]
        print(f"    Step 1-2: removed '{kw}' -> {len(df):,}")

    df = df[(df['TCRA_align'].notna()) & (df['TCRA_align'] != '') &
            (df['TCRB_align'].notna()) & (df['TCRB_align'] != '')]
    print(f"    Step 3: removed empty TCR chains -> {len(df):,}")

    df = df[(df['Gene'] != '-') & (df['Gene'].notna()) &
            (~df['fasta_name'].str.contains('INFO=MISSING', na=False)) &
            (~df['fasta_name'].str.contains('INFO=UNKNOWN', na=False))]
    print(f"    Step 4: removed unannotated peptides -> {len(df):,}")

    peptide_len = df['peptide'].str.len()
    df = df[df['peptide'].isin(set(args.special_peptides)) | (peptide_len >= args.min_peptide_len)]
    print(f"    Step 5: kept special peptides or len >= {args.min_peptide_len} -> {len(df):,}")

    df = df[(df['fold_change'] >= args.fold_change_min) &
            (df['fold_change'] <= args.fold_change_max)]
    print(f"    Step 6: fold change in [{args.fold_change_min}, {args.fold_change_max}] "
          f"-> {len(df):,}")

    allowed_tcr = {f'TCR_{i:03d}' for i in allowed_tcr_numbers}
    df = df[df['TCRA_align'].isin(allowed_tcr) & df['TCRB_align'].isin(allowed_tcr)]
    print(f"    Step 7: TCR in allowed set -> {len(df):,}")

    df = df.sort_values(f'count_{last_round_label}', ascending=False)
    print(f"    Clean filter done: {initial:,} -> {len(df):,}")
    return df


def design_primers(df):
    print("  Primer design...")
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    def reverse_complement(seq):
        if not seq or pd.isna(seq):
            return ''
        return ''.join([comp.get(b, b) for b in str(seq)[::-1]])

    df = df.copy()
    df['forward_primer'] = 'AGGCT' + df['peptideDNA'].astype(str)
    df['reverse_primer'] = df['peptideDNA'].apply(
        lambda x: 'ATCC' + reverse_complement(x) + 'A')
    df['primer_number'] = range(1, len(df) + 1)
    df['forward_name'] = df['primer_number'].astype(str) + '_F'
    df['reverse_name'] = df['primer_number'].astype(str) + '_R'
    print(f"    Primers added for {len(df):,} rows")
    return df


def parse_tcr_numbers(spec):
    """Parse '1-49,101' -> [1..49, 101]."""
    numbers = []
    for part in str(spec).split(','):
        part = part.strip()
        if '-' in part:
            lo, hi = part.split('-')
            numbers.extend(range(int(lo), int(hi) + 1))
        elif part:
            numbers.append(int(part))
    return numbers


def main():
    ap = argparse.ArgumentParser(description="Solo two-round analysis for the PacBio workflow")
    ap.add_argument('qc_round_a', help='final_QC.csv of the earlier round')
    ap.add_argument('qc_round_b', help='final_QC.csv of the later round')
    ap.add_argument('--output_dir', default='analysis')
    ap.add_argument('--sample_name', default='sample')
    ap.add_argument('--round_a_label', default='round2',
                    help='count column label for round A (e.g. round2 for R2 data)')
    ap.add_argument('--round_b_label', default='round3',
                    help='count column label for round B (e.g. round3 for R3 data)')
    ap.add_argument('--human_proteome', required=True,
                    help='human peptide-DNA library FASTA')
    ap.add_argument('--virus_proteome', required=True,
                    help='virome peptide-DNA library FASTA')
    ap.add_argument('--allowed_tcr_numbers', default='1-49,101',
                    help='e.g. "1-49,101"')
    ap.add_argument('--fold_change_min', type=float, default=0.1)
    ap.add_argument('--fold_change_max', type=float, default=60)
    ap.add_argument('--min_count_ratio', type=float, default=0.05,
                    help='per-peptide count threshold ratio (1/20 default)')
    ap.add_argument('--min_peptide_len', type=int, default=8)
    ap.add_argument('--special_peptides', nargs='*', default=[])
    ap.add_argument('--remove_keywords', nargs='*', default=[])
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    start = time.time()
    print(f"\n{'='*80}")
    print(f"Solo analysis - {args.sample_name}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*80}")

    df_a = pd.read_csv(args.qc_round_a)
    df_b = pd.read_csv(args.qc_round_b)
    print(f"  Round A ({args.qc_round_a}): {len(df_a):,} rows")
    print(f"  Round B ({args.qc_round_b}): {len(df_b):,} rows")

    combined = combine_two_rounds(df_a, df_b, args.round_a_label, args.round_b_label)
    combined_file = os.path.join(args.output_dir,
                                 f'{args.sample_name}_analysis_result_combined.csv')
    combined.to_csv(combined_file, index=False)
    del df_a, df_b

    human_map = build_peptidedna_map(args.human_proteome, 'human')
    virus_map = build_peptidedna_map(args.virus_proteome, 'virus')
    df_origin = find_peptide_origin(combined, human_map, virus_map)
    del combined, human_map, virus_map

    df_percent = calculate_percent_fold_change(df_origin, args.round_a_label, args.round_b_label)
    origin_file = os.path.join(args.output_dir,
                               f'{args.sample_name}_analysis_result_with_origin.csv')
    df_percent.to_csv(origin_file, index=False)
    del df_origin

    solo_filtered = solo_dedup_and_filter(df_percent, args.round_b_label, args.min_count_ratio)
    solo_file = os.path.join(args.output_dir,
                             f'{args.sample_name}_analysis_result_solo_filtered.csv')
    solo_filtered.to_csv(solo_file, index=False)
    del df_percent

    allowed = parse_tcr_numbers(args.allowed_tcr_numbers)
    clean_filtered = clean_filter(solo_filtered, args, allowed, args.round_b_label)
    clean_file = os.path.join(args.output_dir, 'filtered_data_clean.csv')
    clean_filtered.to_csv(clean_file, index=False)
    del solo_filtered

    final_result = design_primers(clean_filtered)
    primer_file = os.path.join(args.output_dir, 'filtered_data_clean_with_primer.csv')
    final_result.to_csv(primer_file, index=False)

    print(f"\n{'='*80}")
    print(f"Solo analysis done! Final rows: {len(final_result):,}")
    print(f"Final result: {primer_file}")
    print(f"Total time: {fmt_time(time.time()-start)}")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
