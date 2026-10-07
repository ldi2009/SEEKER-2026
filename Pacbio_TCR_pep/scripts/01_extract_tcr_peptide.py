#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1 (pacbio_tcr_peptide): per-read TCR + peptide extraction from HiFi FASTQ.

For every read (length window configurable, default 1800-2500 bp) the script:
  1. regex-extracts the peptide DNA and the TCRA / TCRB / TCRfull inserts
     (forward orientation first; falls back to reverse complement),
  2. aligns the TCR inserts to the TCR reference with `minimap2 -ax map-hifi`
     and keeps the best hit per insert,
  3. translates the peptide.

Outputs (prefix = output_prefix):
    <prefix>_result.csv       all parsed reads (per-read table)
    <prefix>_result_qc.csv    TCRfull CIGAR all-M and peptide non-empty only
    <prefix>_count_combined_counts.csv            unique-combination counts (all)
    <prefix>_result_qc_count_combined_counts.csv  unique-combination counts (QC)

Usage:
    python 01_extract_tcr_peptide.py <input.fastq> <tcr_reference.fasta> <output_prefix> \
        [--threads 8] [--min_len 1800] [--max_len 2500]
"""

import os
import re
import csv
import sys
import time
import tempfile
import argparse
import subprocess
from multiprocessing import Pool

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# ---- construct anchors (default 2A-based design) ----
FLANK_5PRIME = 'TCGCGCTACTCTCTCTTTCTGGCCTGGAGGCT'
FLANK_3PRIME = 'GGATGCGGAGGGTCCGGCGGGGGAGGTA'

PATTERN_PEPTIDE_POS = re.compile(
    FLANK_5PRIME + r'([ATCG]{0,40}?)' + FLANK_3PRIME)
PATTERN_PEPTIDE_NEG = re.compile(
    'TACCTCCCCCGCCGGACCCTCCGCATCC' + r'([ATCG]{0,40}?)' +
    'AGCCTCCAGGCCAGAAAGAGAGAGTAGCGCGA')
PATTERN_TCR_ALPHA = re.compile(r'AGGAGAACCCCGGCCCC([ATCG]{300,550}?)TCTGTCTGCCTATTCACCGAT')
PATTERN_TCR_BETA = re.compile(r'CGGTGTAGAGCCGCCACC([ATCG]{300,550}?)GAGATCTCCCACACCCAAAAGG')
PATTERN_TCR_FULL = re.compile(r'CGGTGTAGAGCCGCCACC([ATCG]{1300,1600}?)TCTGTCTGCCTATTCACCGAT')

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def reverse_complement(seq):
    return ''.join([COMPLEMENT.get(base, base) for base in reversed(seq)])


def build_minimap2_index(reference_fasta, output_prefix):
    cmd = f"minimap2 -d {output_prefix}.mmi {reference_fasta}"
    subprocess.run(cmd, shell=True, check=True)
    return f"{output_prefix}.mmi"


def align_with_minimap2(query_fasta, index_file, threads=4):
    cmd = f"minimap2 -ax map-hifi -t {threads} {index_file} {query_fasta}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout


def process_alignment_output(alignment_output):
    """Keep the best (highest-identity) hit per query."""
    best_hits = {}
    for line in alignment_output.strip().split('\n'):
        if not line or line.startswith('@'):
            continue
        parts = line.split('\t')
        if len(parts) < 12:
            continue
        query_id = parts[0]
        subject_id = parts[2]
        try:
            cigar = parts[5]
            matches = 0
            for count, op in re.findall(r'(\d+)([M=X])', cigar):
                if op in ['M', '=']:
                    matches += int(count)
            if len(parts) > 10 and parts[10].isdigit():
                alignment_length = int(parts[10])
                identity = (matches / alignment_length) * 100 if alignment_length > 0 else 0
            else:
                identity = 0
        except Exception:
            identity = 0
        if query_id not in best_hits or identity > best_hits[query_id]['identity']:
            best_hits[query_id] = {'subject_id': subject_id, 'identity': identity,
                                   'alignment': line}
    return best_hits


def dump_single_record_seq(record_seq):
    """Extract (tcra, tcrb, tcrfull, peptide_dna) from one read, forward first."""
    forward_count = 0
    tcra_seq = tcrb_seq = tcrfull_seq = peptide_dna = None

    m = PATTERN_PEPTIDE_POS.search(record_seq)
    if m is not None:
        forward_count += 1
        peptide_dna = m[1]
    m = PATTERN_TCR_ALPHA.search(record_seq)
    if m is not None:
        forward_count += 1
        tcra_seq = m[1]
    m = PATTERN_TCR_BETA.search(record_seq)
    if m is not None:
        forward_count += 1
        tcrb_seq = m[1]
    m = PATTERN_TCR_FULL.search(record_seq)
    if m is not None:
        forward_count += 1
        tcrfull_seq = m[1]

    if forward_count >= 3:
        fail = (0 if peptide_dna else 1, 0 if tcra_seq else 1,
                0 if tcrb_seq else 1, 0 if tcrfull_seq else 1)
        if forward_count == 0:
            return None, fail
        return (tcra_seq, tcrb_seq, tcrfull_seq, peptide_dna, fail), fail

    reverse_seq = reverse_complement(record_seq)
    reverse_count = 0
    r_tcra = r_tcrb = r_tcrfull = r_peptide = None

    m = PATTERN_PEPTIDE_NEG.search(record_seq)
    if m is not None:
        reverse_count += 1
        r_peptide = reverse_complement(m[1])
    m = PATTERN_TCR_ALPHA.search(reverse_seq)
    if m is not None:
        reverse_count += 1
        r_tcra = m[1]
    m = PATTERN_TCR_BETA.search(reverse_seq)
    if m is not None:
        reverse_count += 1
        r_tcrb = m[1]
    m = PATTERN_TCR_FULL.search(reverse_seq)
    if m is not None:
        reverse_count += 1
        r_tcrfull = m[1]

    if reverse_count > forward_count:
        tcra_seq, tcrb_seq, tcrfull_seq, peptide_dna = r_tcra, r_tcrb, r_tcrfull, r_peptide
    # tie -> keep forward

    fail = (0 if peptide_dna else 1, 0 if tcra_seq else 1,
            0 if tcrb_seq else 1, 0 if tcrfull_seq else 1)
    total = max(forward_count, reverse_count)
    if total == 0:
        return None, fail
    return (tcra_seq, tcrb_seq, tcrfull_seq, peptide_dna, fail), fail


def align_tcr_sequences(tcra_seq, tcrb_seq, tcrfull_seq, tcr_index, threads):
    with tempfile.TemporaryDirectory() as temp_dir:
        query_fasta = os.path.join(temp_dir, 'query.fasta')
        with open(query_fasta, 'w') as f:
            if tcra_seq:
                f.write(f">tcra\n{tcra_seq}\n")
            if tcrb_seq:
                f.write(f">tcrb\n{tcrb_seq}\n")
            if tcrfull_seq:
                f.write(f">tcrfull\n{tcrfull_seq}\n")
        alignment_output = align_with_minimap2(query_fasta, tcr_index, threads)
        best_hits = process_alignment_output(alignment_output)

    return {
        'tcra': best_hits.get('tcra', {'subject_id': 'NA', 'identity': 0}),
        'tcrb': best_hits.get('tcrb', {'subject_id': 'NA', 'identity': 0}),
        'tcrfull': best_hits.get('tcrfull', {'subject_id': 'NA', 'identity': 0}),
    }


def process_fastq_chunk(args):
    chunk, tcr_index, threads, min_len, max_len = args
    results = []
    combination_counts = {(p, a, b, f): 0
                          for p in (0, 1) for a in (0, 1) for b in (0, 1) for f in (0, 1)}

    for record in chunk:
        seq_len = len(record.seq)
        if not (min_len <= seq_len <= max_len):
            continue

        triplet, fail_combination = dump_single_record_seq(str(record.seq))
        combination_counts[fail_combination] += 1
        if triplet is None:
            continue

        tcra_seq, tcrb_seq, tcrfull_seq, peptide_dna, _ = triplet
        alignments = align_tcr_sequences(tcra_seq, tcrb_seq, tcrfull_seq, tcr_index, threads)
        peptide = str(Seq(peptide_dna).translate()) if peptide_dna else ''

        def cigar_of(align_entry):
            if 'alignment' in align_entry:
                parts = align_entry['alignment'].split('\t')
                if len(parts) > 5:
                    return parts[5]
            return ''

        results.append({
            'TCRA_seq': tcra_seq,
            'TCRA_align': alignments['tcra']['subject_id'],
            'TCRA_identity': alignments['tcra']['identity'],
            'TCRA_cigar': cigar_of(alignments['tcra']),
            'TCRB_seq': tcrb_seq,
            'TCRB_align': alignments['tcrb']['subject_id'],
            'TCRB_identity': alignments['tcrb']['identity'],
            'TCRB_cigar': cigar_of(alignments['tcrb']),
            'TCRfull_seq': tcrfull_seq,
            'TCRfull_align': alignments['tcrfull']['subject_id'],
            'TCRfull_identity': alignments['tcrfull']['identity'],
            'TCRfull_cigar': cigar_of(alignments['tcrfull']),
            'peptideDNA': peptide_dna,
            'peptide': peptide,
        })
    return results, combination_counts


def has_only_M(cigar):
    if not cigar or not isinstance(cigar, str):
        return False
    return bool(re.fullmatch(r'\d+M', cigar))


def generate_combined_counts(df, count_output):
    combined_counts = {}
    for _, row in df.iterrows():
        key = (row.get('TCRA_align', ''), row.get('TCRB_align', ''),
               row.get('TCRfull_align', ''), row.get('peptideDNA', ''),
               row.get('peptide', ''))
        combined_counts[key] = combined_counts.get(key, 0) + 1
    with open(count_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['TCRA_align', 'TCRB_align', 'TCRfull_align',
                         'peptideDNA', 'peptide', 'count'])
        for combo, count in sorted(combined_counts.items(), key=lambda x: x[1], reverse=True):
            writer.writerow([*combo, count])


def main():
    ap = argparse.ArgumentParser(description="Extract TCR + peptide per HiFi read")
    ap.add_argument('input_fastq', help='HiFi FASTQ')
    ap.add_argument('tcr_reference', help='TCR reference FASTA (minimap2 target)')
    ap.add_argument('output_prefix', help='output prefix for the CSV files')
    ap.add_argument('--threads', type=int, default=8, help='minimap2 / worker threads')
    ap.add_argument('--min_len', type=int, default=1800, help='min read length')
    ap.add_argument('--max_len', type=int, default=2500, help='max read length')
    args = ap.parse_args()

    start_time = time.time()
    print(f"Input FASTQ: {args.input_fastq}")
    print(f"TCR reference: {args.tcr_reference}")
    print(f"Output prefix: {args.output_prefix}")

    batch_size = 1600

    print("Building TCR reference index...")
    with tempfile.TemporaryDirectory() as temp_dir:
        tcr_index = build_minimap2_index(args.tcr_reference,
                                         os.path.join(temp_dir, 'tcr_index'))

        print("Reading FASTQ and splitting into batches...")
        chunks, current_chunk, record_count = [], [], 0
        for record in SeqIO.parse(args.input_fastq, 'fastq'):
            current_chunk.append(record)
            record_count += 1
            if len(current_chunk) >= batch_size:
                chunks.append(current_chunk)
                current_chunk = []
                if record_count % (batch_size * 10) == 0:
                    print(f"  read {record_count:,} records, {len(chunks)} batches")
        if current_chunk:
            chunks.append(current_chunk)
        print(f"Total {record_count:,} records in {len(chunks)} batches")

        print("Processing in parallel...")
        pool_args = [(chunk, tcr_index, args.threads, args.min_len, args.max_len)
                     for chunk in chunks]
        with Pool(processes=args.threads) as pool:
            results_list = pool.map(process_fastq_chunk, pool_args)

    results = []
    combination_counts = {}
    for chunk_results, chunk_combos in results_list:
        results.extend(chunk_results)
        for combo, count in chunk_combos.items():
            combination_counts[combo] = combination_counts.get(combo, 0) + count

    print("\nExtraction summary (0=success, 1=fail; peptide, TCRA, TCRB, TCRfull):")
    for combo, count in sorted(combination_counts.items(), key=lambda x: -x[1]):
        if count > 0:
            print(f"  {combo} - {count}")

    result_output = f"{args.output_prefix}_result.csv"
    df = pd.DataFrame(results)
    df.to_csv(result_output, index=False)

    qc_df = df[(df['TCRfull_cigar'].apply(has_only_M)) &
               (df['peptide'].notna()) & (df['peptide'] != '')]
    qc_output = f"{args.output_prefix}_result_qc.csv"
    qc_df.to_csv(qc_output, index=False)

    generate_combined_counts(df, f"{args.output_prefix}_count_combined_counts.csv")
    generate_combined_counts(qc_df, f"{args.output_prefix}_result_qc_count_combined_counts.csv")

    print(f"\nDone! Parsed {len(results):,} reads "
          f"(QC-pass with full-M TCRfull CIGAR + peptide: {len(qc_df):,})")
    print(f"Outputs:\n  {result_output}\n  {qc_output}\n"
          f"  {args.output_prefix}_count_combined_counts.csv\n"
          f"  {args.output_prefix}_result_qc_count_combined_counts.csv")
    print(f"Processing time: {time.time() - start_time:.2f}s")


if __name__ == "__main__":
    main()
