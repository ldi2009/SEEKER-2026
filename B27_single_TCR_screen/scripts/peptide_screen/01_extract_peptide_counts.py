#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1 (peptide_screen): extract peptide-DNA counts from clean FASTQ.

Scans each read for the anchor motif `CTGGAGGCT[peptide DNA]GGATGC`, counts exact
k-mers and writes:
    <outfile_name>          nt, count
    fullresult_<outfile>    nt, count, aa (translated), length
    9mer_<outfile>          27-nt inserts only
    10mer_<outfile>         30-nt inserts only

Files larger than 200M reads are processed in batches; chunks are counted in
parallel (up to 8 workers).

Usage:
    python 01_extract_peptide_counts.py <fastq> <outfile_name> <output_dir> [--max_reads N]
"""

import re
import csv
import mmap
import time
import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

from Bio.Seq import Seq

# peptide DNA sits between the two constant flanks (5' anchor ... 3' anchor)
PATTERN = re.compile(r'CTGGAGGCT([ATCG]*)GGATGC')

BATCH_SIZE = 200_000_000   # reads per batch for very large files
CHUNK_LINES = 2_000_000    # lines per chunk handed to a worker
MAX_WORKERS = 8


def count_total_reads(file_path):
    total_lines = 0
    with open(file_path, 'r') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        while mm.readline():
            total_lines += 1
        mm.close()
    return total_lines // 4  # FASTQ: 4 lines per read


def process_chunk(chunk):
    result = {}
    for line in chunk:
        match = PATTERN.findall(line)
        if match:
            key = match[0]
            result[key] = result.get(key, 0) + 1
    return result


def merge_dicts(dicts):
    merged = {}
    for d in dicts:
        for k, v in d.items():
            merged[k] = merged.get(k, 0) + v
    return merged


def read_file_in_chunks(file_path, chunk_size=CHUNK_LINES, max_reads=None, start_read=0):
    """Yield lists of lines using memory mapping for speed."""
    chunks = []
    total_lines = 0
    start_line = start_read * 4
    max_lines = start_line + (max_reads * 4) if max_reads else None

    with open(file_path, 'r') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        current_line = 0
        while current_line < start_line:
            if not mm.readline():
                break
            current_line += 1

        chunk = []
        line = mm.readline()
        while line:
            chunk.append(line.decode('utf-8'))
            total_lines += 1
            current_line += 1
            if max_lines and current_line >= max_lines:
                break
            if len(chunk) >= chunk_size:
                chunks.append(chunk)
                chunk = []
            line = mm.readline()
        if chunk:
            chunks.append(chunk)
        mm.close()

    if max_reads:
        print(f"Processed {total_lines // 4} reads ({total_lines} lines) from chunk")
    return chunks


def process_batch(input_file, batch_size, start_read):
    print(f"Processing batch starting at read {start_read} with batch size {batch_size}")
    chunks = read_file_in_chunks(input_file, CHUNK_LINES, batch_size, start_read)
    if not chunks:
        return {}

    workers = min(multiprocessing.cpu_count(), MAX_WORKERS)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(process_chunk, chunks))

    batch_result = merge_dicts(results)
    print(f"Batch completed: found {len(batch_result)} unique sequences")
    return batch_result


def main():
    ap = argparse.ArgumentParser(description="Extract peptide-DNA counts from clean FASTQ")
    ap.add_argument('fastq', help='input FASTQ (.clean.fq / decompressed)')
    ap.add_argument('outfile_name', help='output file base name, e.g. TCR59-R2_plasmid_library.csv')
    ap.add_argument('output_dir', help='output directory')
    ap.add_argument('--max_reads', type=int, default=None,
                    help='optional maximum number of reads to process')
    args = ap.parse_args()

    import os
    os.makedirs(args.output_dir, exist_ok=True)

    start = time.time()
    print(f"Processing file: {args.fastq}")
    print(f"Output directory: {args.output_dir}")

    print("Counting total reads...")
    t0 = time.time()
    total_reads = count_total_reads(args.fastq)
    print(f"Total reads in file: {total_reads:,} (counting took {time.time()-t0:.2f}s)")

    effective_total = min(total_reads, args.max_reads) if args.max_reads else total_reads
    print(f"Effective reads to process: {effective_total:,}")

    all_results = []
    if effective_total <= BATCH_SIZE:
        print("Processing as single batch")
        chunks = read_file_in_chunks(args.fastq, CHUNK_LINES, args.max_reads)
        print(f"Split into {len(chunks)} chunks")
        workers = min(multiprocessing.cpu_count(), MAX_WORKERS)
        print(f"Using {workers} workers for parallel processing")
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(process_chunk, chunks))
        all_results.extend(results)
    else:
        num_batches = (effective_total + BATCH_SIZE - 1) // BATCH_SIZE
        print(f"Processing in {num_batches} batches of {BATCH_SIZE:,} reads")
        for i in range(num_batches):
            start_read = i * BATCH_SIZE
            current_batch = min(BATCH_SIZE, effective_total - start_read)
            print(f"\n=== Batch {i+1}/{num_batches} ===")
            print(f"Processing reads {start_read:,} to {start_read + current_batch - 1:,}")
            batch_result = process_batch(args.fastq, current_batch, start_read)
            if batch_result:
                all_results.append(batch_result)

    result_dict = merge_dicts(all_results)
    print(f"\nTotal unique sequences found: {len(result_dict):,}")

    sorted_result = sorted(result_dict.items(), key=lambda s: s[1], reverse=True)

    base_output = os.path.join(args.output_dir, args.outfile_name)
    fullresult_output = os.path.join(args.output_dir, f"fullresult_{args.outfile_name}")
    mer9_output = os.path.join(args.output_dir, f"9mer_{args.outfile_name}")
    mer10_output = os.path.join(args.output_dir, f"10mer_{args.outfile_name}")

    print("Writing output files...")
    buffer_size = 8192 * 10
    sublist9, sublist10 = [], []

    with open(base_output, 'w', newline='', buffering=buffer_size) as f_base, \
         open(fullresult_output, 'w', newline='', buffering=buffer_size) as f_full:
        writer_base = csv.writer(f_base)
        writer_full = csv.writer(f_full)
        for dna_sequence, count in sorted_result:
            writer_base.writerow([dna_sequence, count])
            dna_length = len(dna_sequence)
            if dna_length % 3 == 0:
                protein = str(Seq.translate(dna_sequence))
            else:
                protein = str(Seq.translate(dna_sequence[:-(dna_length % 3)]))
            if dna_length == 27:
                sublist9.append([dna_sequence, count, protein, dna_length])
            if dna_length == 30:
                sublist10.append([dna_sequence, count, protein, dna_length])
            writer_full.writerow([dna_sequence, count, protein, dna_length])

    with open(mer9_output, 'w', newline='', buffering=buffer_size) as f9:
        csv.writer(f9).writerows(sublist9)
    with open(mer10_output, 'w', newline='', buffering=buffer_size) as f10:
        csv.writer(f10).writerows(sublist10)

    print("Processing completed successfully!")
    print(f"- {base_output}")
    print(f"- {fullresult_output}")
    print(f"- {mer9_output}")
    print(f"- {mer10_output}")
    print(f"Total time: {time.time() - start:.2f}s")


if __name__ == "__main__":
    main()
