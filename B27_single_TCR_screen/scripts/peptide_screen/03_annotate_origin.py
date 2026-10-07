#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 3 (peptide_screen): annotate peptide DNA origin by matching against the
human and virome peptide-DNA libraries (exact index match), then resolve gene
names from the UniProt proteome FASTAs.

Adds columns: fasta_name, origin (human/virus/NA), Gene.

Optimized version: pre-builds a peptideDNA -> (fasta_name, origin) index and a
UniProt ID -> gene name index, then annotates the input in one vectorized pass.

Usage:
    python 03_annotate_origin.py <input.csv> <output.csv> <sample_name> \
        [--human_fasta PATH] [--virus_fasta PATH] \
        [--proteome_dir PATH] [--proteome_files F1 F2] \
        [--flank5 SEQ] [--flank3 SEQ]
"""

import os
import time
import argparse
from datetime import datetime

import pandas as pd
from Bio import SeqIO

DEFAULT_FLANK_5 = 'CTGTGCTCGCGCTACTCTCTCTTTCTGGCCTGGAGGCT'
DEFAULT_FLANK_3 = 'GGATGCGGAGGGTCCGGCGGGGGA'


def build_peptide_index(fasta_path, origin_label, flank5, flank3):
    """peptideDNA (sequence between the two flanks) -> (fasta_id, origin)."""
    index = {}
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq).upper()
            if seq.startswith(flank5) and seq.endswith(flank3):
                peptide_seq = seq[len(flank5):len(seq) - len(flank3)]
                index[peptide_seq] = (record.id, origin_label)
    except Exception as e:
        print(f"Error reading FASTA {fasta_path}: {e}")
    return index


def build_gene_index(proteome_files):
    """UniProt ID -> gene name (GN= field) from UniProt proteome FASTAs."""
    gene_index = {}
    total_records = 0
    for fasta_file in proteome_files:
        if not os.path.exists(fasta_file):
            continue
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                total_records += 1
                header = record.description
                parts = header.split('|')
                if len(parts) >= 2:
                    uniprot_id = parts[1]
                    if "GN=" in header:
                        gene_index[uniprot_id] = header.split("GN=")[-1].split()[0]
                    elif len(parts) >= 3:
                        gene_index[uniprot_id] = parts[2].split()[0]
                else:
                    seq_id = record.id
                    if "_" in seq_id:
                        uniprot_id = seq_id.split('_')[1]
                        if "GN=" in header:
                            gene_index[uniprot_id] = header.split("GN=")[-1].split()[0]
        except Exception as e:
            print(f"Error reading FASTA {fasta_file}: {e}")
    return gene_index, total_records


def align_to_reference(input_csv, output_csv, sample_name="",
                       human_fasta=None, virus_fasta=None,
                       proteome_dir=None, proteome_files=None,
                       flank5=DEFAULT_FLANK_5, flank3=DEFAULT_FLANK_3):
    print(f"\n{'='*60}")
    print(f"Reference annotation - {sample_name}")
    print(f"{'='*60}")
    print(f"Input: {input_csv}")
    print(f"Output: {output_csv}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    start_time = time.time()

    print("\nStep 1: building peptide-library index...")
    t0 = time.time()
    human_index = build_peptide_index(human_fasta, 'human', flank5, flank3)
    print(f"  Human peptide index: {len(human_index):,} entries")
    virus_index = build_peptide_index(virus_fasta, 'virus', flank5, flank3)
    print(f"  Virus peptide index: {len(virus_index):,} entries")
    print(f"  Index build time: {time.time()-t0:.1f}s")

    print("\nStep 2: building gene-name index (GN= from UniProt proteomes)...")
    t0 = time.time()
    proteome_paths = [os.path.join(proteome_dir, f)
                      for f in (proteome_files or [])
                      if os.path.exists(os.path.join(proteome_dir, f))]
    gene_index, total_records = build_gene_index(proteome_paths)
    print(f"  Records processed: {total_records:,}")
    print(f"  Gene index entries: {len(gene_index):,}")
    print(f"  Index build time: {time.time()-t0:.1f}s")

    df = pd.read_csv(input_csv)
    total_rows = len(df)
    print(f"\nStep 3: read input - {total_rows:,} records")

    def get_annotation(nt):
        nt_upper = str(nt).upper()
        if nt_upper in human_index:
            matched_name, matched_origin = human_index[nt_upper]
        elif nt_upper in virus_index:
            matched_name, matched_origin = virus_index[nt_upper]
        else:
            return "INFO=MISSING", "NA", "INFO=UNKNOWN"
        try:
            uniprot_id = matched_name.split('_')[1]
        except IndexError:
            uniprot_id = matched_name
        gene_name = gene_index.get(uniprot_id, "INFO=UNKNOWN")
        return matched_name, matched_origin, gene_name

    print("\nStep 4: annotating...")
    results = df['nt'].apply(get_annotation)
    df['fasta_name'] = results.apply(lambda x: x[0])
    df['origin'] = results.apply(lambda x: x[1])
    df['Gene'] = results.apply(lambda x: x[2])

    matched_human = int((df['origin'] == 'human').sum())
    matched_virus = int((df['origin'] == 'virus').sum())
    unmatched = int((df['origin'] == 'NA').sum())

    df.to_csv(output_csv, index=False)

    elapsed = time.time() - start_time
    elapsed_str = (f"{elapsed/3600:.2f}h" if elapsed > 3600 else
                   f"{elapsed/60:.2f}min" if elapsed > 60 else f"{elapsed:.2f}s")

    print(f"\n{'='*60}")
    print(f"Annotation complete!")
    print(f"Human matches: {matched_human:,} ({matched_human/total_rows*100:.1f}%)")
    print(f"Virus matches: {matched_virus:,} ({matched_virus/total_rows*100:.1f}%)")
    print(f"Unmatched:     {unmatched:,} ({unmatched/total_rows*100:.1f}%)")
    print(f"Total time: {elapsed_str}")
    print(f"Output: {output_csv}")
    print(f"{'='*60}")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Annotate peptide DNA origin + gene name")
    ap.add_argument('input_csv', help='input CSV (must contain an `nt` column)')
    ap.add_argument('output_csv', help='annotated output CSV')
    ap.add_argument('sample_name', nargs='?', default='', help='sample name for display')
    ap.add_argument('--human_fasta', required=True, help='human peptide-DNA library FASTA')
    ap.add_argument('--virus_fasta', required=True, help='virome peptide-DNA library FASTA')
    ap.add_argument('--proteome_dir', default='.', help='directory with UniProt proteome FASTAs')
    ap.add_argument('--proteome_files', nargs='+',
                    default=['10239_viruses_reviewed.fasta', 'human_reviewed_unreviewed.fasta'],
                    help='UniProt proteome FASTA file names')
    ap.add_argument('--flank5', default=DEFAULT_FLANK_5, help="5' constant flank")
    ap.add_argument('--flank3', default=DEFAULT_FLANK_3, help="3' constant flank")
    args = ap.parse_args()

    align_to_reference(args.input_csv, args.output_csv, args.sample_name,
                       human_fasta=args.human_fasta, virus_fasta=args.virus_fasta,
                       proteome_dir=args.proteome_dir, proteome_files=args.proteome_files,
                       flank5=args.flank5, flank3=args.flank3)
