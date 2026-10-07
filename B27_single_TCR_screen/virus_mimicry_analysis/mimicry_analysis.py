# -*- coding: utf-8 -*-
"""
Molecular mimicry analysis: pair activation-positive virus peptides with
activation-positive human (self) peptides of the same TCR screen.

Inputs (relative to this script):
  virus_peptides_annotated_v4.csv   virus peptides + UniProt annotation + single_pct
  data/TCR59_peptide_validation_summary.csv
  data/TCR63_peptide_validation_summary.csv

Activation flags in the inputs (addflow validation, > 5% positive / > 30% strong):
  data/TCR59_peptide_validation_summary.csv : single- OR high-copy addflow
  data/TCR63_peptide_validation_summary.csv : single-copy addflow (no high-copy data)
  virus_peptides_annotated_v4.csv           : single-copy addflow

A virus-self pair is called when the peptide lengths differ by <= 6 and either
  combined score > 0.35  (0.6 x Levenshtein similarity + 0.4 x 3-mer Jaccard)
  or the two peptides share a contiguous motif of >= 4 amino acids.

Output: mimicry_pairs.json (sorted by combined score, descending).
"""
import csv, json, os

BASE = os.path.dirname(os.path.abspath(__file__))
POS_THR = 5.0
MOTIF_MIN = 4
COMBINED_THR = 0.35
LEN_DIFF_MAX = 6

def read_csv(path):
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def to_bool(x):
    return str(x).strip().lower() == 'true'

# ---- self (human) positive peptides, TCR59 first then TCR63 ----
human_peptides = []
for tcr, fname in [('TCR59', 'TCR59_peptide_validation_summary.csv'),
                   ('TCR63', 'TCR63_peptide_validation_summary.csv')]:
    for row in read_csv(os.path.join(BASE, 'data', fname)):
        if row.get('Origin', '').strip().lower() != 'human':
            continue
        if not to_bool(row.get('Positive', 'False')):
            continue
        seq = row.get('Sequence', '').strip()
        if not seq:
            continue
        human_peptides.append({
            'tcr': tcr,
            'sequence': seq,
            'gene': row.get('Gene', ''),
            'peptide_number': row.get('Peptide#', ''),
            'single_pct': float(row.get('Single(%)', 0) or 0),
            'strong': to_bool(row.get('Strong', 'False')),
        })

print(f"Human positive peptides: {len(human_peptides)}")
for tcr in ('TCR59', 'TCR63'):
    n = sum(1 for p in human_peptides if p['tcr'] == tcr)
    print(f"  {tcr}: {n}")

# ---- virus positive peptides (single > 5) ----
virus_pos = []
for row in read_csv(os.path.join(BASE, 'virus_peptides_annotated_v4.csv')):
    if not to_bool(row.get('positive', 'False')):
        continue
    virus_pos.append({
        'tcr': row['tcr'],
        'sequence': row['sequence'],
        'gene': row.get('gene', ''),
        'organism': row.get('organism', ''),
        'virus_family': row.get('virus_family', ''),
        'single_pct': float(row.get('single_pct', 0) or 0),
    })
print(f"Virus positive peptides: {len(virus_pos)}")
for tcr in ('TCR59', 'TCR63'):
    n = sum(1 for p in virus_pos if p['tcr'] == tcr)
    print(f"  {tcr}: {n}")

# ---- sequence metrics ----
def shared_kmers(seq1, seq2, k=3):
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1)-k+1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2)-k+1))
    return len(kmers1 & kmers2), len(kmers1 | kmers2)

def levenshtein(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

def similarity_score(seq1, seq2):
    dist = levenshtein(seq1, seq2)
    max_len = max(len(seq1), len(seq2))
    return 1.0 - dist / max_len if max_len > 0 else 0

def find_common_motif(seq1, seq2, min_len=3):
    m, n = len(seq1), len(seq2)
    longest = ""
    for i in range(m):
        for j in range(n):
            k = 0
            while i+k < m and j+k < n and seq1[i+k] == seq2[j+k]:
                k += 1
            if k >= min_len and k > len(longest):
                longest = seq1[i:i+k]
    return longest

# ---- pair every virus positive with every human positive ----
mimicry_pairs = []
for vp in virus_pos:
    v_seq = vp['sequence']
    for hp in human_peptides:
        h_seq = hp['sequence']
        if abs(len(v_seq) - len(h_seq)) > LEN_DIFF_MAX:
            continue
        sim = similarity_score(v_seq, h_seq)
        shared, union = shared_kmers(v_seq, h_seq, k=3)
        jaccard = shared / union if union > 0 else 0
        motif = find_common_motif(v_seq, h_seq, min_len=3)
        combined = sim * 0.6 + jaccard * 0.4
        if combined > COMBINED_THR or (motif and len(motif) >= MOTIF_MIN):
            mimicry_pairs.append({
                'virus_tcr': vp['tcr'],
                'virus_seq': v_seq,
                'virus_gene': vp['gene'],
                'virus_organism': vp['organism'][:40],
                'virus_family': vp['virus_family'],
                'human_tcr': hp['tcr'],
                'human_seq': h_seq,
                'human_gene': hp['gene'],
                'similarity': round(sim, 3),
                'jaccard': round(jaccard, 3),
                'combined': round(combined, 3),
                'shared_kmers': shared,
                'common_motif': motif,
                'virus_single': vp['single_pct'],
                'human_single': hp['single_pct'],
            })

mimicry_pairs.sort(key=lambda x: -x['combined'])

print(f"\nMimicry pairs (combined > {COMBINED_THR} or motif >= {MOTIF_MIN}): {len(mimicry_pairs)}")
for tcr in ('TCR59', 'TCR63'):
    n = sum(1 for p in mimicry_pairs if p['virus_tcr'] == tcr)
    print(f"  with {tcr} virus peptides: {n}")

out = os.path.join(BASE, 'mimicry_pairs.json')
with open(out, 'w', encoding='utf-8') as f:
    json.dump(mimicry_pairs, f, indent=2, ensure_ascii=False)
print(f"Saved to {out}")
