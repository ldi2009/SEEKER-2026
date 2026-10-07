# -*- coding: utf-8 -*-
"""
Generate legible SVG versions of Fig11 and Fig14.

Fig11 (simplified): Top 3 TCR59 + top 6 TCR63 mimicry pairs.
       Minimal layout - TCR label, two aligned sequences with
       match/motif coloring, gene names, score. Nothing else.
Fig14 (large fonts): Mimicry network heatmap (virus family x human gene).

Cell-style: Arial, TCR59=#66C2A5, TCR63=#FC8D62.
"""
import json, csv, re, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from collections import Counter, defaultdict

plt.rcParams.update({
    'font.family': 'Arial',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'svg.fonttype': 'none',
})

COLOR_TCR59 = '#66C2A5'
COLOR_TCR63 = '#FC8D62'
COLOR_MOTIF = '#FF5722'
COLOR_MATCH_BG = '#F5F5F5'

BASE = os.path.dirname(os.path.abspath(__file__))
data_file = os.path.join(BASE, 'mimicry_pairs.json')
with open(data_file, 'r', encoding='utf-8') as f:
    pairs = json.load(f)

# Full organism names from annotated CSV (lookup by tcr + virus sequence)
csv_rows = list(csv.DictReader(
    open(os.path.join(BASE, 'virus_peptides_annotated_v4.csv'),
         encoding='utf-8-sig')))
organism_lookup = {}
for r in csv_rows:
    organism_lookup.setdefault((r['tcr'], r['sequence']), r['organism'])

def virus_species(pair):
    org = organism_lookup.get((pair['virus_tcr'], pair['virus_seq']),
                               pair.get('virus_organism', ''))
    return re.sub(r'\s*\(.*$', '', org).strip()

plot_dir = BASE
os.makedirs(plot_dir, exist_ok=True)

# ============================================================
# Figure 11 (simplified): Top 5 TCR59 + Top 5 TCR63 pairs
# Common human-infecting viruses only, one pair per virus species
# ============================================================
print("=== Fig11 Simplified ===")

# annotated CSV lookup: (tcr, sequence) -> row
csv_idx = {}
for r in csv_rows:
    csv_idx.setdefault((r['tcr'], r['sequence']), r)

COMMON_HUMAN_VIRUSES = {
    'Rotavirus A', 'Human herpesvirus 2', 'Human herpesvirus 1',
    'Human cytomegalovirus', 'Human papillomavirus 30', 'Echovirus 12',
    'Torque teno virus', 'Torque teno midi virus 1',
    'Hepatitis delta virus genotype I', 'Yellow fever virus',
}

def full_gene_name(pair):
    r = csv_idx.get((pair['virus_tcr'], pair['virus_seq']))
    if r:
        return r['gene'] if '_' in r['gene'] else r['uniprot_entry_id']
    return pair['virus_gene']

best_per_species = {}
for p in pairs:
    r = csv_idx.get((p['virus_tcr'], p['virus_seq']))
    if not r:
        continue
    if p['human_tcr'] != p['virus_tcr']:   # exclude cross-TCR pairs
        continue
    org = r['organism'].split(' (')[0]
    if org not in COMMON_HUMAN_VIRUSES:
        continue
    key = (p['virus_tcr'], org)
    if key not in best_per_species or p['combined'] > best_per_species[key]['combined']:
        best_per_species[key] = p

# user-selected pairs: TCR59 #1-2, TCR63 #1-4
SELECTED_KEYS = [
    ('TCR59', 'VP2_ROTHS', 'MYO5A'),
    ('TCR59', 'HEPA_HHV2H', 'PRODH'),
    ('TCR63', 'SHDAG_HDVS2', 'ADNP'),
    ('TCR63', 'CAPSD_TTVV9', 'CACNA1H'),
    ('TCR63', 'VL2_HPV30', 'MRPS25'),
    ('TCR63', 'UL87_HCMVA', 'CYP26C1'),
]
selected = []
for tcr, vg, hg in SELECTED_KEYS:
    cands = [p for p in pairs if p['virus_tcr'] == tcr and p['human_gene'] == hg
             and p['human_tcr'] == p['virus_tcr']
             and full_gene_name(p) == vg]
    cands.sort(key=lambda x: -x['combined'])
    selected.append(cands[0])

def virus_label(pair):
    r = csv_idx.get((pair['virus_tcr'], pair['virus_seq']))
    if r:
        gene_full = r['gene'] if '_' in r['gene'] else r['uniprot_entry_id']
        base = gene_full.split('_')[0]
        org = r['organism'].split(' (')[0]
        if ' genotype' in org:
            org = org[:org.index(' genotype')]
        return f"{base} ({org})"
    return pair['virus_gene']

n = len(selected)
row_h = 3.4        # height per pair block
seq_w = 1.0        # width per residue box
max_len = max(max(len(p['virus_seq']), len(p['human_seq'])) for p in selected)
gene_x = max_len * seq_w + 1.2   # gene labels start after longest sequence

max_label = max(max(len(virus_label(p)), len(p['human_gene'])) for p in selected)
fig_w = (gene_x + max_label * 0.95 + 1.6) * 0.36
fig_h = n * row_h * 0.19 + 1.0

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_xlim(-1.8, gene_x + max_label * 0.95 + 1.2)
ax.set_ylim(-1.2, n * row_h + 0.6)
ax.axis('off')

def draw_seq(seq, other, y, motif, tcr_color):
    m_start = seq.find(motif) if motif else -1
    for j, aa in enumerate(seq):
        in_motif = m_start >= 0 and m_start <= j < m_start + len(motif)
        is_match = j < len(other) and aa == other[j]
        if in_motif or is_match:
            color, tcol = tcr_color, 'white'
        else:
            color, tcol = COLOR_MATCH_BG, '#555'
        ax.add_patch(Rectangle((j * seq_w, y - 0.5), seq_w - 0.12, 1.0,
                               facecolor=color, edgecolor='none'))
        ax.text(j * seq_w + (seq_w - 0.12) / 2, y, aa, ha='center',
                va='center', fontsize=11, fontweight='bold', color=tcol)

for i, pair in enumerate(selected):
    y_top = n * row_h - i * row_h - 1.0
    y_bot = y_top - 1.55
    v_seq, h_seq = pair['virus_seq'], pair['human_seq']
    motif = pair.get('common_motif', '') or ''
    tcr_color = COLOR_TCR59 if pair['virus_tcr'] == 'TCR59' else COLOR_TCR63

    # TCR label
    ax.text(-1.8, (y_top + y_bot) / 2, pair['virus_tcr'], fontsize=10,
            fontweight='bold', va='center', ha='left', color=tcr_color)

    draw_seq(v_seq, h_seq, y_top, motif, tcr_color)
    draw_seq(h_seq, v_seq, y_bot, motif, tcr_color)

    # Virus gene (blue) and human gene (colored by self-antigen TCR origin)
    ax.text(gene_x, y_top + 0.12, virus_label(pair), fontsize=10,
            va='center', ha='left', color='#1565C0', fontstyle='italic',
            fontweight='bold')
    h_tcr_color = COLOR_TCR59 if pair['human_tcr'] == 'TCR59' else COLOR_TCR63
    ax.text(gene_x, y_bot, pair['human_gene'], fontsize=10,
            va='center', ha='left', color=h_tcr_color, fontweight='bold')

# Legend
ly = -0.8
ax.add_patch(Rectangle((0, ly - 0.4), 0.88, 0.8, facecolor=COLOR_TCR59,
                       edgecolor='none'))
ax.text(1.15, ly, 'TCR59 match / motif', fontsize=9, va='center')
ax.add_patch(Rectangle((6.6, ly - 0.4), 0.88, 0.8, facecolor=COLOR_TCR63,
                       edgecolor='none'))
ax.text(7.75, ly, 'TCR63 match / motif', fontsize=9, va='center')
ax.add_patch(Rectangle((13.2, ly - 0.4), 0.88, 0.8, facecolor=COLOR_MATCH_BG,
                       edgecolor='none'))
ax.text(14.35, ly, 'Mismatch', fontsize=9, va='center')

plt.savefig(os.path.join(plot_dir, 'Fig11_Mimicry_Top_Pairs_Selected.svg'),
            facecolor='white')
plt.savefig(os.path.join(plot_dir, 'Fig11_Mimicry_Top_Pairs_Selected.png'),
            dpi=300, facecolor='white')
plt.close()
print(f"  Saved simplified Fig11 (n={n})")

# ============================================================
# Figure 14 (extra large fonts)
# ============================================================
print("=== Fig14 Extra Large ===")

# Per-TCR top genes: 5 from TCR59 + 5 from TCR63
def top_genes(tcr, k):
    c = Counter(p['human_gene'] for p in pairs
                if (p.get('human_tcr') or p['virus_tcr']) == tcr)
    return [g for g, _ in c.most_common(k)]

top_hg = top_genes('TCR59', 5) + top_genes('TCR63', 5)

# Only families with >= 1 pair among selected genes, ordered by total pairs
fam_in_sel = Counter(p['virus_family'] for p in pairs
                     if p['human_gene'] in top_hg)
top_vf = [f for f, _ in fam_in_sel.most_common()]

matrix = np.zeros((len(top_vf), len(top_hg)))
for p in pairs:
    vf, hg = p['virus_family'], p['human_gene']
    if vf in top_vf and hg in top_hg:
        matrix[top_vf.index(vf), top_hg.index(hg)] += 1

# TCR origin of each human gene (self-antigen) and virus family
gene_tcr = defaultdict(set)
fam_tcr = defaultdict(set)
for p in pairs:
    gene_tcr[p['human_gene']].add(p.get('human_tcr') or p['virus_tcr'])
    fam_tcr[p['virus_family']].add(p['virus_tcr'])

def tcr_label_color(tcrs):
    if len(tcrs) > 1:
        return '#333333'          # shared by both TCRs
    return COLOR_TCR59 if 'TCR59' in tcrs else COLOR_TCR63

fig, ax = plt.subplots(figsize=(16, 12))
cmap = LinearSegmentedColormap.from_list(
    'whites', ['#FFFFFF', '#FFF3E0', '#FFE0B2', '#FFCC80', '#FFA726', '#E65100'])
im = ax.imshow(matrix, aspect='auto', cmap=cmap, interpolation='nearest')

ax.set_xticks(np.arange(len(top_hg)))
ax.set_xticklabels(top_hg, fontsize=22, rotation=45, ha='right',
                   fontweight='bold')
ax.set_yticks(np.arange(len(top_vf)))
ax.set_yticklabels(top_vf, fontsize=22, fontweight='bold')

# Color gene / family labels by TCR origin
for lbl, gene in zip(ax.get_xticklabels(), top_hg):
    lbl.set_color(tcr_label_color(gene_tcr[gene]))
for lbl, fam in zip(ax.get_yticklabels(), top_vf):
    lbl.set_color(tcr_label_color(fam_tcr[fam]))

for i in range(len(top_vf)):
    for j in range(len(top_hg)):
        if matrix[i, j] > 0:
            ax.text(j, i, int(matrix[i, j]), ha='center', va='center',
                    fontsize=22, fontweight='bold',
                    color='black' if matrix[i, j] < 5 else 'white')

# Separator between TCR59 gene block and TCR63 gene block
n59 = 5
ax.axvline(n59 - 0.5, color='#333333', linewidth=3)

ax.set_xlabel('Human Gene', fontsize=22, labelpad=14)
ax.set_ylabel('Virus Family', fontsize=22, labelpad=14)
ax.set_title('Mimicry Network: Virus Family vs Human Gene',
             fontsize=24, fontweight='bold', pad=18)
ax.tick_params(length=6, width=2.2)

ax.set_xticks(np.arange(-0.5, len(top_hg), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(top_vf), 1), minor=True)
ax.grid(which='minor', color='white', linewidth=4)
ax.tick_params(which='minor', length=0)

cbar = plt.colorbar(im, ax=ax, shrink=0.65, pad=0.02)
cbar.set_label('Number of Pairs', fontsize=20)
cbar.ax.tick_params(labelsize=18, width=2.2)
cbar.outline.set_linewidth(2)

# TCR-origin legends (below the heatmap, outside data area)
leg_gene = ax.legend(
    handles=[Patch(facecolor=COLOR_TCR59, label='TCR59'),
             Patch(facecolor=COLOR_TCR63, label='TCR63')],
    title='Self-antigen origin (x-axis)',
    loc='upper left', bbox_to_anchor=(0.0, -0.22), ncol=2,
    fontsize=15, title_fontsize=16, framealpha=0.95,
    edgecolor='#666', borderpad=0.8)
ax.add_artist(leg_gene)
ax.legend(
    handles=[Patch(facecolor=COLOR_TCR59, label='TCR59'),
             Patch(facecolor=COLOR_TCR63, label='TCR63'),
             Patch(facecolor='#333333', label='Both')],
    title='Virus peptide origin (y-axis)',
    loc='upper right', bbox_to_anchor=(1.0, -0.22), ncol=3,
    fontsize=15, title_fontsize=16, framealpha=0.95,
    edgecolor='#666', borderpad=0.8)

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'Fig14_Mimicry_Network_Large.svg'),
            facecolor='white')
plt.savefig(os.path.join(plot_dir, 'Fig14_Mimicry_Network_Large.png'),
            dpi=300, facecolor='white')
plt.close()
print("  Saved extra-large Fig14")

print("\nDone.")
