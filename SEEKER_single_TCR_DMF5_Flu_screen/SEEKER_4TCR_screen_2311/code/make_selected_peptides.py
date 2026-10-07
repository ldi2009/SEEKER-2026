# -*- coding: utf-8 -*-
"""Build Selected_Peptides.csv (2026-10-07).

Criteria = the verified criteria of the four results/seqlogos/ PDFs (see README table):
  DMF5: top200, r3_r2_ratio>=4, deduplicated -> 147
  1G4:  top150, r3_r2_ratio>=4               -> 105
  EBV:  top50,  r3_r2_ratio>=10              -> 48
  Flu:  top25,  r3_r2_ratio>=40              -> 6
306 peptides in total, no overlap between TCRs.
"""
import os
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE, "data")

TCR_CONFIG = {
    "DMF5": ("top10000_combine_9mer_2311-4TCR-DMF5-SCD-round123.csv", 200, 4, True),
    "1G4":  ("top10000_combine_9mer_2311-4TCR-1G4-SCD-round123.csv", 150, 4, False),
    "EBV":  ("top10000_combine_9mer_2311-4TCR-EBV-SCD-round123.csv", 50, 10, False),
    "Flu":  ("top10000_combine_9mer_2311-4TCR-Flu-SCD-round123.csv", 25, 40, False),
}


def select_peptides(csv_path, top, ratio, dedup):
    df = pd.read_csv(csv_path)
    s = df.iloc[:top]
    s = s[~s["aa"].str.contains(r"\*", na=False)]
    s = s[s["r3_r2_ratio"] >= ratio]
    seqs = s["aa"].dropna().astype(str).tolist()
    seqs = [x for x in seqs if x]
    mcl = max(set(map(len, seqs)), key=lambda L: sum(1 for x in seqs if len(x) == L))
    out = [x for x in seqs if len(x) == mcl]
    if dedup:
        seen, dd = set(), []
        for x in out:
            if x not in seen:
                seen.add(x)
                dd.append(x)
        out = dd
    return out


if __name__ == "__main__":
    rows = []
    for tcr, (fname, top, ratio, dedup) in TCR_CONFIG.items():
        seqs = select_peptides(os.path.join(DATA_DIR, fname), top, ratio, dedup)
        print(f"[{tcr}] top{top}/r>={ratio}{'/dedup' if dedup else ''} -> {len(seqs)} peptides")
        rows.extend((s, tcr) for s in seqs)
    out = pd.DataFrame(rows, columns=["aa", "tcr_name"])
    assert out["aa"].is_unique, "duplicate peptides across TCRs"
    out.to_csv(os.path.join(DATA_DIR, "Selected_Peptides.csv"), index=False)
    print(f"total {len(out)} peptides ({out['aa'].nunique()} unique), written to data/Selected_Peptides.csv")
