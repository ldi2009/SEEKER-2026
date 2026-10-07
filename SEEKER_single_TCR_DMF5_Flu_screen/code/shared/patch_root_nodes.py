# -*- coding: utf-8 -*-
"""Flu network nodes.csv patch: fix singleton clusters to gray + add the is_anchor
column, matching the batch pipeline output format.
(run after extract_network_data.py)"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
p = os.path.join(REPO, "results", "flu_network", "nodes.csv")

df = pd.read_csv(p)
df.loc[df["cluster"] >= 3, "color"] = "#9E9E9E"
df["is_anchor"] = (df["name"] == "GILGFVFTL").astype(int)
df.to_csv(p, index=False)
print(df.columns.tolist())
print(df["is_anchor"].value_counts().to_dict())
