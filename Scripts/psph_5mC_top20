#!/usr/bin/env python3
"""
# Pseudomonas top-20 5mC sites 
"""
import pandas as pd

# Loading and correcting column meanings
df = pd.read_csv("psph_top20_5mC_sites.tsv", sep="\t", header=None,
                 names=["contig", "start", "end", "percent_mod", "coverage"])

# Sanity check
assert (df["end"] - df["start"] == 1).all()

# taking selective fields
df["frac_mod"] = df["percent_mod"] / 100.0
df["n_mod"]    = (df["frac_mod"] * df["coverage"]).round().astype(int)

# Showing clean table
cols = ["contig","start","end","coverage","n_mod","frac_mod","percent_mod"]
print(df[cols].to_string(index=False))

df[cols].to_csv("psph_top20_5mC_sites_fixed.csv", index=False)
df[["contig","start","end"]].to_csv("psph_top20_5mC_sites.bed",
                                    sep="\t", header=False, index=False)
