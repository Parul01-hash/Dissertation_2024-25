#!/usr/bin/env python3
"""
Input:  psph_5mC_cov30_pm10.fa  
"""

# Validate centered 5mC windows
# emit PFM + annotated tables, and split FASTAs by joint motif. 
import re
from collections import Counter
import pandas as pd

FA_PATH = "psph_5mC_cov30_pm10.fa"   

#  Reading FASTA 
records = []
with open(FA_PATH) as f:
    header = None
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            header = line[1:]
        else:
            seq = line.upper()
            # parse header fields
            m = re.match(
                r'^(?P<contig>[^:]+):(?P<start>\d+)-(?P<end>\d+)\|ref=(?P<ref>[^|]+)\|cov=(?P<cov>[\d.]+)\|pct=(?P<pct>[\d.]+)\|center=(?P<center>[ACGT])$',
                header
            )
            if not m:
                raise ValueError(f"Unexpected FASTA header format: {header}")
            rec = m.groupdict()
            rec["start"] = int(rec["start"])
            rec["end"]   = int(rec["end"])
            rec["cov"]   = float(rec["cov"])
            rec["pct"]   = float(rec["pct"])
            rec["seq"]   = seq
            records.append(rec)

assert records, "No sequences read from FASTA."

# Sanity check
L = len(records[0]["seq"])
assert all(len(r["seq"]) == L for r in records), "Not all windows have equal length."
center = L // 2
center_rate = sum(r["seq"][center] == "C" for r in records) / len(records)
print(f"Windows: {len(records)} | length: {L} | center 'C' rate: {center_rate:.3f} (expect ~1.000)")

# Context counters
left1  = Counter(r["seq"][center-1] for r in records)
right1 = Counter(r["seq"][center+1] for r in records)
right23 = Counter(r["seq"][center+1:center+4] for r in records)

leftC_count  = left1["C"]
rightGAG_count = sum(r["seq"][center+1:center+4] == "GAG" for r in records)
both_count = sum((r["seq"][center-1] == "C") and (r["seq"][center+1:center+4] == "GAG") for r in records)
n = len(records)

print("\nContext counts:")
print(f"  C at -1:            {leftC_count}/{n}  ({leftC_count/n:.1%})")
print(f"  GAG at +1..+3:      {rightGAG_count}/{n}  ({rightGAG_count/n:.1%})")
print(f"  BOTH (C@-1 & GAG):  {both_count}/{n}  ({both_count/n:.1%})")

#  Position by nucleotide PFM 
bases = "ACGT"
pfm_rows = []
for i in range(L):
    c = Counter(r["seq"][i] for r in records)
    pfm_rows.append({b: c.get(b, 0) for b in bases})
pfm_df = pd.DataFrame(pfm_rows)
pfm_df.index.name = "pos"  
pfm_df.to_csv("psph_5mC_cov30_PFM.csv")
print("\nSaved PFM → psph_5mC_cov30_PFM.csv (columns A,C,G,T; center index =", center, ")")

# Annotate each window and save a tidy table
df = pd.DataFrame(records)
df["leftC"] = df["seq"].str.slice(center-1, center) == "C"
df["rightGAG"] = df["seq"].str.slice(center+1, center+4) == "GAG"
df["both"] = df["leftC"] & df["rightGAG"]
df_out = df[["contig","start","end","cov","pct","leftC","rightGAG","both","seq"]].copy()
df_out.rename(columns={"cov":"coverage","pct":"percent_mod"}, inplace=True)
df_out.to_csv("psph_5mC_windows_annotated.csv", index=False)
print("Saved annotated table → psph_5mC_windows_annotated.csv")

#  Writing two FASTAs: those matching both vs the rest 
def write_fa(subdf, path):
    with open(path, "w") as out:
        for _, r in subdf.iterrows():
            hdr = f"{r['contig']}:{r['start']}-{r['end']}|cov={int(r['coverage'])}|pct={r['percent_mod']:.2f}"
            out.write(f">{hdr}\n{r['seq']}\n")

df_both = df_out[df_out["both"]]
df_other = df_out[~df_out["both"]]
write_fa(df_both,  "psph_5mC_both_Cat-1_GAG.fa")
write_fa(df_other, "psph_5mC_not_both.fa")
print(f"Wrote FASTAs: psph_5mC_both_Cat-1_GAG.fa ({len(df_both)})  |  psph_5mC_not_both.fa ({len(df_other)})")
