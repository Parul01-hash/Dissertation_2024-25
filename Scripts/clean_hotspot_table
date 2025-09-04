#!/usr/bin/env python3
"""
# Maps BED contigs to bc_reference.fasta and fetches ref base per site.
# Normalizes frac_mod, runs strand/base mismatch QC.
# Calls hotspots , computes binomial p-values and BH q-values.

Inputs 
# bc_reference.fasta  
"""

# Map to FASTA
# normalizing fractions, strand/base QC, compute FDR, preview, and save. 
import re, numpy as np, pandas as pd
from pyfaidx import Fasta
from IPython.display import display


assert 'df' in globals() and isinstance(df, pd.DataFrame), "I can't find a DataFrame named `df`."

#  Map BED contigs -> FASTA keys
ref = Fasta("bc_reference.fasta")
REF_KEYS = set(ref.keys())

def to_ref_name(name: str):
    base = name.split()[0]
    nobase = re.sub(r'\.\d+$', '', base)            
    for cand in (base, nobase, f"NZ_{base}", f"NZ_{nobase}"):
        if cand in REF_KEYS: return cand
    return np.nan

if "chrom_ref" not in df.columns:
    df["chrom_ref"] = df["chrom"].map(to_ref_name)
assert df["chrom_ref"].notna().all(), "Some contigs didn't map to the FASTA."

# Normalize frac_mod 
if (df["frac_mod"] > 1).any():        
    df["frac_mod"] = df["frac_mod"] / 100.0

df["percent_mod"] = (df["frac_mod"] * 100).round(2)
df["n_mod_check"] = (df["frac_mod"] * df["valid_cov"]).round().astype(int)

# Reference base at 0-based BED start 
def fetch_ref_base(chrom_ref, start_0b):
    s = int(start_0b)
    return ref[chrom_ref][s:s+1].seq.upper()
if "ref_base" not in df.columns:
    df["ref_base"] = df.apply(lambda r: fetch_ref_base(r["chrom_ref"], r["position"]), axis=1)

# Expected base from mod + strand; strand-aware mismatch 
def expected_base(mod, strand):
    if mod == "6mA": return "A" if strand != "-" else "T"
    if mod == "5mC": return "C" if strand != "-" else "G"
    return "N"

df["expected_base"] = df.apply(lambda r: expected_base(r["modification"], r["strand"]), axis=1)
is_nt = df["ref_base"].isin(list("ACGT"))
mismatch = (df.loc[is_nt, "ref_base"] != df.loc[is_nt, "expected_base"]).mean()
print(f"Base/mod mismatch (strand-aware): {mismatch:.3%}")

# QC subset & useful intervals 
keep = ((df["modification"]=="6mA") & df["ref_base"].isin(["A","T"])) | \
       ((df["modification"]=="5mC") & df["ref_base"].isin(["C","G"]))
df_qc = df[keep].copy()

# high-confidence hotspots 
hot = df_qc.query("valid_cov >= 30 and frac_mod >= 0.10").copy()

p0 = {"6mA": 0.0038, "5mC": 0.0021}
from math import comb
def binom_sf(k, n, p):
    s = 0.0
    for i in range(int(k), int(n)+1):
        s += comb(int(n), i) * (p**i) * ((1-p)**(int(n)-i))
    return min(1.0, s)

if len(hot):
    hot["pval"] = hot.apply(lambda r: binom_sf(r["n_mod"], r["valid_cov"], p0.get(r["modification"], 1e-9)), axis=1)
    m = len(hot)
    hot = hot.sort_values("pval").assign(
        rank=lambda d: np.arange(1, len(d)+1),
        qval=lambda d: (d["pval"] * m / d["rank"]).clip(upper=1)
    ).sort_values(["qval","frac_mod"])
    sig = hot.query("qval <= 0.05").copy()
else:
    sig = hot.copy()

# Showing & save a clean preview 
cols = ["chrom","position","strand","modification","valid_cov","n_mod","frac_mod","percent_mod","ref_base","expected_base"]
cols = [c for c in cols if c in df_qc.columns]
preview = df_qc[cols].head(50).style.format({
    "frac_mod": "{:.3f}",
    "percent_mod": "{:.2f}",
    "valid_cov": "{:,}",
    "n_mod": "{:,}",
})
display(preview)

df_qc.to_csv("hotspots_clean.csv", index=False)
print("Saved cleaned table to hotspots_clean.csv")
print(f"Rows (QC): {len(df_qc):,} | Hot (>=30 cov & >=10%): {len(hot):,} | Significant (q<=0.05): {len(sig):,}")
