#!/usr/bin/env python3
"""
# Cross-species region stack
# Inputs 
# hotspots_hot.csv
# psph_6mA_hotspots_annot_fixed.csv  
"""

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt

MIN_PCT, MIN_COV = 5.0, 30.0

def _coerce(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def _norm_region(s):
    return s.replace({"gene":"gene_body"}) if s is not None else s

# Bacillus 
assert os.path.exists("hotspots_hot.csv"), "hotspots_hot.csv not found"
bc = pd.read_csv("hotspots_hot.csv")
bc = _coerce(bc, ["percent_mod","coverage"])
bc = bc[(bc["coverage"]>=MIN_COV) & (bc["percent_mod"]>=MIN_PCT)].copy()
bc["region"] = _norm_region(bc["region"])
bc5 = bc[bc["modification"].astype(str).str.upper().eq("5MC")].copy()
bc6 = bc[bc["modification"].astype(str).str.upper().eq("6MA")].copy()
bc5["species"] = "Bacillus"; bc6["species"] = "Bacillus"

# Pseudomonas 6mA 
ps6_path = "psph_6mA_hotspots_annot_fixed.csv" if os.path.exists("psph_6mA_hotspots_annot_fixed.csv") \
           else "psph_6mA_hotspots_annot.csv"
assert os.path.exists(ps6_path), "Pseudomonas 6mA annotated file not found"
ps6 = pd.read_csv(ps6_path)
ps6 = _coerce(ps6, ["percent_mod","coverage"])

if "coverage" not in ps6.columns:    ps6["coverage"] = MIN_COV
if "percent_mod" not in ps6.columns: ps6["percent_mod"] = MIN_PCT
ps6 = ps6[(ps6["coverage"]>=MIN_COV) & (ps6["percent_mod"]>=MIN_PCT)].copy()
if "region" not in ps6.columns:
    raise ValueError("PS 6mA file has no 'region' column. Re-run the annotation step.")
ps6["region"] = _norm_region(ps6["region"])
ps6["species"] = "Pseudomonas"; ps6["modification"] = "5mC"  # temp to keep columns
ps6["modification"] = "6mA"

# Pseudomonas 5mC 
ps5 = pd.DataFrame()
for p in ["psph_5mC_cov30_annot.csv", "psph_5mC_hotspots_full.csv",
          "psph_5mC_windows_annotated.csv", "psph_5mC_windows_annot.csv"]:
    if os.path.exists(p):
        ps5 = pd.read_csv(p); found = p; break

if not ps5.empty:
    
    if "contig" in ps5.columns: ps5.rename(columns={"contig":"chrom"}, inplace=True)
    if "cov" in ps5.columns and "coverage" not in ps5.columns: ps5.rename(columns={"cov":"coverage"}, inplace=True)
    if "percent_mod" not in ps5.columns:
        if "mod_frac" in ps5.columns:
            ps5["percent_mod"] = pd.to_numeric(ps5["mod_frac"], errors="coerce") * 100.0
        elif {"NumModified","Coverage"} <= set(ps5.columns):
            ps5["percent_mod"] = (pd.to_numeric(ps5["NumModified"], errors="coerce") /
                                  pd.to_numeric(ps5["Coverage"],     errors="coerce")) * 100.0
            if "coverage" not in ps5.columns:
                ps5["coverage"] = ps5["Coverage"]
    
    if "coverage" not in ps5.columns:    ps5["coverage"] = MIN_COV
    if "percent_mod" not in ps5.columns: ps5["percent_mod"] = MIN_PCT

    ps5 = _coerce(ps5, ["percent_mod","coverage"])
    ps5 = ps5[(ps5["coverage"]>=MIN_COV) & (ps5["percent_mod"]>=MIN_PCT)].copy()
    if "region" in ps5.columns: ps5["region"] = _norm_region(ps5["region"])
    else:                       ps5["region"] = np.nan
    ps5["species"] = "Pseudomonas"; ps5["modification"] = "5mC"
else:
    found = "(none)"
    print("Note: no Pseudomonas 5mC file found; PS-5mC bar may be missing.")

# plotting 
keep_cols = ["species","modification","region"]
frames = [bc5[keep_cols], bc6[keep_cols], ps6[keep_cols]]
if not ps5.empty: frames.append(ps5[keep_cols])
plot_df = pd.concat(frames, ignore_index=True)
plot_df = plot_df.dropna(subset=["region"])
plot_df.to_csv("region_plot_input.csv", index=False)

print("\nFiles used:")
print(f"  Bacillus: hotspots_hot.csv")
print(f"  Pseudomonas 6mA: {ps6_path}")
print(f"  Pseudomonas 5mC: {found}")

print("\nRows with region by group:")
plot_df["group"] = list(zip(plot_df["species"], plot_df["modification"]))
print(plot_df.groupby("group").size().to_string())

# Plotting 
order = [("Bacillus","5mC"), ("Bacillus","6mA"),
         ("Pseudomonas","5mC"), ("Pseudomonas","6mA")]
plot_df = plot_df[plot_df["group"].isin(order)]
cts = (plot_df.groupby(["group","region"])
       .size().rename("count").reset_index())
piv = cts.pivot(index="group", columns="region", values="count").fillna(0).reindex(order)
frac = piv.div(piv.sum(axis=1), axis=0)

plt.rcParams.update({"figure.figsize":(5.4,3.6), "figure.dpi":300, "savefig.dpi":600, "font.size":9})
x = np.arange(len(frac.index)); bottom = np.zeros(len(frac.index))
for col in ["CDS","gene_body","intergenic","promoter"]:
    if col in frac.columns:
        vals = frac[col].values
        plt.bar(x, vals, bottom=bottom, label=col)
        bottom += vals
plt.xticks(x, [f"{s}\n{m}" for (s,m) in frac.index])
plt.ylabel("Fraction of sites")
plt.title("Region Distribution per Species × Modification (hotspots)")
plt.legend()
plt.tight_layout()
plt.savefig("Cross_Fig3_region_stack_with_PS6mA_v2.png", bbox_inches="tight")
plt.savefig("Cross_Fig3_region_stack_with_PS6mA_v2.pdf", bbox_inches="tight")
print("\nSaved: Cross_Fig3_region_stack_with_PS6mA_v2.[png|pdf]")
