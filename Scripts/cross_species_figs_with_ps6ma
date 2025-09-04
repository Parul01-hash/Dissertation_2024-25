#!/usr/bin/env python3
"""
# Cross-species figures with Pseudomonas 6mA.

# Inputs 
# psph_hotspots_6mA.bed                
# hotspots_hot.csv                     
# psph_5mC_windows_annotated.csv  
# psph_5mC_windows_annot.csv           
# psph_5mC_cov30_annot.csv             

"""

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from pathlib import Path

# Building Pseudomonas 6mA table from BED 
BED_6MA = "psph_hotspots_6mA.bed"   # change if your filename differs
CSV_6MA = "psph_6mA_hotspots.csv"

def read_bed_auto(path, mod_label, species_label):
    df = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False)
    df.columns = [f"c{i}" for i in range(df.shape[1])]
    df = df.rename(columns={"c0":"chrom","c1":"start","c2":"end"})
    num_cols = [c for c in df.columns if c != "chrom"]

    best = None
    samp = df.head(50000)
    for cov_col in num_cols:
        cov = pd.to_numeric(samp[cov_col], errors="coerce")
        if cov.median() < 5:  # unlikely to be coverage
            continue
        for n_col in num_cols:
            if n_col == cov_col: 
                continue
            n = pd.to_numeric(samp[n_col], errors="coerce")
            mask = (~cov.isna()) & (~n.isna()) & (cov > 0) & (n >= 0)
            if mask.mean() < 0.5:
                continue
            ok = ((n[mask]/cov[mask]).between(0,1)).mean()
            cand = (ok, cov.median(), cov_col, n_col)
            if best is None or cand > best:
                best = cand
    if best is None:
        raise ValueError(f"Couldn't infer coverage/n_mod from {path}")
    _, _, cov_col, n_col = best

    out = df[["chrom","start","end"]].copy()
    out["coverage"]     = pd.to_numeric(df[cov_col], errors="coerce")
    out["n_mod"]        = pd.to_numeric(df[n_col],   errors="coerce")
    out["percent_mod"]  = (out["n_mod"]/out["coverage"])*100.0
    out["modification"] = mod_label
    out["species"]      = species_label
    return out[["species","modification","percent_mod","coverage","chrom","start","end"]]

if not Path(CSV_6MA).exists():
    assert Path(BED_6MA).exists(), f"Can't find {BED_6MA}. If your 6mA file has another name, update BED_6MA."
    ps6 = read_bed_auto(BED_6MA, "6mA", "Pseudomonas")
    ps6.to_csv(CSV_6MA, index=False)
    print(f"Built {CSV_6MA} with {len(ps6)} rows")
else:
    print(f"Found existing {CSV_6MA}")

# Loading all datasets 
# Bacillus
bc = pd.read_csv("hotspots_hot.csv")[["percent_mod","coverage","region","modification"]].copy()
bc["species"] = "Bacillus"

# Pseudomonas 5mC
ps5_path = "psph_5mC_windows_annotated.csv" if Path("psph_5mC_windows_annotated.csv").exists() else "psph_5mC_windows_annot.csv"
ps5 = pd.read_csv(ps5_path)
if "contig" in ps5.columns: ps5.rename(columns={"contig":"chrom"}, inplace=True)

if "region" not in ps5.columns and Path("psph_5mC_cov30_annot.csv").exists():
    ann = pd.read_csv("psph_5mC_cov30_annot.csv")[["contig","start","end","region"]].rename(columns={"contig":"chrom"})
    ps5 = ps5.merge(ann, on=["chrom","start","end"], how="left")
ps5 = ps5[["percent_mod","coverage","region"]].copy()
ps5["species"] = "Pseudomonas"; ps5["modification"] = "5mC"

# Pseudomonas 6mA 
ps6 = pd.read_csv(CSV_6MA)
ps6["region"] = np.nan  

# Combine
comp = pd.concat([
    bc,
    ps5[["percent_mod","coverage","region","modification","species"]],
    ps6[["percent_mod","coverage","region","modification","species"]],
], ignore_index=True)

# If percent_mod was stored as 0–1 anywhere, convert to 0–100
if comp["percent_mod"].dropna().quantile(0.9) <= 1.0:
    comp["percent_mod"] = comp["percent_mod"] * 100.0

# Making compact figures
plt.rcParams.update({"figure.figsize":(4.6,3.4),"figure.dpi":300,"savefig.dpi":600,"font.size":9})

# Fig1: counts 
MIN_PCT, MIN_COV = 5.0, 30
hot = comp[(comp["coverage"]>=MIN_COV) & (comp["percent_mod"]>=MIN_PCT)]
cnt = (hot.groupby(["species","modification"]).size().rename("count").reset_index())
if len(cnt):
    piv = cnt.pivot(index="species", columns="modification", values="count").fillna(0)
    x = np.arange(len(piv.index)); mods = list(piv.columns); w = 0.8/max(1,len(mods))
    plt.figure()
    for i,m in enumerate(mods):
        vals = piv[m].values
        bars = plt.bar(x+i*w-(len(mods)-1)*w/2, vals, width=w, label=m)
        for b,v in zip(bars, vals): plt.text(b.get_x()+b.get_width()/2, v, f"{int(v)}", ha="center", va="bottom", fontsize=8)
    plt.xticks(x, piv.index); plt.ylabel(f"Hotspot count (≥{MIN_PCT:.0f}%, ≥{MIN_COV}×)")
    plt.title("Hotspots per Species × Modification")
    if len(mods)>1: plt.legend()
    plt.tight_layout(); plt.savefig("Cross_Fig1_counts_relaxed.png", bbox_inches="tight"); plt.close()

# Fig2: percent-mod distributions 
groups, labels = [], []
for sp in sorted(comp["species"].unique()):
    for mod in sorted(comp.loc[comp["species"]==sp,"modification"].unique()):
        vals = comp[(comp["species"]==sp)&(comp["modification"]==mod)]["percent_mod"].dropna().values
        if len(vals): groups.append(vals); labels.append(f"{sp} {mod} (n={len(vals)})")
if groups:
    plt.figure(); plt.boxplot(groups, labels=labels, showmeans=True)
    plt.ylabel("Percent modified"); plt.title("Percent Modified per Species × Modification")
    plt.xticks(rotation=15); plt.ylim(0,100); plt.tight_layout(); plt.savefig("Cross_Fig2_percent_mod_box.png", bbox_inches="tight"); plt.close()

# Fig3: region fractions 
if "region" in comp.columns and comp["region"].notna().any():
    reg = (comp.dropna(subset=["region"])
           .groupby(["species","modification","region"]).size().rename("count").reset_index())
    if len(reg):
        reg_piv = reg.pivot_table(index=["species","modification"], columns="region", values="count", fill_value=0)
        frac = reg_piv.div(reg_piv.sum(axis=1), axis=0)
        idx = np.arange(len(frac.index)); bottom = np.zeros(len(idx))
        plt.figure()
        for col in frac.columns:
            vals = frac[col].values
            plt.bar(idx, vals, bottom=bottom, label=col); bottom += vals
        plt.xticks(idx, [f"{a}\n{b}" for a,b in frac.index]); plt.ylabel("Fraction of sites")
        plt.title("Region Distribution per Species × Modification")
        plt.legend(); plt.tight_layout(); plt.savefig("Cross_Fig3_region_stack_small.png", bbox_inches="tight"); plt.close()

print("Done. You should now see Pseudomonas 6mA included in the counts and boxplot.")
