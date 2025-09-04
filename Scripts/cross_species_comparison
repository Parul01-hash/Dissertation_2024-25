#!/usr/bin/env python3
"""
# cross-species methylation tables and plot.
# input
# bc_hotspots_6mA.bed
# bc_hotspots_5mC.bed
# psph_hotspots_6mA.bed
# psph_5mC_windows_annotated.csv  OR  psph_5mC_windows_annot.csv
# bc_modifications_with_region.csv
# psph_5mC_cov30_annot.csv
"""

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from pathlib import Path

BED_COLS = [
    "Chrom","Start","End","ModType","Coverage","Strand",
    "ThickStart","ThickEnd","ItemRgb","BlockCount","ModFrac","NumModified",
    "Extra1","Extra2","Extra3","Extra4","Extra5","Extra6"
]

def load_pacbio_bed(path, mod_label, species_label):
    df = pd.read_csv(path, sep="\t", header=None, names=BED_COLS, low_memory=False)
    # Compute % correctly
    cov = pd.to_numeric(df["Coverage"], errors="coerce")
    nmd = pd.to_numeric(df["NumModified"], errors="coerce")
    out = pd.DataFrame({
        "species": species_label,
        "modification": mod_label,
        "percent_mod": (nmd / cov) * 100.0,
        "coverage": cov,
        "chrom": df["Chrom"],
        "start": pd.to_numeric(df["Start"], errors="coerce"),
        "end":   pd.to_numeric(df["End"],   errors="coerce"),
    })
    # simple QC prints
    print(f"{Path(path).name}: n={len(out)}, cov median={np.nanmedian(out['coverage']) :.1f}, "
          f"%mod median={np.nanmedian(out['percent_mod']) :.2f}")
    return out

def main():
    # 1) Bacillus from its two BEDs 
    bc6_path = "bc_hotspots_6mA.bed"
    bc5_path = "bc_hotspots_5mC.bed"
    bc_parts = []
    if Path(bc6_path).exists(): bc_parts.append(load_pacbio_bed(bc6_path, "6mA", "Bacillus"))
    if Path(bc5_path).exists(): bc_parts.append(load_pacbio_bed(bc5_path, "5mC", "Bacillus"))
    assert bc_parts, "Could not find bc_hotspots_6mA.bed / bc_hotspots_5mC.bed"
    bc = pd.concat(bc_parts, ignore_index=True)

    if Path("bc_modifications_with_region.csv").exists():
        reg = pd.read_csv("bc_modifications_with_region.csv")
        if "start" not in reg.columns and "position" in reg.columns:
            reg = reg.rename(columns={"position":"start"})
        reg_key = reg[["chrom","start","region"]].drop_duplicates()
        bc = bc.merge(reg_key, on=["chrom","start"], how="left")

    bc.to_csv("hotspots_hot.csv", index=False)
    print("Wrote corrected Bacillus → hotspots_hot.csv")

    # Pseudomonas 6mA from BED 
    ps6_bed = "psph_hotspots_6mA.bed"
    assert Path(ps6_bed).exists(), f"Missing {ps6_bed}"
    ps6 = load_pacbio_bed(ps6_bed, "6mA", "Pseudomonas")
    # annotating regions 
    ps6.to_csv("psph_6mA_hotspots.csv", index=False)
    print("Wrote corrected Pseudomonas 6mA → psph_6mA_hotspots.csv")

    # Pseudomonas 5mC 
    ps5_path = "psph_5mC_windows_annotated.csv" if Path("psph_5mC_windows_annotated.csv").exists() else "psph_5mC_windows_annot.csv"
    ps5 = pd.read_csv(ps5_path)
    if "contig" in ps5.columns: ps5 = ps5.rename(columns={"contig":"chrom"})
    if "region" not in ps5.columns and Path("psph_5mC_cov30_annot.csv").exists():
        ann = pd.read_csv("psph_5mC_cov30_annot.csv")[["contig","start","end","region"]].rename(columns={"contig":"chrom"})
        ps5 = ps5.merge(ann, on=["chrom","start","end"], how="left")
    ps5 = ps5[["percent_mod","coverage","region"]].copy()
    ps5["species"] = "Pseudomonas"; ps5["modification"] = "5mC"

    # Combining and making the three figures
    comp = pd.concat([
        bc[["percent_mod","coverage","region","modification","species"]],
        ps5[["percent_mod","coverage","region","modification","species"]],
        ps6[["percent_mod","coverage","modification","species"]].assign(region=np.nan)
    ], ignore_index=True)

    # Ensuring % in 0–100
    if comp["percent_mod"].dropna().quantile(0.9) <= 1.0:
        comp["percent_mod"] *= 100.0

    plt.rcParams.update({"figure.figsize":(4.6,3.4),"figure.dpi":300,"savefig.dpi":600,"font.size":9})

    # Fig1
    MIN_PCT, MIN_COV = 5.0, 30
    hot = comp[(comp["coverage"]>=MIN_COV) & (comp["percent_mod"]>=MIN_PCT)]
    cnt = hot.groupby(["species","modification"]).size().rename("count").reset_index()
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

    # Fig2 percent-mod distributions 
    groups, labels = [], []
    for sp in sorted(comp["species"].unique()):
        for mod in sorted(comp.loc[comp["species"]==sp,"modification"].unique()):
            vals = comp[(comp["species"]==sp)&(comp["modification"]==mod)]["percent_mod"].dropna().values
            if len(vals): groups.append(vals); labels.append(f"{sp} {mod} (n={len(vals)})")
    if groups:
        plt.figure(); plt.boxplot(groups, labels=labels, showmeans=True)
        plt.ylabel("Percent modified"); plt.title("Percent Modified per Species × Modification")
        plt.xticks(rotation=15); plt.ylim(0,100); plt.tight_layout(); plt.savefig("Cross_Fig2_percent_mod_box.png", bbox_inches="tight"); plt.close()

    # Fig3 region fractions 
    if "region" in comp.columns and comp["region"].notna().any():
        reg = (comp.dropna(subset=["region"]).groupby(["species","modification","region"]).size()
               .rename("count").reset_index())
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

    # sanity check
    print("\nNew hotspot counts (≥5%, ≥30×):")
    print(hot.groupby(['species','modification']).size().rename('n'))
    print("\nNew percent-mod medians:")
    print(comp.groupby(['species','modification'])['percent_mod'].median())

if __name__ == "__main__":
    main()
