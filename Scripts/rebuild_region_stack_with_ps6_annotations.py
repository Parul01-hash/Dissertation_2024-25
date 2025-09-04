#!/usr/bin/env python3
"""
# Rebuilding the region-distribution figure 
# Inputs 
# Cross_species_mod_summary.csv
# psph_6mA_hotspots_annot.csv
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Loading the combined summary 
comp = pd.read_csv("Cross_species_mod_summary.csv")

# Loading PS 6mA hotspot annotations and merge by chrom/start/end only
ps6_ann = pd.read_csv("psph_6mA_hotspots_annot.csv")
keys = [k for k in ["chrom","start","end"] if k in ps6_ann.columns and k in comp.columns]
if keys:
    comp = comp.merge(
        ps6_ann[keys + ["region"]].dropna(subset=keys).drop_duplicates(keys),
        on=keys, how="left", suffixes=("", "_ps6")
    )
    
    comp["region"] = comp["region"].fillna(comp.get("region_ps6"))
    if "region_ps6" in comp.columns: comp = comp.drop(columns=["region_ps6"])

# ploting
plt.rcParams.update({"figure.figsize":(4.6,3.4),"figure.dpi":300,"savefig.dpi":600,"font.size":9})
has_region = comp["region"].notna()
reg = (comp[has_region]
       .groupby(["species","modification","region"])
       .size().rename("count").reset_index())
if len(reg):
    piv = reg.pivot_table(index=["species","modification"], columns="region",
                          values="count", fill_value=0)
    frac = piv.div(piv.sum(axis=1), axis=0)
    import numpy as np
    x = np.arange(len(frac.index)); bottom = np.zeros(len(x))
    plt.figure()
    for col in frac.columns:
        vals = frac[col].values
        plt.bar(x, vals, bottom=bottom, label=col)
        bottom += vals
    plt.xticks(x, [f"{a}\n{b}" for a,b in frac.index])
    plt.ylabel("Fraction of sites"); plt.title("Region Distribution per Species × Modification")
    plt.legend(); plt.tight_layout()
    plt.savefig("Cross_Fig3_region_stack_small.png", bbox_inches="tight")
    plt.savefig("Cross_Fig3_region_stack_small.pdf", bbox_inches="tight")
    plt.close()
    print("Rewrote Cross_Fig3_region_stack_small.[png|pdf] with Pseudomonas included (if annotated).")
else:
    print("No rows with region labels to plot.")

import pandas as pd

comp = pd.read_csv(
    "Cross_species_mod_summary.csv",
    low_memory=False,
    dtype={
        "percent_mod": "float64",
        "coverage": "float64",
        "species": "category",
        "modification": "category",
        "region": "category"
    }
)
