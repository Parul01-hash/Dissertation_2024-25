#!/usr/bin/env python3
"""

# Input
# region_plot_input.csv
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    
    df = pd.read_csv("region_plot_input.csv")

    
    df["region"] = df["region"].replace({"gene_body": "CDS"})

    # species and modification
    df["group"] = list(zip(df["species"], df["modification"]))
    counts = (df.groupby(["group","region"])
                .size().rename("count").reset_index())

    piv = counts.pivot(index="group", columns="region", values="count").fillna(0)
    fractions = piv.div(piv.sum(axis=1), axis=0)
 
    fractions.to_csv("region_fractions_CDSonly.csv")
    print(fractions)

    # plotting (CDS, intergenic, promoter)
    order = [("Bacillus","5mC"),("Bacillus","6mA"),("Pseudomonas","5mC"),("Pseudomonas","6mA")]
    fractions = fractions.reindex(order)

    plt.rcParams.update({"figure.figsize":(5.4,3.6), "figure.dpi":300, "savefig.dpi":600, "font.size":9})
    x = np.arange(len(fractions.index))
    bottom = np.zeros(len(x))

    for col in ["CDS", "intergenic", "promoter"]:
        if col in fractions.columns:
            vals = fractions[col].values
            plt.bar(x, vals, bottom=bottom, label=col)
            bottom += vals

    plt.xticks(x, [f"{a}\n{b}" for a,b in fractions.index])
    plt.ylabel("Fraction of sites")
    plt.title("Region Distribution per Species × Modification (CDS-only)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("Cross_Fig3_region_stack_CDSonly.png", bbox_inches="tight")
    plt.savefig("Cross_Fig3_region_stack_CDSonly.pdf", bbox_inches="tight")
    print("Saved: Cross_Fig3_region_stack_CDSonly.[png|pdf] and region_fractions_CDSonly.csv")

if __name__ == "__main__":
    main()
