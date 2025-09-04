#!/usr/bin/env python3
"""
# Visualize Pseudomonas 5mC motif 
# Input
# psph_5mC_cov30_PFM.csv   
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt

def main():
    pfm = pd.read_csv("psph_5mC_cov30_PFM.csv", index_col=0)[["A","C","G","T"]]
    L = len(pfm)
    center = L//2                      
    pwm = pfm.div(pfm.sum(1), axis=0)  # probabilities per position
    pwm.index = np.arange(-center, L-center)  

    plt.figure(figsize=(10,3))
    bottom = np.zeros(L)
    for base, color in zip(["A","C","G","T"], ["#1f77b4","#2ca02c","#ff7f0e","#d62728"]):
        vals = pwm[base].values
        plt.bar(pwm.index, vals, bottom=bottom, label=base, color=color, width=0.9)
        bottom += vals

    plt.axvline(0, color="k", ls="--", lw=1.5)
    plt.text(0.1, 1.02, "[mC]", transform=plt.gca().get_xaxis_transform())
    plt.ylim(0,1); plt.xlim(-center-0.5, center+0.5)
    plt.xlabel("Position relative to methylated C")
    plt.ylabel("Base probability")
    plt.title("Pseudomonas 5mC motif (nsites=17) — consensus ≈ CC[mC]GAG")
    plt.legend(ncol=4, frameon=False)
    plt.tight_layout()
    plt.savefig("psph_5mC_stackedbars.png", dpi=300)
    print("Saved: psph_5mC_stackedbars.png")

if __name__ == "__main__":
    main()
