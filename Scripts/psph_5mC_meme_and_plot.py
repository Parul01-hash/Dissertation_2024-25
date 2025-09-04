#!/usr/bin/env python3
"""
# Building PWM,consensus, MEME file and quick plot from PFM.
# Loading a PFM CSV, 
# converting to a PWM, 
# printing a naïve consensus

# Input
# psph_5mC_cov30_PFM.csv 

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PFM_CSV = "psph_5mC_cov30_PFM.csv"
bases = ["A","C","G","T"]

# loading PFM 
pfm_df = pd.read_csv(PFM_CSV, index_col=0)
# ensuring column order A,C,G,T
pfm_df = pfm_df.reindex(columns=bases)

L = len(pfm_df)
center = L // 2
print(f"Loaded PFM: length={L}, center index={center}")

# PWM 
row_sums = pfm_df.sum(axis=1).replace(0, np.nan)
pwm = (pfm_df.div(row_sums, axis=0)).fillna(0)

# naive consensus
cons = "".join(bases[int(np.argmax(row))] for row in pwm.to_numpy())
print("Naive consensus:", cons)

# writing a MEME motif file 
with open("psph_5mC.meme","w") as out:
    out.write("MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\n")
    out.write(f"MOTIF psph_5mC_core\nletter-probability matrix: alength= 4 w= {L} nsites= {int(row_sums.iloc[0]) if pd.notna(row_sums.iloc[0]) else 0} E= 0\n")
    for i in range(L):
        a,c,g,t = pwm.iloc[i][bases].tolist()
        out.write(f" {a:.6f} {c:.6f} {g:.6f} {t:.6f}\n")
print("Wrote MEME motif → psph_5mC.meme")

# quick visualization
plt.figure(figsize=(10,3))
plt.plot(pfm_df["A"], label="A")
plt.plot(pfm_df["C"], label="C")
plt.plot(pfm_df["G"], label="G")
plt.plot(pfm_df["T"], label="T")
plt.axvline(center, linestyle="--")
plt.title("Base counts by position (center = methylated C)")
plt.xlabel("position (0-based)")
plt.ylabel("count")
plt.legend()
plt.tight_layout()
plt.savefig("psph_5mC_pf_counts.png", dpi=150)
print("Saved plot → psph_5mC_pf_counts.png")
