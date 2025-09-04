#!/usr/bin/env python3
"""
# Pseudomonas 5mC motif enrichment counters.
Input: psph_5mC_cov30_pm10.fa  
"""

# Counting GAG(+1..+3) and left-neighbor C(−1) around centered 5mC sites to quantify motif enrichment.
fa = "psph_5mC_cov30_pm10.fa"

# loading sequences
seqs = []
with open(fa) as f:
    for line in f:
        if line.startswith(">"):
            s = next(f).strip().upper()
            seqs.append(s)

L = len(seqs[0]); c = L // 2  

core_GAG = sum(s[c+1:c+4] == "GAG" for s in seqs)
left_C   = sum(s[c-1] == "C" for s in seqs)

print(f"Total windows: {len(seqs)}")
print(f"GAG at +1..+3: {core_GAG}/{len(seqs)}")
print(f"C at -1:       {left_C}/{len(seqs)}")

# joint motif enrichment.
both = sum((s[c-1]=="C") and (s[c+1:c+4]=="GAG") for s in seqs)
print(f"C at -1 AND GAG at +1..+3: {both}/{len(seqs)}")
