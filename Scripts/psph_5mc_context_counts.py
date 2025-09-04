#!/usr/bin/env python3
"""
# Counting immediate sequence context around centered 5mC windows.
# Input: psph_5mC_cov30_pm10.fa
"""

# Counting left/right nearby base frequencies 
from collections import Counter

right1, left1 = Counter(), Counter()
with open("psph_5mC_cov30_pm10.fa") as f:
    for line in f:
        if line.startswith(">"): 
            seq = next(f).strip().upper()
            center = len(seq)//2
            left1[seq[center-1]]  += 1
            right1[seq[center+1]] += 1

print("Left of C (−1):", left1)
print("Right of C (+1):", right1)
