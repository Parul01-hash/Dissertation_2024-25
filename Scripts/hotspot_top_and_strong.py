#!/usr/bin/env python3
"""
# Top-25 & strong hotspots 
"""

# Top-25 high-confidence hotspots (cov ≥30)
top = (
    df.query("valid_cov >= 30")
      .sort_values(["frac_mod", "valid_cov"], ascending=False)
      .head(25)
)
top[[
    "chrom","position","strand","valid_cov","n_mod",
    "frac_mod","percent_mod","ref_base","expected_base"
]]
top.to_csv("top25_hotspots.csv", index=False)

# Strong hotspots
strong = (
    df.query("valid_cov >= 50 and frac_mod >= 0.20")
      .sort_values(["frac_mod", "valid_cov"], ascending=False)
)
len(strong), strong.head(10)

# Save strong hotspots 
strong.to_csv("strong_hotspots_cov50_frac20.csv", index=False)
