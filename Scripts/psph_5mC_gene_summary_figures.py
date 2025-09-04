#!/usr/bin/env python3
"""
# Generate compact P. syringae 5mC figures  and a small gene summary.
# Inputs 
# psph_5mC_cov30_annot.csv           
# psph_5mC_windows_annotated.csv     
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ann = pd.read_csv("psph_5mC_cov30_annot.csv")
win = pd.read_csv("psph_5mC_windows_annotated.csv") if os.path.exists("psph_5mC_windows_annotated.csv") \
      else pd.read_csv("psph_5mC_windows_annot.csv")

# Merging to attach percent_mod & coverage to region annotations
merge_cols = ["contig","start","end"]
df = ann.merge(win[merge_cols + ["coverage","percent_mod"]], on=merge_cols, how="left")

plt.rcParams.update({
    "figure.figsize": (4.2, 3.2),   
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
})

# Fig 1: Sites by genomic region (bar)
region_counts = df["region"].value_counts().rename_axis("Region").reset_index(name="Count")
fig1_png, fig1_pdf = "Fig1_5mC_by_region_small.png", "Fig1_5mC_by_region_small.pdf"

plt.figure()
plt.bar(region_counts["Region"], region_counts["Count"])
plt.title("5mC Sites by Region")
plt.xlabel("Region"); plt.ylabel("Count")
plt.gca().margins(x=0.1)
plt.tight_layout()
plt.savefig(fig1_png, bbox_inches="tight")
plt.savefig(fig1_pdf, bbox_inches="tight")
plt.close()

# Fig 2: Percent modified by region (boxplot)
order = region_counts["Region"].tolist()
groups = [df.loc[df["region"]==r, "percent_mod"].dropna().values for r in order]
fig2_png, fig2_pdf = "Fig2_5mC_percent_by_region_small.png", "Fig2_5mC_percent_by_region_small.pdf"

plt.figure()
bp = plt.boxplot(groups, labels=order, showmeans=True)
plt.title("Percent Modified by Region")
plt.ylabel("Percent modified"); plt.xlabel("Region")
plt.gca().set_ylim(bottom=0)        # start at 0%
plt.tight_layout()
plt.savefig(fig2_png, bbox_inches="tight")
plt.savefig(fig2_pdf, bbox_inches="tight")
plt.close()

# Fig 3: Functional categories (horizontal bar)
KEYWORDS = {
    "Regulation": ["transcription", "sigma", "regulator", "two-component", "response regulator", "sensor"],
    "Transport":  ["transporter", "permease", "secretin", "mfs", "abc transporter", "phosphonate", "phnD"],
    "Secretion":  ["secretion", "secretin", "type iii", "type vi"],
    "DNA/mobile": ["integrase", "recombinase", "restriction", "methyltransferase", "endonuclease"],
    "Metabolism": ["dehydrogenase", "kinase", "synthase", "transferase", "reductase", "aminoglycoside"],
}
def tag_product(prod):
    if not isinstance(prod, str): return []
    p = prod.lower()
    return [k for k, words in KEYWORDS.items() if any(w in p for w in words)]

gene_df = df.dropna(subset=["locus"]).drop_duplicates(subset=["locus"]).copy()
gene_df["tags"] = gene_df["product"].apply(tag_product)

from collections import Counter
tc = Counter()
for tags in gene_df["tags"]:
    tc.update(tags)

tag_df = (pd.DataFrame({"Tag": list(tc.keys()), "Count": list(tc.values())})
          .sort_values("Count", ascending=False))

fig3_png, fig3_pdf = "Fig3_5mC_functional_cats_small.png", "Fig3_5mC_functional_cats_small.pdf"
plt.figure()
plt.barh(tag_df["Tag"], tag_df["Count"])
plt.title("Functional Categories of Genes with 5mC")
plt.xlabel("Gene count")
plt.tight_layout()
plt.savefig(fig3_png, bbox_inches="tight")
plt.savefig(fig3_pdf, bbox_inches="tight")
plt.close()

print("Saved compact figures:")
print(" -", fig1_png, "|", fig1_pdf)
print(" -", fig2_png, "|", fig2_pdf)
print(" -", fig3_png, "|", fig3_pdf)

# creating table 
summary = (df.sort_values("percent_mod", ascending=False)
             .drop_duplicates(subset=["locus"])
             [["locus","gene_name","product","region","percent_mod","coverage","contig","start","end"]])
summary.to_csv("Table_S1_psph_5mC_gene_summary.csv", index=False)
print(" - Table_S1_psph_5mC_gene_summary.csv")
