#!/usr/bin/env python3
"""
# Cross-species, all-mods figures with Pseudomonas 6mA region annotation

# inputs 
# hotspots_hot.csv
# psph_6mA_hotspots.csv
# psph_5mC_hotspots_full.csv  OR  psph_5mC_windows_annotated.csv 
# GFF: psph_annotation.gff
"""

import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def find_gff():
    
    for name in ["psph_annotation.gff", "psph.gff"]:
        if Path(name).exists(): return name
    gffs = sorted([p for p in Path(".").glob("*.gff")])
    return str(gffs[0]) if gffs else None

def parse_gff_rows(gff_path, keep_types=("gene","CDS")):
    rows = []
    with open(gff_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype not in keep_types: continue
            try:
                start = int(start); end = int(end)
            except:
                continue
            # parse attributes
            a = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k,v = kv.split("=",1)
                    a[k.strip()] = v.strip()
            locus = a.get("locus_tag") or a.get("ID") or a.get("Name")
            gene  = a.get("gene") or a.get("Name")
            prod  = a.get("product")
            rows.append((seqid, ftype, start, end, strand, locus, gene, prod))
    gff = pd.DataFrame(rows, columns=["seqid","type","start","end","strand","locus","gene","product"])
    return gff

def annotate_ps6ma_regions(ps6_df, gff_path, min_pct=5.0, min_cov=30, promoter_bp=200):
    """Annotate ONLY the hotspots (≥min_pct & ≥min_cov) to speed up."""
    use = ps6_df[(ps6_df["coverage"]>=min_cov) & (ps6_df["percent_mod"]>=min_pct)].copy()
    if use.empty:
        print("No PS-6mA hotspots pass the filter; skipping region annotation.")
        use["region"] = np.nan
        return use

    
    if not {"chrom","start","end"} <= set(use.columns):
    
        print("PS-6mA table lacks chrom/start/end; cannot annotate. Keeping region=NaN.")
        use["region"] = np.nan
        return use

    gff = parse_gff_rows(gff_path, keep_types=("gene","CDS"))
    # Split by type
    cds = gff[gff["type"]=="CDS"].copy()
    genes = gff[gff["type"]=="gene"].copy()

    # building promoter table from genes
    prom_rows = []
    for _,r in genes.iterrows():
        if r["strand"] == "+":
            p_start = max(1, r["start"] - promoter_bp)
            p_end   = r["start"] - 1
        else:
            p_start = r["end"] + 1
            p_end   = r["end"] + promoter_bp
        prom_rows.append((r["seqid"], p_start, p_end, r["strand"], r["locus"], r["gene"], r["product"]))
    prom = pd.DataFrame(prom_rows, columns=["seqid","start","end","strand","locus","gene","product"])

    # function to classify a single site 
    def classify_row(ch, pos):
        hit = cds[(cds["seqid"]==ch) & (cds["start"]<=pos) & (cds["end"]>=pos)]
        if len(hit):
            h = hit.iloc[0]
            return "CDS", h["locus"], h["gene"], h["product"], h["strand"]
        # within gene 
        ghit = genes[(genes["seqid"]==ch) & (genes["start"]<=pos) & (genes["end"]>=pos)]
        if len(ghit):
            g = ghit.iloc[0]
            return "gene", g["locus"], g["gene"], g["product"], g["strand"]
        # promoter
        phit = prom[(prom["seqid"]==ch) & (prom["start"]<=pos) & (prom["end"]>=pos)]
        if len(phit):
            p = phit.iloc[0]
            return "promoter", p["locus"], p["gene"], p["product"], p["strand"]
        return "intergenic", None, None, None, None

    # annotation
    out_cols = ["region","locus","gene_name","product","gene_strand"]
    for col in out_cols: use[col] = None
    # process per chromosome
    for ch, dfc in use.groupby("chrom"):
        mask_cds   = (cds["seqid"]==ch)
        mask_gene  = (genes["seqid"]==ch)
        mask_prom  = (prom["seqid"]==ch)
        
        if not (mask_cds.any() or mask_gene.any() or mask_prom.any()):
            continue
        idxs = dfc.index
        for i in idxs:
            pos = int(use.at[i,"start"])
            region, locus, gname, prod, gstr = classify_row(ch, pos)
            use.at[i,"region"] = region
            use.at[i,"locus"] = locus
            use.at[i,"gene_name"] = gname
            use.at[i,"product"] = prod
            use.at[i,"gene_strand"] = gstr
    return use

def main():
    #  loading inputs
    assert Path("hotspots_hot.csv").exists(), "Missing Bacillus table: hotspots_hot.csv"
    assert Path("psph_6mA_hotspots.csv").exists(), "Missing PS 6mA table: psph_6mA_hotspots.csv"

    bc  = pd.read_csv("hotspots_hot.csv")  # Bacillus 
    ps6 = pd.read_csv("psph_6mA_hotspots.csv")  # Pseudomonas 6mA from BED
    # ensuring numeric
    for c in ["percent_mod","coverage","start","end"]:
        if c in bc.columns:  bc[c]  = pd.to_numeric(bc[c],  errors="coerce")
        if c in ps6.columns: ps6[c] = pd.to_numeric(ps6[c], errors="coerce")

    # Pseudomonas 5mC
    if Path("psph_5mC_hotspots_full.csv").exists():
        ps5 = pd.read_csv("psph_5mC_hotspots_full.csv")
        if "percent_mod" not in ps5.columns and {"NumModified","Coverage"} <= set(ps5.columns):
            ps5["percent_mod"] = ps5["NumModified"]/ps5["Coverage"]*100.0
            ps5.rename(columns={"Coverage":"coverage"}, inplace=True)
        keep = ["percent_mod","coverage","chrom","start","end"]
        have = [k for k in keep if k in ps5.columns]
        ps5 = ps5[have].copy()
    else:
        win = "psph_5mC_windows_annotated.csv" if Path("psph_5mC_windows_annotated.csv").exists() else \
              ("psph_5mC_windows_annot.csv" if Path("psph_5mC_windows_annot.csv").exists() else None)
        assert win is not None, "No Pseudomonas 5mC file found"
        ps5 = pd.read_csv(win)
        if "contig" in ps5.columns: ps5.rename(columns={"contig":"chrom"}, inplace=True)
        # try to standardize numeric cols
        if "percent_mod" not in ps5.columns and "mod_frac" in ps5.columns:
            ps5["percent_mod"] = pd.to_numeric(ps5["mod_frac"], errors="coerce")*100.0
        if "coverage" not in ps5.columns and "cov" in ps5.columns:
            ps5.rename(columns={"cov":"coverage"}, inplace=True)

    # annotating PS-6mA hotspots via GFF 
    gff_path = find_gff()
    if gff_path:
        ps6_hot_annot = annotate_ps6ma_regions(ps6, gff_path, min_pct=5.0, min_cov=30, promoter_bp=200)
        
        merge_cols = ["chrom","start","end","percent_mod","coverage"]
        common = [c for c in merge_cols if c in ps6.columns and c in ps6_hot_annot.columns]
        ps6 = ps6.merge(ps6_hot_annot[common + ["region","locus","gene_name","product","gene_strand"]],
                        on=common, how="left")
        ps6.to_csv("psph_6mA_hotspots_annot.csv", index=False)
        print(f"Saved PS 6mA hotspot annotations → psph_6mA_hotspots_annot.csv (rows annotated={ps6_hot_annot.shape[0]})")
    else:
        print("No GFF found. Skipping PS-6mA region annotation (region will be NaN).")

    # building combined table 
    bc_use = bc[["percent_mod","coverage","region","modification"]].copy()
    bc_use["species"] = "Bacillus"

    ps6_use = ps6.copy()
    ps6_use["species"] = "Pseudomonas"
    ps6_use["modification"] = "6mA"
    if "region" not in ps6_use.columns: ps6_use["region"] = np.nan

    ps5_use = ps5.copy()
    ps5_use["species"] = "Pseudomonas"
    ps5_use["modification"] = "5mC"
    if "region" not in ps5_use.columns: ps5_use["region"] = np.nan

    # keeping common columns
    def keep_cols(df):
        cols = ["species","modification","percent_mod","coverage","region"]
        return df[[c for c in cols if c in df.columns]].copy()

    comp = pd.concat([keep_cols(bc_use), keep_cols(ps6_use), keep_cols(ps5_use)], ignore_index=True)

   
    if comp["percent_mod"].dropna().quantile(0.9) <= 1.0:
        comp["percent_mod"] = comp["percent_mod"] * 100.0

    comp.to_csv("Cross_species_mod_summary.csv", index=False)

    #  plots
    plt.rcParams.update({
        "figure.figsize": (4.6, 3.4),
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8
    })

    # Fig 1: counts 
    MIN_PCT, MIN_COV = 5.0, 30
    hot = comp[(comp["coverage"]>=MIN_COV) & (comp["percent_mod"]>=MIN_PCT)]
    counts = (hot.groupby(["species","modification"]).size().rename("count").reset_index())
    if len(counts):
        piv = counts.pivot(index="species", columns="modification", values="count").fillna(0)
        x = np.arange(len(piv.index)); mods = list(piv.columns); w = 0.8/max(1,len(mods))
        plt.figure()
        for i,m in enumerate(mods):
            vals = piv[m].values
            bars = plt.bar(x+i*w-(len(mods)-1)*w/2, vals, width=w, label=m)
            for b,v in zip(bars, vals):
                plt.text(b.get_x()+b.get_width()/2, v, f"{int(v)}", ha="center", va="bottom", fontsize=8)
        plt.xticks(x, piv.index); plt.ylabel(f"Hotspot count (≥{MIN_PCT:.0f}%, ≥{MIN_COV}×)")
        plt.title("Hotspots per Species × Modification")
        if len(mods)>1: plt.legend()
        plt.tight_layout()
        plt.savefig("Cross_Fig1_counts_relaxed.png", bbox_inches="tight")
        plt.savefig("Cross_Fig1_counts_relaxed.pdf", bbox_inches="tight")
        plt.close()

    # Fig 2 percent-mod boxplot 
    groups, labels = [], []
    for sp in sorted(comp["species"].dropna().unique()):
        for mod in sorted(comp.loc[comp["species"]==sp, "modification"].dropna().unique()):
            vals = comp[(comp["species"]==sp)&(comp["modification"]==mod)]["percent_mod"].dropna().values
            if len(vals):
                groups.append(vals); labels.append(f"{sp} {mod} (n={len(vals)})")
    if groups:
        plt.figure()
        plt.boxplot(groups, labels=labels, showmeans=True)
        plt.ylabel("Percent modified")
        plt.title("Percent Modified per Species × Modification")
        plt.xticks(rotation=15); plt.ylim(0,100)
        plt.tight_layout()
        plt.savefig("Cross_Fig2_percent_mod_box.png", bbox_inches="tight")
        plt.savefig("Cross_Fig2_percent_mod_box.pdf", bbox_inches="tight")
        plt.close()

    # Fig 3 region distribution 
    fig3_saved = False
    if "region" in comp.columns and comp["region"].notna().any():
        reg = (comp.dropna(subset=["region"])
               .groupby(["species","modification","region"]).size()
               .rename("count").reset_index())
        if len(reg):
            reg_piv = reg.pivot_table(index=["species","modification"], columns="region", values="count", fill_value=0)
            frac = reg_piv.div(reg_piv.sum(axis=1), axis=0)
            idx = np.arange(len(frac.index)); bottom = np.zeros(len(idx))
            plt.figure()
            for col in frac.columns:
                vals = frac[col].values
                plt.bar(idx, vals, bottom=bottom, label=col); bottom += vals
            plt.xticks(idx, [f"{a}\n{b}" for a,b in frac.index])
            plt.ylabel("Fraction of sites")
            plt.title("Region Distribution per Species × Modification")
            plt.legend(); plt.tight_layout()
            plt.savefig("Cross_Fig3_region_stack_small.png", bbox_inches="tight")
            plt.savefig("Cross_Fig3_region_stack_small.pdf", bbox_inches="tight")
            plt.close(); fig3_saved = True

    print("\n=== Figure-ready numbers ===")
    print("Counts (≥5%, ≥30×):")
    print(counts.sort_values(['species','modification']).to_string(index=False))
    print("\nPercent-mod medians:")
    print(comp.groupby(['species','modification'])['percent_mod'].median())
    print("\nSaved:")
    print("  Cross_Fig1_counts_relaxed.[png|pdf]")
    print("  Cross_Fig2_percent_mod_box.[png|pdf]")
    print("  Cross_Fig3_region_stack_small.[png|pdf]" + ("" if fig3_saved else " (skipped: no region labels)"))
    print("  Cross_species_mod_summary.csv")

if __name__ == "__main__":
    main()
