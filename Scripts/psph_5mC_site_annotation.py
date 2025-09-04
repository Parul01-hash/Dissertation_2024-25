#!/usr/bin/env python3
"""
Annotation 5mC sites to GFF features (CDS/gene/promoter)

Inputs 
# psph_5mC_cov30.bed      
# psph_annotation.gff     
# psph_reference.fasta    
"""

import re
import pandas as pd
from pyfaidx import Fasta

# paths
BED = "psph_5mC_cov30.bed"          
GFF = "psph_annotation.gff"
FASTA = "psph_reference.fasta"

# loading sites  
sites = pd.read_csv(BED, sep="\t", header=None, names=["contig","start","end"])
sites["pos0"] = sites["start"].astype(int)     
sites["pos1"] = sites["pos0"] + 1              

# parse GFF 
rows = []
with open(GFF) as fh:
    for ln in fh:
        if not ln.strip() or ln.startswith("#"): 
            continue
        seqid, source, ftype, start, end, score, strand, phase, attrs = ln.rstrip("\n").split("\t")
        if ftype not in {"gene","CDS","rRNA","tRNA","ncRNA"}:
            continue
        
        s1, e1 = int(start), int(end)
        s0, e0 = s1-1, e1
        
        def pick(*keys):
            for k in keys:
                m = re.search(rf"{k}=([^;]+)", attrs)
                if m: return m.group(1)
            return None
        locus = pick("locus_tag","ID","gene")
        name  = pick("gene","Name","product") or locus
        prod  = pick("product","note")
        rows.append([seqid, ftype, s0, e0, s1, e1, strand, locus, name, prod])
gff = pd.DataFrame(rows, columns=["seqid","type","start0","end0","start1","end1","strand","locus","name","product"])

assert len(gff), "No features parsed from GFF. Check the path/format."

# contig name mapping
ref = Fasta(FASTA)
ref_ids = set(ref.keys())
contig_map = {}
for c in sites["contig"].unique():
    if c in ref_ids or c in gff["seqid"].unique():
        contig_map[c] = c
    elif c.startswith("CP000058") and "NC_005773.3" in ref_ids:
        contig_map[c] = "NC_005773.3"
    elif c.startswith("CP000060") and "NC_007274.1" in ref_ids:
        contig_map[c] = "NC_007274.1"
    else:
        contig_map[c] = c
sites["seqid"] = sites["contig"].map(contig_map)

# building promoter intervals per gene 
PROM_UP = 250   
PROM_DN = 50   

genes = gff[gff["type"]=="gene"].copy()
prom_rows = []
for _, g in genes.iterrows():
    if g["strand"] == "+":
        s0 = max(0, g["start0"] - PROM_UP)
        e0 = g["start0"] + PROM_DN
    else:  # '-' strand
        s0 = g["end0"] - PROM_DN
        e0 = g["end0"] + PROM_UP
    prom_rows.append([g["seqid"], s0, e0, g["strand"], g["locus"], g["name"], g["product"]])
prom = pd.DataFrame(prom_rows, columns=["seqid","start0","end0","strand","locus","name","product"])
prom["type"] = "promoter"

# finding overlapping feature for a single site
def overlap_row(sub, pos0):
    hit = sub[(sub["start0"] <= pos0) & (pos0 < sub["end0"])]
    return hit.iloc[0] if len(hit) else None

ann = []
for _, s in sites.iterrows():
    seqid, p0 = s["seqid"], int(s["pos0"])
    sub_all = gff[gff["seqid"]==seqid]
    sub_prom = prom[prom["seqid"]==seqid]

    hit = overlap_row(sub_all[sub_all["type"].isin(["CDS","rRNA","tRNA","ncRNA"])], p0)
    region = None
    if hit is not None:
        region = hit["type"]
    else:
        phit = overlap_row(sub_prom, p0)
        if phit is not None:
            region = "promoter"; hit = phit
        else:
            hit = overlap_row(sub_all[sub_all["type"]=="gene"], p0)
            region = "gene" if hit is not None else "intergenic"

    locus = hit["locus"] if hit is not None else None
    name  = hit["name"]  if hit is not None else None
    prod  = hit["product"] if hit is not None else None
    strand= hit["strand"] if hit is not None else None
    ann.append([s["contig"], s["start"], s["end"], seqid, p0, region, locus, name, prod, strand])

annot = pd.DataFrame(ann, columns=["contig","start","end","seqid","pos0","region","locus","gene_name","product","gene_strand"])

annot.to_csv("psph_5mC_cov30_annot.csv", index=False)
print("Saved per-site annotations → psph_5mC_cov30_annot.csv")
print(annot.head(12).to_string(index=False))

# summarising by region & gene  
by_region = annot["region"].value_counts().rename_axis("region").reset_index(name="n_sites")
by_gene = (annot.groupby(["locus","gene_name","product"])
                 .size().reset_index(name="n_sites")
                 .sort_values("n_sites", ascending=False))

print("\nSites by region:")
print(by_region.to_string(index=False))

print("\nTop genes by # 5mC sites:")
print(by_gene.head(15).to_string(index=False))
 
KEYWORDS = {
    "Regulation": ["transcription", "sigma", "regulator", "two-component", "response regulator", "sensor"],
    "Virulence":  ["secretion", "type iii", "type vi", "pilus", "flagell", "motility", "exopolysaccharide", "toxin"],
    "DNA/RM":     ["restriction", "methyltransferase", "endonuclease", "DNA repair", "recA", "uvr", "mut"],
    "Metabolism": ["dehydrogenase", "kinase", "synthase", "transferase", "transport", "permease"],
}
def score_product(prod):
    if not isinstance(prod, str): return []
    prod_l = prod.lower()
    hits = [k for k, words in KEYWORDS.items() if any(w in prod_l for w in words)]
    return hits

by_gene["tags"] = by_gene["product"].apply(score_product)
interesting = by_gene[by_gene["tags"].map(len) > 0].copy()
print("\nGenes with 'interesting' tags (keyword hit in product):")
print(interesting.head(20).to_string(index=False))
interesting.to_csv("psph_5mC_top_genes_tagged.csv", index=False)
print("\nSaved tagged gene list → psph_5mC_top_genes_tagged.csv")
