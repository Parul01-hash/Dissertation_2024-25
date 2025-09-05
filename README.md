# Dissertation
Detection and Comparative Analysis of 5mC and 6mA Modifications in Bacillus cereus and Pseudomonas

# REQUIREMENTS
  - Python 3.9+
  - pip install: pandas, numpy, matplotlib, pyfaidx
# Data
  # Raw
  - bc_hotspots_6mA.bed
  - bc_hotspots_5mC.bed
  - psph_hotspots_6mA.bed
  - psph_hotspots_5mC.bed
  # Annotation
  - bc_annotation.gff
  - psph_annotation.gff
  # References
  - bc_reference.fasta
  - bc_reference.fasta.fai
  - psph_reference.fasta
  - psph_reference.fasta.fai
# Analysis
  
# STEP 1 — Hotspot summary
  - Run: Scripts/Hotspot_summary.py
  - Required:
     - df_6ma_std (Bacillus 6mA)
     - df_5mc_std (Bacillus 5mC)
     - df_pse_6ma_std (Pseudomonas 6mA)
     - df_pse_5mc_std (Pseudomonas 5mC)
  - Output: cross_species_global_summary_hotspots.tsv
# STEP 2 — Loading genome wide summary
  - Quick load/preview of the global summary table for plots/QC.
  - Run: Scripts/genome_wide_summary.py
  - Required: cross_species_global_summary.tsv
# STEP 3 — Bacillus BED → ref-base QC
  - Read BED18, map contigs to FASTA, fetch reference base, compute strand-aware mismatch.
  - Run: Scripts/bc_modkit_refbase_qc.py
  - Required:
      - Data/raw/bc_hotspots_6mA.bed
      - Data/raw/bc_hotspots_5mC.bed
      - Data/references/bc_reference.fasta
  - Output: hotspots_clean.csv, hotspots_preview.csv, give  mismatch rate
# Step 4 - Hotspot table
- Run: Scripts/clean_hotspot_table.py
- Required: bc_reference.fasta
- Output: hotspots_clean.csv , mismatch rate, preview
# STEP 5 — Top & Strong hotspots
  - Export top-25 (cov ≥30) and strong (cov ≥50 & frac ≥0.20) hotspot tables.
  - Run: Scripts/hotspot_top_and_strong.py
  - Required: hotspots_clean.csv 
  - Output: top25_hotspots.csv, strong_hotspots_cov50_frac20.csv
# Step 6 - Pseudomonas top-20 5mC sites
- computing frac_mod/n_mod, creating CSV and IGV BED
- Run: Scripts/psph_5mC_top20.py
- Required: psph_top20_5mC_sites.tsv
- Output : psph_top20_5mC_sites_fixed.csv, psph_top20_5mC_sites.bed
# Step 7 - Pseudomonas 5mC ±10bp motif windows
- Run: Scripts/psph_5mc_cov30_prepare_motif_inputs.py
- Required: psph_top20_5mC_sites_fixed.csv, psph_reference.fasta
- Output: psph_5mC_cov30.bed, psph_5mC_cov30_pm10.fa
# Step 8 - Counting immediate sequence context around centered 5mC
- Reads oriented ±10 bp FASTA windows and counts bases at -1  and +1.
- Run: Scripts/psph_5mc_context_counts.py
- Required: psph_5mC_cov30_pm10.fa
- Output: prints two Counter dicts for left and right bases
# Step 9 - Counting nucleotide frequencies
- Counts base frequencies at positions −3..+3 around centered 5mC windows.
- Run: scripts/psph_5mc_context_window_counts.py
- Input: psph_5mC_cov30_pm10.fa
- Output: prints six Counter dicts
# Step 10 - Pseudomonas 5mC motif enrichment 
- Run: Scripts/psph_5mc_motif_counts.py
- Required: psph_5mC_cov30_pm10.fa
- Output: Prints total windows and the three counts.
# STEP 11 — Motif context & window validation
- Count C@−1 and GAG@+1..+3, build PFM, annotate each window.
- Run: Scripts/psph_5mC_analyze_windows.py
- Required: psph_5mC_cov30_pm10.fa
- Output: psph_5mC_cov30_PFM.csv, psph_5mC_windows_annotated.csv, psph_5mC_both_Cat-1_GAG.fa, psph_5mC_not_both.fa
# Step 12 - writes a MEME-format motif, and saves a base-count plot.
 - Run: Scripts/psph_5mC_meme_and_plot.py
 - Required :psph_5mC_cov30_PFM.csv  
 - Outputs:
    - psph_5mC.meme
    - psph_5mC_pf_counts.png
    - prints PFM length, naïve consensus
# Step 13 - Visualizing the Pseudomonas 5mC motif
 - Run: Scripts/psph_5mC_plot_stackedbars.py
 - Required:psph_5mC_cov30_PFM.csv
 - Output:psph_5mC_stackedbars.png
# Step 14 - psph_5mC_site_annotation
 - Run:Scripts/psph_5mC_site_annotation.py
 - Required:
    - psph_5mC_cov30.bed
    - psph_annotation.gff
    - psph_reference.fasta   
 - Output:
    - psph_5mC_cov30_annot.csv
    - summaries 
    - psph_5mC_top_genes_tagged.csv
# Step15 - Generating  P. syringae 5mC figures and gene summary
 - Run:Scripts/psph_5mC_gene_summary_figures.py
 - Inputs 
    - psph_5mC_cov30_annot.csv           
    - psph_5mC_windows_annotated.csv
 - Outputs: 
- Fig1_5mC_by_region_small.
- Fig2_5mC_percent_by_region_small
- Fig3_5mC_functional_cats_small
- Table_S1_psph_5mC_gene_summary.csv
# Steps 16 - cross species figure
- Run: Scripts/cross_species_figs_with_ps6ma.py
- Required:
   - psph_hotspots_6mA.bed
   - hotspots_hot.csv
   - psph_5mC_windows_annotated.csv
   - psph_5mC_windows_annot.csv
   - psph_5mC_cov30_annot.csv             
- Outputs
   -  psph_6mA_hotspots.csv
   -  Cross_Fig1_counts_relaxed.png
   -  Cross_Fig2_percent_mod_box.png
   -  Cross_Fig3_region_stack_small.png
# Step 17 - cross-species methylation
- Run: Scripts:cross_species_figs_with_ps6ma.py
- Required:
  - bc_hotspots_6mA.bed
  - bc_hotspots_5mC.bed
  - psph_hotspots_6mA.bed
  - psph_5mC_windows_annotated.csv
  - bc_modifications_with_region.csv
  - psph_5mC_cov30_annot.csv
- Output:
  - hotspots_hot.csv
  - psph_6mA_hotspots.csv
  - Cross_Fig1_counts_relaxed
  - Cross_Fig2_percent_mod_box
  - Cross_Fig3_region_stack_small
  - hotspot counts and median percent_mod
# Step 18 all-mods figures with Pseudomonas 6mA region annotation
Run: Scripts/cross_species_figures_with_ps6ma_annotation.py
Required:
 - hotspots_hot.csv
 - psph_6mA_hotspots.csv
 - psph_5mC_hotspots_full.csv 
- GFF: psph_annotation.gff
Output: 
 - Cross_Fig1_counts_relaxed
 -  Cross_Fig2_percent_mod_box
 -  Cross_Fig3_region_stack_small
 -  psph_6mA_hotspots_annot.csv
 -  Cross_species_mod_summary.csv
# Step 19 Cross-species region stack
Run: Scripts/ cross_region_stack_with_PS6mA.py
- Inputs
- hotspots_hot.csv
- psph_6mA_hotspots_annot_fixed.csv  OR  psph_6mA_hotspots_annot.csv
- One of:
    psph_5mC_cov30_annot.csv
    psph_5mC_hotspots_full.csv
    psph_5mC_windows_annotated.csv
    psph_5mC_windows_annot.csv
- Outputs
- region_plot_input.csv
- Cross_Fig3_region_stack_with_PS6mA_v2.png
- Cross_Fig3_region_stack_with_PS6mA_v2.pdf
# Step: 20 cross_fig3_region_stack_CDSonly
 - Run: Scripts/cross_region_CDSonly.py
 - Required:region_plot_input.csv
 - Outputs: 
  - region_fractions_CDSonly.csv
  - Cross_Fig3_region_stack_CDSonly.png
  - Cross_Fig3_region_stack_CDSonly.pdf

