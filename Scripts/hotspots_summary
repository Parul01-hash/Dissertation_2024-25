#Building a cross-species methylation summary table.

# Inputs 
# df_6ma_std (Bacillus 6mA)
# df_5mc_std (Bacillus 5mC)
# df_pse_6ma_std (Pseudomonas 6mA)
# df_pse_5mc_std (Pseudomonas 5mC)

import pandas as pd

def build_row(df, species, mod):
    tmp = df.dropna(subset=["cov","percent_mod"]).copy()
    tmp["percent_mod"] = tmp["percent_mod"].clip(0, 100)
    tmp["mod_reads"] = tmp["cov"] * (tmp["percent_mod"] / 100.0)
    cov_sum = float(tmp["cov"].sum())
    mod_sum = float(tmp["mod_reads"].sum())
    return pd.DataFrame([{
        "Species": species,
        "mod": mod,  # 'a' 6mA, 'm' 5mC
        "Sites": int(tmp.shape[0]),
        "Sum_valid_cov": int(round(cov_sum)),
        "Sum_mod": int(round(mod_sum)),
        "Global_frac_modified": (mod_sum / cov_sum) if cov_sum else 0.0
    }])

rows = []
if 'df_6ma_std'     in locals(): rows.append(build_row(df_6ma_std,     "B_cereus",                   "a"))
if 'df_5mc_std'     in locals(): rows.append(build_row(df_5mc_std,     "B_cereus",                   "m"))
if 'df_pse_6ma_std' in locals(): rows.append(build_row(df_pse_6ma_std, "P_syringae_pv_phaseolicola", "a"))
if 'df_pse_5mc_std' in locals() and df_pse_5mc_std[["cov","percent_mod"]].notna().all(axis=1).any():
    rows.append(build_row(df_pse_5mc_std, "P_syringae_pv_phaseolicola", "m"))

summary_hotspots = pd.concat(rows, ignore_index=True)
summary_hotspots.to_csv("cross_species_global_summary_hotspots.tsv", sep="\t", index=False)
summary_hotspots93355	71331811	0.367135
1	B_cereus	m	23906	5776817	2323846	0.402271
2	P_syringae_pv_phaseolicola	a	4094521	349386838	349380207	0.999981
