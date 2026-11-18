import numpy as np, pandas as pd, celloracle as co, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# ---- supplementary all-84 csv with the same columns as the curated csv ----
FIG="phase3_figs"
df=pd.read_csv(f"{FIG}/lipid_TF_KO_alignment_all.csv")
cur=pd.read_csv(f"{FIG}/lipid_TF_KO_main_curated.csv")[["TF","category"]]
df=df.merge(cur,on="TF",how="left")
df["in_main_figure"]=df.category.notna()
df["category"]=df.category.fillna("supplementary only")
cols=["TF","in_main_figure","category","source_group","class","cell_cos_mean","cell_cos_median","cell_frac_opposite","grid_cos_mean","grid_frac_opposite","ko_mag_mean","proj_on_natural_mean","random_cos_mean","n_cells","grid_n","outdeg_Q1","outdeg_Q2","outdeg_Q3","outdeg_Q4","outdeg_Q5","batch"]
df=df.sort_values("cell_cos_mean")
df[cols].to_csv(f"{FIG}/lipid_TF_KO_supp_all84.csv",index=False)
print("saved supp csv:", df.shape, " main-figure TFs:", int(df.in_main_figure.sum()))
