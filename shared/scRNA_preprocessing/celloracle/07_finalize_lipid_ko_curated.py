"""Curated main-figure table (narrow lipid-metabolism TFs) + supplementary (all 84), two-panel layout."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
FIG="phase3_figs"
df=pd.read_csv(f"{FIG}/lipid_TF_KO_alignment_all.csv")

MAIN={  # narrow-sense lipid-metabolism TFs (edit here)
 "Sterol / fatty-acid synthesis": ["Srebf2","Srebf1","Usf1","Usf2","Nfya","Nfyc","Mlx"],
 "Lipid-sensing nuclear receptor": ["Ppard","Rxra","Nr1d1","Nr2f2","Thra","Esrra"],
 "Lipogenic C/EBP-KLF":            ["Cebpa","Cebpb","Klf15"],
 "ER-stress / lipid homeostasis":  ["Xbp1","Atf6","Nfe2l1"],
 "Lipophagy / FOXO":               ["Tfe3","Foxo1"],
}
cat={t:c for c,ts in MAIN.items() for t in ts}
main=df[df.TF.isin(cat)].copy(); main["category"]=main.TF.map(cat)
missing=[t for t in cat if t not in set(df.TF)]
if missing: print("WARNING not in results:", missing)
main=main.sort_values("cell_cos_mean")
cols=["TF","category","class","cell_cos_mean","cell_frac_opposite","grid_cos_mean","ko_mag_mean","random_cos_mean","outdeg_Q1","outdeg_Q3","outdeg_Q5"]
main[cols].to_csv(f"{FIG}/lipid_TF_KO_main_curated.csv",index=False)
print(main[cols].round(3).to_string(index=False))

colmap={"reverses rejuvenation":"#b2182b","reverses (weak KO)":"#ef8a62",
        "accelerates rejuvenation":"#2166ac","accelerates (weak KO)":"#67a9cf","neutral":"#bdbdbd"}
rmin,rmax=df.random_cos_mean.min(),df.random_cos_mean.max()

def two_panel(d, title, out, row_h=0.32, label_cat=False):
    d=d.sort_values("cell_cos_mean").reset_index(drop=True); n=len(d)
    fig,(a,b)=plt.subplots(1,2,figsize=(10.5,row_h*n+1.6),sharey=True,gridspec_kw={"width_ratios":[3,1.4]})
    cols_=[colmap[x] for x in d["class"]]
    a.barh(range(n),d.cell_cos_mean,color=cols_,edgecolor="none")
    a.axvline(0,color="k",lw=0.8); a.axvspan(rmin,rmax,color="grey",alpha=0.3)
    a.axvline(-0.2,color="grey",ls="--",lw=0.7); a.axvline(0.2,color="grey",ls="--",lw=0.7)
    a.set_yticks(range(n)); a.set_yticklabels(d.TF,fontsize=9)
    for lab in a.get_yticklabels():
        if lab.get_text()=="Srebf2": lab.set_fontweight("bold"); lab.set_color("#b2182b")
    if label_cat:
        for i,r in d.iterrows():
            a.text(a.get_xlim()[0] if False else (min(d.cell_cos_mean.min(),-0.25)-0.02), i, r["category"], fontsize=6.5, color="dimgrey", ha="right", va="center")
    a.set_xlabel("Alignment with rejuvenation flow\nmean cosine(KO shift, natural flow), Q1–Q3 cells")
    a.set_title("A  Direction of KO effect",loc="left",fontsize=11,fontweight="bold")
    b.barh(range(n),d.ko_mag_mean,color=cols_,edgecolor="none")
    b.set_xlabel("KO effect size\nmean |shift| (Q1–Q3 cells)")
    b.set_title("B  Magnitude",loc="left",fontsize=11,fontweight="bold")
    a.set_ylim(-0.8,n-0.2)
    handles=[Patch(color=v,label=k) for k,v in colmap.items()]+[Patch(color="grey",alpha=0.3,label=f"random control ({rmin:+.2f}..{rmax:+.2f})")]
    a.legend(handles=handles,loc="lower right",fontsize=7,frameon=False)
    fig.suptitle(title,fontsize=12,fontweight="bold",y=0.995)
    fig.tight_layout()
    fig.savefig(f"{FIG}/{out}.pdf",bbox_inches="tight"); fig.savefig(f"{FIG}/{out}.png",dpi=200,bbox_inches="tight"); plt.close(fig)
    print("saved",out)

two_panel(main, f"In-silico KO of {len(main)} lipid-metabolism TFs (CellOracle, R0–R9 GRN)", "lipid_TF_KO_main_curated", row_h=0.36)
two_panel(df, f"Supplementary: in-silico KO of all {len(df)} lipid-associated regulators", "lipid_TF_KO_supp_all84", row_h=0.24)
