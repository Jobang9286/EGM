"""Merge lipid/lipid2a/lipid2b KO batches -> one annotated table + ranking figure + top-hit panel."""
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.image as mpimg

BASE="."; FIG=f"{BASE}/phase3_figs"
dfs=[]
for tag in ["lipid","lipid2a","lipid2b"]:
    d=pd.read_csv(f"{FIG}/batch_{tag}/ko_alignment_scores.csv"); d["batch"]=tag; dfs.append(d)
df=pd.concat(dfs).drop_duplicates("TF").reset_index(drop=True)

# --- annotation: GO category (batch2) or curated group (batch1) ---
go=pd.read_csv(f"{BASE}/lipid_go_genes.csv").set_index("SYMBOL")["category"]
manual={
 "SREBP/cholesterol":["Srebf2","Srebf1","Nfya","Nfyb","Nfyc","Sp1"],
 "nuclear receptor":["Ppard","Rxra","Nr2f2","Nr1d1","Nr1d2","Rora","Thra","Esrra","Nr3c1"],
 "ChREBP":["Mlx"],
 "adipogenic/CEBP":["Cebpa","Cebpb","Cebpd","Klf15","Klf5","Ebf1","Bhlhe40"],
 "ER stress/UPR":["Xbp1","Atf4","Atf6","Creb3l1","Creb3l2","Ddit3"],
 "lysosome/lipophagy":["Tfe3","Mitf"],
 "oxidative/FOXO":["Nfe2l1","Nfe2l2","Foxo1","Foxo3"],
}
m2g={t:g for g,ts in manual.items() for t in ts}
df["source_group"]=[m2g.get(t, "GO:"+str(go.get(t,""))) for t in df.TF]

# --- GRN out-degree per phase ---
c=pd.read_csv(f"{BASE}/phase2_centrality.csv")
piv=c[c.degree_centrality_out>0].pivot_table(index="gene",columns="cluster",values="degree_centrality_out").fillna(0).astype(int)
for q in ["Q1","Q2","Q3","Q4","Q5"]:
    df[f"outdeg_{q}"]=df.TF.map(piv[q]).fillna(0).astype(int)

# --- classification ---
def cls(r):
    strong = r.ko_mag_mean>=0.05
    if r.cell_cos_mean<=-0.2: return "reverses rejuvenation" if strong else "reverses (weak KO)"
    if r.cell_cos_mean>=0.2:  return "accelerates rejuvenation" if strong else "accelerates (weak KO)"
    return "neutral"
df["class"]=df.apply(cls,axis=1)
df=df.sort_values("cell_cos_mean")
cols=["TF","class","cell_cos_mean","cell_cos_median","cell_frac_opposite","grid_cos_mean","grid_frac_opposite","ko_mag_mean","proj_on_natural_mean","random_cos_mean","n_cells","grid_n","outdeg_Q1","outdeg_Q2","outdeg_Q3","outdeg_Q4","outdeg_Q5","source_group","batch"]
df[cols].to_csv(f"{FIG}/lipid_TF_KO_alignment_all.csv",index=False)
print(df["class"].value_counts().to_string())

# --- ranking figure ---
colmap={"reverses rejuvenation":"#b2182b","reverses (weak KO)":"#ef8a62",
        "accelerates rejuvenation":"#2166ac","accelerates (weak KO)":"#67a9cf","neutral":"#bdbdbd"}
d=df.sort_values("cell_cos_mean",ascending=True).reset_index(drop=True)
fig,ax=plt.subplots(figsize=(9,0.27*len(d)+1.5))
ax.barh(range(len(d)), d.cell_cos_mean, color=[colmap[x] for x in d["class"]], edgecolor="none")
ax.scatter(d.cell_cos_mean, range(len(d)), s=8+900*d.ko_mag_mean, facecolors="none", edgecolors="black", linewidths=0.6, zorder=3)
ax.set_yticks(range(len(d))); ax.set_yticklabels(d.TF, fontsize=8)
for lab in ax.get_yticklabels():
    if lab.get_text()=="Srebf2": lab.set_fontweight("bold")
ax.axvline(0,color="k",lw=0.8)
rmin,rmax=df.random_cos_mean.min(),df.random_cos_mean.max()
ax.axvspan(rmin,rmax,color="grey",alpha=0.3,label=f"random-control range ({rmin:+.2f}..{rmax:+.2f})")
ax.axvline(-0.2,color="grey",ls="--",lw=0.7); ax.axvline(0.2,color="grey",ls="--",lw=0.7)
ax.set_xlabel("mean cosine(KO shift, natural rejuvenation flow)  in Q1–Q3 cells\n<0: KO pushes cells back toward aged state   >0: KO accelerates rejuvenation")
ax.set_title(f"In-silico KO of {len(d)} lipid-associated regulators (CellOracle, R0–R9 GRN)", fontsize=12, fontweight="bold")
handles=[Patch(color=v,label=k) for k,v in colmap.items()]+[Patch(color="grey",alpha=0.3,label="random-control range")]
ax.legend(handles=handles,loc="lower right",fontsize=8,frameon=False)
ax.text(0.01,0.995,"circle size = mean |KO shift|",transform=ax.transAxes,ha="left",va="top",fontsize=7,color="dimgrey")
ax.set_ylim(-1,len(d)); ax.margins(x=0.05)
fig.tight_layout()
fig.savefig(f"{FIG}/lipid_TF_KO_ranking.pdf",bbox_inches="tight"); fig.savefig(f"{FIG}/lipid_TF_KO_ranking.png",dpi=200,bbox_inches="tight"); plt.close(fig)

# --- top-hit panel composed from existing per-TF figures (KO half only) ---
def ko_half(tf,batch):
    img=mpimg.imread(f"{FIG}/batch_{batch}/{tf}_Q1Q3_natural_vs_KO.png")
    h,w,_=img.shape; return img[int(0.08*h):, int(0.52*w):]
top_rev=df[df["class"]=="reverses rejuvenation"].sort_values("cell_cos_mean").head(9)
top_acc=df[df["class"]=="accelerates rejuvenation"].sort_values("cell_cos_mean",ascending=False).head(3)
sel=pd.concat([top_rev,top_acc])
nat=mpimg.imread(f"{FIG}/batch_lipid/Srebf2_Q1Q3_natural_vs_KO.png"); h,w,_=nat.shape; nat=nat[int(0.08*h):, :int(0.49*w)]
panels=[("Natural rejuvenation flow",nat,"black")]+[(f"{r.TF} KO  (cos {r.cell_cos_mean:+.2f}, {r.cell_frac_opposite*100:.0f}% opp.)",ko_half(r.TF,r.batch),"darkred" if r.cell_cos_mean<0 else "navy") for r in sel.itertuples()]
ncol=4; nrow=int(np.ceil(len(panels)/ncol))
fig,axes=plt.subplots(nrow,ncol,figsize=(4.2*ncol,4.4*nrow)); axes=np.array(axes).flatten()
for ax,(t,im,col) in zip(axes,panels):
    ax.imshow(im); ax.set_title(t,fontsize=10,fontweight="bold",color=col); ax.axis("off")
for ax in axes[len(panels):]: ax.axis("off")
fig.suptitle("Q1→Q3 region: strongest reversing (red) and accelerating (blue) lipid-regulator KOs",fontsize=13,fontweight="bold",y=1.0)
fig.tight_layout()
fig.savefig(f"{FIG}/lipid_TF_KO_top_panel.pdf",bbox_inches="tight"); fig.savefig(f"{FIG}/lipid_TF_KO_top_panel.png",dpi=150,bbox_inches="tight"); plt.close(fig)
print("saved:", f"{FIG}/lipid_TF_KO_alignment_all.csv", f"{FIG}/lipid_TF_KO_ranking.png", f"{FIG}/lipid_TF_KO_top_panel.png")
print(df[df["class"]!="neutral"][["TF","class","cell_cos_mean","cell_frac_opposite","grid_cos_mean","ko_mag_mean","outdeg_Q1","outdeg_Q5","source_group"]].round(3).to_string(index=False))
