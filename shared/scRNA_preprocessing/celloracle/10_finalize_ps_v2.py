"""Final figures/tables from v2 (official CellOracle PS) batches: curated main + supplementary + top-hit PS panel."""
import os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from statsmodels.stats.multitest import multipletests

BASE="."; FIG=f"{BASE}/phase3_figs"
dfs=[]
for tag in ["v2a","v2b"]:
    p=f"{FIG}/ps_{tag}/ps_scores.csv"
    if os.path.exists(p):
        d=pd.read_csv(p); d["batch"]=tag; dfs.append(d)
df=pd.concat(dfs).drop_duplicates("TF").reset_index(drop=True)
print("TFs scored:", len(df))

# --- multiple testing (BH across TFs) ---
for k in ["p_neg_early","p_pos_early","p_neg_all","p_pos_all"]:
    df[k]=df[k].fillna(1.0)
    df["q"+k[1:]]=multipletests(df[k],method="fdr_bh")[1]
df["negPS_fold_early"]=df.neg_ps_early/df.neg_ps_rand_early.replace(0,np.nan)
df["posPS_fold_early"]=df.pos_ps_early/df.pos_ps_rand_early.replace(0,np.nan)

PS_THR=0.05   # |mean PS| effect-size threshold (~50% of the Srebf2 effect, natural break in the distribution)
def cls(r):
    if r.q_neg_early<0.05 and r.mean_ps_early<=-PS_THR: return "reverses rejuvenation"
    if r.q_pos_early<0.05 and r.mean_ps_early>= PS_THR: return "accelerates rejuvenation"
    return "not significant"
df["class"]=df.apply(cls,axis=1)
def stars(q): return "***" if q<0.001 else "**" if q<0.01 else "*" if q<0.05 else ""
df["sig_neg"]=df.q_neg_early.map(stars); df["sig_pos"]=df.q_pos_early.map(stars)
df["cosine_agrees"]=((np.sign(df.mean_ps_early)==np.sign(df.cell_cos_mean))&(df.cell_cos_mean.abs()>=0.2)).map({True:"yes",False:"no"})
df["mean_abs_ps_early"]=(df.neg_ps_early+df.pos_ps_early)/df.n_grid_early

# --- annotation (category / out-degree) from batch-1 tables ---
old=pd.read_csv(f"{FIG}/lipid_TF_KO_supp_all84.csv")[["TF","in_main_figure","category","source_group","outdeg_Q1","outdeg_Q2","outdeg_Q3","outdeg_Q4","outdeg_Q5"]]
df=df.merge(old,on="TF",how="left")
df["in_main_figure"]=df.in_main_figure.fillna(False).astype(bool)
df=df.sort_values("mean_ps_early")
cols=["TF","in_main_figure","category","source_group","class",
      "mean_ps_early","mean_ps_rand_early","neg_ps_early","neg_ps_rand_early","negPS_fold_early","p_neg_early","q_neg_early","sig_neg",
      "pos_ps_early","pos_ps_rand_early","posPS_fold_early","p_pos_early","q_pos_early","sig_pos","n_grid_early","early_bins",
      "mean_ps_all","neg_ps_all","neg_ps_rand_all","q_neg_all","pos_ps_all","pos_ps_rand_all","q_pos_all","n_grid_all",
      "mean_abs_ps_early","cosine_agrees","cell_cos_mean","cell_frac_opposite","ko_mag_mean","outdeg_Q1","outdeg_Q2","outdeg_Q3","outdeg_Q4","outdeg_Q5","batch"]
df[cols].to_csv(f"{FIG}/PS_TF_KO_all.csv",index=False)
main=df[df.in_main_figure].copy(); main[cols].to_csv(f"{FIG}/PS_TF_KO_main_curated.csv",index=False)
print(df["class"].value_counts().to_string())

colmap={"reverses rejuvenation":"#b2182b","accelerates rejuvenation":"#2166ac","not significant":"#bdbdbd"}
def two_panel(d,title,out,row_h=0.34):
    d=d.sort_values("mean_ps_early").reset_index(drop=True); n=len(d); cols_=[colmap[x] for x in d["class"]]
    fig,(a,b)=plt.subplots(1,2,figsize=(11,row_h*n+1.8),sharey=True,gridspec_kw={"width_ratios":[3,1.6]})
    a.barh(range(n),d.mean_ps_early,color=cols_,edgecolor="none")
    a.axvline(0,color="k",lw=0.8)
    rmin,rmax=d.mean_ps_rand_early.min(),d.mean_ps_rand_early.max(); a.axvspan(rmin,rmax,color="grey",alpha=0.3)
    a.set_yticks(range(n)); a.set_yticklabels(d.TF,fontsize=9)
    for lab in a.get_yticklabels():
        if lab.get_text()=="Srebf2": lab.set_fontweight("bold"); lab.set_color("#b2182b")
    xr=max(abs(d.mean_ps_early.min()),abs(d.mean_ps_early.max()))*1.25
    for i,r in d.iterrows():
        s=r.sig_neg if r.mean_ps_early<0 else r.sig_pos
        if s and r["class"]!="not significant":
            a.text(r.mean_ps_early+(-0.01 if r.mean_ps_early<0 else 0.01)*xr, i, s, va="center", ha="right" if r.mean_ps_early<0 else "left", fontsize=8)
    a.set_xlim(-xr,xr)
    a.set_xlabel("Perturbation score (PS)\nmean inner product of KO shift with developmental flow, Q1–Q3 grid")
    a.set_title("A  Direction of KO effect",loc="left",fontsize=11,fontweight="bold")
    # B: negative / positive PS sum vs randomized
    neg=-d.neg_ps_early.values; pos=d.pos_ps_early.values
    b.barh(range(n),neg,color=[colmap["reverses rejuvenation"] if c=="reverses rejuvenation" else "#e0b4b4" for c in d["class"]],edgecolor="none")
    b.barh(range(n),pos,color=[colmap["accelerates rejuvenation"] if c=="accelerates rejuvenation" else "#b4c6e0" for c in d["class"]],edgecolor="none")
    b.scatter(-d.neg_ps_rand_early,range(n),marker="|",color="k",s=40,lw=1,zorder=3)
    b.scatter(d.pos_ps_rand_early,range(n),marker="|",color="k",s=40,lw=1,zorder=3)
    b.axvline(0,color="k",lw=0.8)
    b.set_xlabel("Σ negative PS (left) / Σ positive PS (right)\nblack tick = randomized simulation")
    b.set_title("B  Magnitude vs random",loc="left",fontsize=11,fontweight="bold")
    a.set_ylim(-0.8,n-0.2)
    handles=[Patch(color=v,label=k) for k,v in colmap.items()]+[Patch(color="grey",alpha=0.3,label="random-control PS range")]
    a.legend(handles=handles,loc="lower right",fontsize=7,frameon=False)
    a.text(0.01,0.995,"*, **, *** : BH q < 0.05, 0.01, 0.001 (paired Wilcoxon vs randomized simulation)",transform=a.transAxes,fontsize=6.5,va="top",color="dimgrey")
    fig.suptitle(title,fontsize=12,fontweight="bold",y=0.995); fig.tight_layout()
    fig.savefig(f"{FIG}/{out}.pdf",bbox_inches="tight"); fig.savefig(f"{FIG}/{out}.png",dpi=200,bbox_inches="tight"); plt.close(fig); print("saved",out)

if len(main): two_panel(main,f"In-silico KO of {len(main)} lipid-metabolism TFs — CellOracle perturbation score","PS_TF_KO_main_curated",row_h=0.38)
two_panel(df,f"Supplementary: in-silico KO of all {len(df)} lipid-associated regulators — perturbation score","PS_TF_KO_supp_all",row_h=0.25)

# --- top-hit PS-on-grid panel from saved vectors ---
ref=None
for tag in ["v2a","v2b"]:
    p=f"{FIG}/ps_{tag}/reference_flow.npz"
    if os.path.exists(p): ref=np.load(p,allow_pickle=True); break
if ref is not None:
    grid=ref["grid"]; emb=ref["emb"]; rf=ref["ref_flow"]; mr=~ref["mass_filter"]
    top_rev=df[df["class"]=="reverses rejuvenation"].sort_values("mean_ps_early").head(8)
    top_acc=df[df["class"]=="accelerates rejuvenation"].sort_values("mean_ps_early",ascending=False).head(4)
    sel=pd.concat([top_rev,top_acc]); sel=pd.concat([df[df.TF=="Srebf2"],sel]).drop_duplicates("TF")
    panels=[]
    for r in sel.itertuples():
        p=f"{FIG}/ps_{r.batch}/vectors/{r.TF}.npz"
        if os.path.exists(p): panels.append((r,np.load(p)))
    if panels:
        vm=np.nanpercentile(np.concatenate([np.abs(v["inner_product"][~v["mass_filter_sim"]]) for _,v in panels]),98)
        ncol=4; nrow=int(np.ceil((len(panels)+1)/ncol))
        fig,axes=plt.subplots(nrow,ncol,figsize=(4.6*ncol,4.6*nrow)); axes=np.array(axes).flatten()
        ax=axes[0]; ax.scatter(emb[:,0],emb[:,1],c="#dddddd",s=1,edgecolors="none")
        ax.quiver(grid[mr,0],grid[mr,1],rf[mr,0],rf[mr,1],color="black",angles="xy",scale_units="xy",width=0.004)
        ax.set_title("Reference flow (Palantir pseudotime)",fontsize=10,fontweight="bold")
        for ax,(r,v) in zip(axes[1:],panels):
            m=~v["mass_filter_sim"]
            ax.scatter(emb[:,0],emb[:,1],c="#dddddd",s=1,edgecolors="none")
            sc_=ax.scatter(grid[m,0],grid[m,1],c=v["inner_product"][m],cmap="RdBu_r",vmin=-vm,vmax=vm,s=22,edgecolors="none")
            ax.quiver(grid[m,0],grid[m,1],v["flow"][m,0],v["flow"][m,1],color="black",alpha=0.5,angles="xy",scale_units="xy",width=0.003)
            col="#b2182b" if r.mean_ps_early<0 else "#2166ac"
            ax.set_title(f"{r.TF} KO  PS={r.mean_ps_early:+.3f} {r.sig_neg if r.mean_ps_early<0 else r.sig_pos}",fontsize=10,fontweight="bold",color=col)
        for ax in axes: ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for ax in axes[len(panels)+1:]: ax.axis("off")
        fig.colorbar(sc_,ax=axes.tolist(),fraction=0.015,pad=0.01,label="perturbation score (inner product)")
        fig.suptitle("KO flow (arrows) and perturbation score on grid: strongest reversing (red) and accelerating (blue) KOs",fontsize=12,fontweight="bold")
        fig.savefig(f"{FIG}/PS_TF_KO_top_panel.pdf",bbox_inches="tight"); fig.savefig(f"{FIG}/PS_TF_KO_top_panel.png",dpi=150,bbox_inches="tight"); plt.close(fig); print("saved PS_TF_KO_top_panel")
print(df[["TF","class","mean_ps_early","negPS_fold_early","q_neg_early","posPS_fold_early","q_pos_early","cell_cos_mean"]].round(4).to_string(index=False))
