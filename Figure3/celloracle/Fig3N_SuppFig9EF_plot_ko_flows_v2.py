"""
plot_ko_flows_v2.py — KO perturbation flow plots on the diffusion-map trajectory embedding (from saved v2 vectors)
usage: python plot_ko_flows_v2.py Srebf1 Srebf2 Foxo1 ...   [--tag lipid7]
outputs -> phase3_figs/ps_flows_<tag>/
  <TF>_flow.png/pdf            natural (reference) flow | <TF> KO flow, full diffusion map, phase colours
  <TF>_flow_Q1Q3.png/pdf       same, zoomed to Q1–Q3 region
  panel_<tag>_full.png/pdf     natural + all KOs (full view)
  panel_<tag>_Q1Q3.png/pdf     natural + all KOs (Q1–Q3 zoom)
All panels share one arrow scale (set from the reference flow) so arrow lengths are comparable.
"""
import os, sys, glob, argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ap=argparse.ArgumentParser(); ap.add_argument("tfs",nargs="+"); ap.add_argument("--tag",default="sel")
ap.add_argument("--arrow",type=float,default=1.6,help="90th-percentile reference arrow length in grid spacings")
ap.add_argument("--scale-ref",default="Srebf2",help="KO arrow scale fixed to this TF ('pool' = pooled plotted KOs)")
args=ap.parse_args()
BASE="EGM_multiflex/oracle"; FIG=f"{BASE}/phase3_figs"
OUT=f"{FIG}/ps_flows_{args.tag}"; os.makedirs(OUT,exist_ok=True)

ref=np.load(f"{FIG}/ps_v2a/reference_flow.npz",allow_pickle=True)
grid=ref["grid"]; rf=ref["ref_flow"]; mref=~ref["mass_filter"]; emb=ref["emb"]; phases=ref["phases"].astype(str)
scores=pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{FIG}/ps_*/ps_scores.csv")) if "/ps_test/" not in f]).drop_duplicates("TF",keep="last").set_index("TF")
phase_colors={"Q1":"#440154","Q2":"#3b528b","Q3":"#21918c","Q4":"#5ec962","Q5":"#fde725"}
bg=np.array([phase_colors.get(p,"#888") for p in phases]); in_early=np.isin(phases,["Q1","Q2","Q3"])
legend=[Patch(color=c,label=p) for p,c in phase_colors.items()]

# common arrow scale from the reference flow
gsp=np.diff(np.unique(grid[:,0]))[0]
p90=np.percentile(np.linalg.norm(rf[mref],axis=1),90)
SCALE_REF=p90/(args.arrow*gsp)      # quiver scale for the reference (pseudotime) flow
SCALE_KO=None                        # set below from the pooled KO flows (shared by all KO panels)

def load_vec(tf):
    for p in sorted(glob.glob(f"{FIG}/ps_*/vectors/{tf}.npz")):
        if "/ps_test/" not in p: return np.load(p)
    return None

def draw(ax,U,V,mask,color,title,zoom=False,bgcol=None,scale=None):
    c=bg if bgcol is None else bgcol
    ax.scatter(emb[:,0],emb[:,1],c=c,s=2,alpha=0.35,edgecolors="none")
    ax.quiver(grid[mask,0],grid[mask,1],U[mask],V[mask],color=color,angles="xy",scale_units="xy",scale=scale,
              width=0.005,headwidth=4,headlength=5)
    ax.set_aspect("equal"); ax.set_xlabel("DC1"); ax.set_ylabel("DC2"); ax.set_title(title,fontsize=12,fontweight="bold")
    if zoom:
        q=emb[in_early]; pad=0.0015
        ax.set_xlim(q[:,0].min()-pad,q[:,0].max()+pad); ax.set_ylim(q[:,1].min()-pad,q[:,1].max()+pad)

vecs={tf:load_vec(tf) for tf in args.tfs}
if args.scale_ref!="pool" and load_vec(args.scale_ref) is not None:
    rv=load_vec(args.scale_ref); pool=np.linalg.norm(rv["flow"][~rv["mass_filter_sim"]],axis=1)   # fixed, subset-independent
else:
    pool=np.concatenate([np.linalg.norm(v["flow"][~v["mass_filter_sim"]],axis=1) for v in vecs.values() if v is not None])
SCALE_KO=np.percentile(pool,90)/(args.arrow*gsp)
print(f"scales: ref={SCALE_REF:.3g}  ko={SCALE_KO:.3g}  (ref/ko = {SCALE_REF/SCALE_KO:.1f}x)")
panels_full=[("Natural rejuvenation flow\n(Palantir pseudotime gradient)",rf[:,0],rf[:,1],mref,"black",SCALE_REF)]
panels_zoom=list(panels_full)
for tf in args.tfs:
    v=vecs[tf]
    if v is None: print("missing vectors:",tf); continue
    fl=v["flow"]; m=~v["mass_filter_sim"]; ps=scores.loc[tf,"mean_ps_early"]; q=scores.loc[tf,"q_neg_early"] if "q_neg_early" in scores.columns else np.nan
    lab=f"{tf} KO simulation flow\n(mean PS = {ps:+.3f})"
    for zoom,suffix in [(False,"flow"),(True,"flow_Q1Q3")]:
        fig,axes=plt.subplots(1,2,figsize=(14,6.5))
        draw(axes[0],rf[:,0],rf[:,1],mref,"black",panels_full[0][0],zoom,scale=SCALE_REF)
        draw(axes[1],fl[:,0],fl[:,1],m,"darkred",lab,zoom,scale=SCALE_KO)
        fig.legend(handles=legend,loc="lower center",ncol=5,bbox_to_anchor=(0.5,-0.02),frameon=False)
        fig.suptitle(f"{tf} in-silico KO vs natural rejuvenation — diffusion map"+("  (Q1–Q3 region)" if zoom else ""),fontsize=14,fontweight="bold",y=0.99)
        fig.tight_layout(); fig.savefig(f"{OUT}/{tf}_{suffix}.png",dpi=200,bbox_inches="tight"); fig.savefig(f"{OUT}/{tf}_{suffix}.pdf",bbox_inches="tight"); plt.close(fig)
    panels_full.append((f"{tf} KO  (PS {ps:+.3f})",fl[:,0],fl[:,1],m,"darkred",SCALE_KO))
    panels_zoom.append((f"{tf} KO  (PS {ps:+.3f})",fl[:,0],fl[:,1],m,"darkred",SCALE_KO))
    print(f"{tf}: PS={ps:+.4f}  grid pts={m.sum()}")

for panels,zoom,name in [(panels_full,False,"full"),(panels_zoom,True,"Q1Q3")]:
    ncol=4; nrow=int(np.ceil(len(panels)/ncol))
    fig,axes=plt.subplots(nrow,ncol,figsize=(5.2*ncol,5.4*nrow)); axes=np.array(axes).flatten()
    for ax,(t,U,V,m,c,sc) in zip(axes,panels): draw(ax,U,V,m,c,t,zoom,scale=sc)
    for ax in axes[len(panels):]: ax.axis("off")
    fig.legend(handles=legend,loc="lower center",ncol=5,bbox_to_anchor=(0.5,-0.01),frameon=False)
    fig.suptitle("In-silico KO perturbation flows on the rejuvenation trajectory"+("  (Q1–Q3 region)" if zoom else "")+"\n(reference and KO panels use separate arrow scales; all KO panels share one scale)",fontsize=14,fontweight="bold",y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/panel_{args.tag}_{name}.png",dpi=170,bbox_inches="tight"); fig.savefig(f"{OUT}/panel_{args.tag}_{name}.pdf",bbox_inches="tight"); plt.close(fig)
print("saved ->",OUT)

# ============================================================
# clean style: full DC map, no axes/box, grey cells, arrows only on Q1–Q3-dominant grid points
# ============================================================
from scipy.spatial import cKDTree
tree=cKDTree(emb); r=gsp*1.0
early_frac=np.zeros(len(grid))
for i,(gx,gy) in enumerate(grid):
    idx=tree.query_ball_point([gx,gy],r)
    early_frac[i]=np.mean(in_early[idx]) if len(idx) else 0
q13=early_frac>=0.5

def draw_clean(ax,U,V,mask,color,title,scale,arrow_w=0.0045):
    ax.scatter(emb[:,0],emb[:,1],c="#d9d9d9",s=3,alpha=0.6,edgecolors="none")
    ax.scatter(emb[in_early,0],emb[in_early,1],c="#a6a6a6",s=3,alpha=0.7,edgecolors="none")
    mm=mask&q13
    ax.quiver(grid[mm,0],grid[mm,1],U[mm],V[mm],color=color,angles="xy",scale_units="xy",scale=scale,
              width=arrow_w,headwidth=4,headlength=5,minlength=0.5)
    ax.set_aspect("equal"); ax.axis("off")
    if title: ax.set_title(title,fontsize=13,fontweight="bold",color=color if color!="black" else "black")

# per-TF single panels
for tf in args.tfs:
    v=vecs[tf]
    if v is None: continue
    fl=v["flow"]; m=~v["mass_filter_sim"]; ps=scores.loc[tf,"mean_ps_early"]
    fig,ax=plt.subplots(figsize=(6.5,6.5))
    draw_clean(ax,fl[:,0],fl[:,1],m,"#b2182b",f"{tf} KO",SCALE_KO)
    fig.tight_layout(); fig.savefig(f"{OUT}/{tf}_clean.png",dpi=250,bbox_inches="tight",facecolor="white"); fig.savefig(f"{OUT}/{tf}_clean.pdf",bbox_inches="tight"); plt.close(fig)
    # natural vs KO pair
    fig,axes=plt.subplots(1,2,figsize=(12.5,6.5))
    draw_clean(axes[0],rf[:,0],rf[:,1],mref,"black","Natural rejuvenation flow",SCALE_REF)
    draw_clean(axes[1],fl[:,0],fl[:,1],m,"#b2182b",f"{tf} KO",SCALE_KO)
    fig.tight_layout(); fig.savefig(f"{OUT}/{tf}_clean_pair.png",dpi=250,bbox_inches="tight",facecolor="white"); fig.savefig(f"{OUT}/{tf}_clean_pair.pdf",bbox_inches="tight"); plt.close(fig)

# combined clean panel: natural + KOs
items=[("Natural rejuvenation flow",rf[:,0],rf[:,1],mref,"black",SCALE_REF)]
for tf in args.tfs:
    v=vecs[tf]
    if v is None: continue
    items.append((f"{tf} KO",v["flow"][:,0],v["flow"][:,1],~v["mass_filter_sim"],"#b2182b",SCALE_KO))
ncol=4; nrow=int(np.ceil(len(items)/ncol))
fig,axes=plt.subplots(nrow,ncol,figsize=(5*ncol,5.2*nrow)); axes=np.array(axes).flatten()
for ax,(t,U,V,m,c,sc) in zip(axes,items): draw_clean(ax,U,V,m,c,t,sc,arrow_w=0.005)
for ax in axes[len(items):]: ax.axis("off")
fig.tight_layout(); fig.savefig(f"{OUT}/panel_{args.tag}_clean.png",dpi=200,bbox_inches="tight",facecolor="white"); fig.savefig(f"{OUT}/panel_{args.tag}_clean.pdf",bbox_inches="tight"); plt.close(fig)
print("clean style saved:", f"{OUT}/panel_{args.tag}_clean.png", f"(Q1-Q3 grid points: {int(q13.sum())})")
