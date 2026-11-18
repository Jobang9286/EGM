"""
celloracle_ko_batch_v2.py — in-silico TF KO scored with CellOracle's official perturbation score (PS)

Reference flow : Gradient_calculator on palantir_pseudotime (X_diffmap2d), knn transfer, n_grid=40
KO flow        : simulate_shift(n_propagation=3) -> estimate_transition_prob(n_neighbors=200, knn_random) -> calculate_embedding_shift(sigma_corr=0.05)
Score          : Oracle_development_module inner product (PS) per grid point; negative/positive PS sums + Wilcoxon p-value
                 vs randomized simulation (get_negative_PS_p_value / get_positive_PS_p_value), computed for
                 (a) early region  = grid pseudotime bins covering Q1–Q3  (palantir pt <= --early-pt-max)
                 (b) whole trajectory
Also saved     : per-cell delta_embedding (+random), grid flows, inner products  -> vectors/<TF>.npz
                 4-panel figure per TF (reference flow / KO flow / PS on grid / random PS) -> figs/<TF>_PS.png
usage: python celloracle_ko_batch_v2.py --tag A Srebf2 Klf4 ...   |  --tf-file list.txt  |  --no-ref
"""
import os, sys, time, argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import celloracle as co
from celloracle.applications import Gradient_calculator, Oracle_development_module

t0=time.time()
def stamp(m): print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

ap=argparse.ArgumentParser()
ap.add_argument("tfs", nargs="*")
ap.add_argument("--tf-file")
ap.add_argument("--tag", default=time.strftime("%Y%m%d_%H%M"))
ap.add_argument("--no-ref", action="store_true")
ap.add_argument("--n-jobs", type=int, default=16)
ap.add_argument("--n-grid", type=int, default=40)
ap.add_argument("--min-mass", type=float, default=2000.0, help="total_p_mass threshold; on the diffmap scale ~2000 hugs the cell cloud")
ap.add_argument("--early-pt-max", type=float, default=0.572, help="palantir pt upper bound of Q3")
args=ap.parse_args()

TFS=list(args.tfs)
if args.tf_file: TFS+=[l.strip() for l in open(args.tf_file) if l.strip() and not l.startswith("#")]
if not args.no_ref and "Srebf2" not in TFS: TFS=["Srebf2"]+TFS
if not TFS: sys.exit("no TFs")

BASE="."
OUT=f"{BASE}/phase3_figs/ps_{args.tag}"
for d in ["vectors","figs"]: os.makedirs(f"{OUT}/{d}", exist_ok=True)

# ---- 1. oracle + GRN ----
stamp(f"celloracle {co.__version__}  tag={args.tag}  TFs={TFS}")
oracle=co.load_hdf5(f"{BASE}/r1r9.celloracle.oracle")
links=co.load_hdf5(f"{BASE}/r1r9.celloracle.links")
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
stamp("GRN fit done")
emb=oracle.adata.obsm["X_diffmap2d"]; pt=oracle.adata.obs["palantir_pseudotime"].values
phases=oracle.adata.obs["phase_q"].astype(str).values
in_early=np.isin(phases,["Q1","Q2","Q3"])
reg_union=set()
for cl,df_ in links.filtered_links.items(): reg_union|=set(df_["source"].unique())
usable=[]
for tf in TFS:
    if tf not in oracle.adata.var_names: stamp(f"SKIP {tf}: not in adata")
    elif tf not in reg_union: stamp(f"SKIP {tf}: not a regulator in filtered GRN")
    else: usable.append(tf)
if not usable: sys.exit("no usable TFs")

# ---- 2. reference (developmental) flow ----
g=Gradient_calculator(oracle_object=oracle, pseudotime_key="palantir_pseudotime")
g.calculate_p_mass(smooth=0.8, n_grid=args.n_grid, n_neighbors=200)
g.calculate_mass_filter(min_mass=args.min_mass, plot=False)
g.transfer_data_into_grid(args={"method":"knn","n_knn":50}, plot=False)
g.calculate_gradient()
np.savez(f"{OUT}/reference_flow.npz", grid=g.gridpoints_coordinates, ref_flow=g.ref_flow,
         mass_filter=g.mass_filter, pt_grid=g.pseudotime_on_grid, emb=emb, pt=pt, phases=phases)
stamp(f"reference flow done: {int((~g.mass_filter).sum())} grid points pass mass filter")
try:
    fig,axes=plt.subplots(1,2,figsize=(14,6))
    sc_=axes[0].scatter(emb[:,0],emb[:,1],c=pt,cmap="viridis",s=2,edgecolors="none"); axes[0].set_aspect("equal"); plt.colorbar(sc_,ax=axes[0],fraction=0.04); axes[0].set_title("Palantir pseudotime")
    g.plot_reference_flow_on_grid(ax=axes[1], scale=30, s=3); axes[1].set_title("Reference (developmental) flow")
    fig.tight_layout(); fig.savefig(f"{OUT}/reference_flow.png", dpi=150, bbox_inches="tight"); plt.close(fig)
except Exception as e:
    stamp(f"reference plot failed: {e}")

# ---- 3. natural KNN flow (legacy cosine metric, for continuity with batch 1) ----
K=50; tree=cKDTree(emb); _,idx=tree.query(emb,k=K+1); idx=idx[:,1:]
nat=np.zeros_like(emb)
for i in range(len(emb)):
    w=np.clip(pt[idx[i]]-pt[i],0,None)
    if w.sum()<1e-8: continue
    nat[i]=((emb[idx[i]]-emb[i])*w[:,None]).sum(0)/w.sum()
nn_=np.linalg.norm(nat,axis=1); nat_u=np.divide(nat,nn_[:,None],out=np.zeros_like(nat),where=nn_[:,None]>0)

def cell_cos(vec):
    vn=np.linalg.norm(vec,axis=1); ok=in_early&(vn>1e-12)&(nn_>1e-12)
    c=(vec[ok]*nat_u[ok]).sum(1)/vn[ok]
    return float(c.mean()), float((c<0).mean()), float(vn[in_early].mean())

# ---- 4. per-TF ----
rows=[]
phase_colors={"Q1":"#440154","Q2":"#3b528b","Q3":"#21918c","Q4":"#5ec962","Q5":"#fde725"}
bg=np.array([phase_colors.get(p,"#888") for p in phases])
for tf in usable:
    stamp(f"=== {tf} ===")
    oracle.simulate_shift(perturb_condition={tf:0.0}, n_propagation=3, ignore_warning=True)
    oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1, n_jobs=args.n_jobs)
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    de=np.asarray(oracle.delta_embedding,dtype=np.float32); der=np.asarray(oracle.delta_embedding_random,dtype=np.float32)

    dev=Oracle_development_module()
    dev.load_differentiation_reference_data(gradient_object=g)
    dev.load_perturb_simulation_data(oracle_object=oracle)
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)
    df=dev.inner_product_df
    early_ids=sorted(df.loc[df.pseudotime<=args.early_pt_max,"pseudotime_id"].unique().tolist())
    early_str=",".join(str(int(i)) for i in early_ids)
    p_neg_e,neg_e,neg_re=dev.get_negative_PS_p_value(pseudotime=early_str, return_ps_sum=True)
    p_pos_e,pos_e,pos_re=dev.get_positive_PS_p_value(pseudotime=early_str, return_ps_sum=True)
    p_neg_a,neg_a,neg_ra=dev.get_negative_PS_p_value(return_ps_sum=True)
    p_pos_a,pos_a,pos_ra=dev.get_positive_PS_p_value(return_ps_sum=True)
    de_=df[df.pseudotime_id.isin(early_ids)]
    cc,fo,mag=cell_cos(de)
    row=dict(TF=tf, early_bins=early_str, n_grid_early=int(len(de_)), n_grid_all=int(len(df)),
             mean_ps_early=float(de_.score.mean()), mean_ps_rand_early=float(de_.score_randomized.mean()),
             neg_ps_early=float(neg_e), neg_ps_rand_early=float(neg_re), p_neg_early=float(p_neg_e),
             pos_ps_early=float(pos_e), pos_ps_rand_early=float(pos_re), p_pos_early=float(p_pos_e),
             mean_ps_all=float(df.score.mean()), mean_ps_rand_all=float(df.score_randomized.mean()),
             neg_ps_all=float(neg_a), neg_ps_rand_all=float(neg_ra), p_neg_all=float(p_neg_a),
             pos_ps_all=float(pos_a), pos_ps_rand_all=float(pos_ra), p_pos_all=float(p_pos_a),
             cell_cos_mean=cc, cell_frac_opposite=fo, ko_mag_mean=mag)
    rows.append(row)
    stamp(f"  {tf}: mean PS(early)={row['mean_ps_early']:+.2e} (rand {row['mean_ps_rand_early']:+.2e})  "
          f"negPS={neg_e:.2e} vs rand {neg_re:.2e} p={p_neg_e:.1e} | posPS={pos_e:.2e} vs {pos_re:.2e} p={p_pos_e:.1e} | cos={cc:+.3f}")

    np.savez_compressed(f"{OUT}/vectors/{tf}.npz", delta_embedding=de, delta_embedding_random=der,
        flow=dev.flow, flow_rndm=dev.flow_rndm, mass_filter_sim=dev.mass_filter_simulation,
        inner_product=dev.inner_product, inner_product_random=dev.inner_product_random,
        pt_grid=dev.pseudotime_on_grid, early_ids=np.array(early_ids))

    # 4-panel figure
    try:
        grid=g.gridpoints_coordinates; m=~dev.mass_filter_simulation
        ip=dev.inner_product; vm=np.nanpercentile(np.abs(ip[m]),98) if m.any() else 1
        fig,axes=plt.subplots(1,4,figsize=(24,6))
        for ax in axes:
            ax.scatter(emb[:,0],emb[:,1],c="#dddddd",s=1,alpha=0.5,edgecolors="none"); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        mr=~g.mass_filter
        axes[0].quiver(grid[mr,0],grid[mr,1],g.ref_flow[mr,0],g.ref_flow[mr,1],color="black",angles="xy",scale_units="xy",width=0.004)
        axes[0].set_title("Reference flow (Palantir pseudotime)")
        axes[1].quiver(grid[m,0],grid[m,1],dev.flow[m,0],dev.flow[m,1],color="darkred",angles="xy",scale_units="xy",width=0.004)
        axes[1].set_title(f"{tf} KO simulation flow")
        s1=axes[2].scatter(grid[m,0],grid[m,1],c=ip[m],cmap="RdBu_r",vmin=-vm,vmax=vm,s=45,edgecolors="none"); axes[2].set_title(f"{tf} KO: perturbation score (inner product)")
        plt.colorbar(s1,ax=axes[2],fraction=0.04)
        s2=axes[3].scatter(grid[m,0],grid[m,1],c=dev.inner_product_random[m],cmap="RdBu_r",vmin=-vm,vmax=vm,s=45,edgecolors="none"); axes[3].set_title("Randomized simulation: PS")
        plt.colorbar(s2,ax=axes[3],fraction=0.04)
        fig.suptitle(f"{tf} KO — early negPS={neg_e:.2e} (rand {neg_re:.2e}), p={p_neg_e:.1e}",fontsize=13,fontweight="bold")
        fig.tight_layout(); fig.savefig(f"{OUT}/figs/{tf}_PS.png",dpi=130,bbox_inches="tight"); plt.close(fig)
    except Exception as e:
        stamp(f"  figure failed for {tf}: {e}")
    pd.DataFrame(rows).to_csv(f"{OUT}/ps_scores.csv",index=False)

pd.DataFrame(rows).to_csv(f"{OUT}/ps_scores.csv",index=False)
stamp(f"DONE — {len(rows)} TFs -> {OUT}")
