"""
celloracle_ko_batch.py — batch in-silico KO for arbitrary TFs on the R0–R9 CellOracle object
(same GRN / embedding / natural-flow definition as celloracle_srebf2_Q1Q3.py)

usage:
  python celloracle_ko_batch.py Srebf1 Xbp1 Cebpb            # KO these TFs (+ Srebf2 reference)
  python celloracle_ko_batch.py --tag lipid --no-ref Nr2f2   # custom output tag, skip Srebf2 ref
  python celloracle_ko_batch.py --fast Atf4                  # n_neighbors=100, sampled_fraction=0.5 (~2x faster)

outputs -> phase3_figs/batch_<tag>/
  <TF>_Q1Q3_natural_vs_KO.{png,pdf}   zoomed Q1->Q3 region, natural (black) vs KO (darkred)
  <TF>_KO_fullview.png                whole diffusion map, KO quiver, phase colours
  batch_<tag>_panel.{png,pdf}         natural + every KO, zoomed
  ko_alignment_scores.csv             per-TF quantitative summary (see columns below)

score columns (all restricted to Q1–Q3 cells unless noted):
  cell_cos_mean / cell_cos_median     cosine(KO delta_embedding, natural flow) per cell, averaged
  cell_frac_opposite                  fraction of cells with cosine < 0  (KO pushes against rejuvenation)
  grid_cos_mean / grid_frac_opposite  same on the Q1–Q3 grid used for the figure
  ko_mag_mean                         mean |delta_embedding| (how strongly the KO moves cells)
  proj_on_natural_mean                mean of (KO · natural_unit) -> negative = net push backwards
  random_cos_mean                     same cosine for CellOracle's random-permutation control
"""
import os, sys, time, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.spatial import cKDTree
import celloracle as co

t0 = time.time()
def stamp(m): print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

ap = argparse.ArgumentParser()
ap.add_argument("tfs", nargs="*", help="TF symbols to knock out")
ap.add_argument("--tf-file", help="text file, one TF per line")
ap.add_argument("--tag", default=time.strftime("%Y%m%d_%H%M"))
ap.add_argument("--no-ref", action="store_true", help="do not add Srebf2 reference KO")
ap.add_argument("--fast", action="store_true", help="n_neighbors=100, sampled_fraction=0.5")
ap.add_argument("--n-jobs", type=int, default=8)
args = ap.parse_args()

TFS = list(args.tfs)
if args.tf_file:
    TFS += [l.strip() for l in open(args.tf_file) if l.strip() and not l.startswith("#")]
if not args.no_ref and "Srebf2" not in TFS:
    TFS = ["Srebf2"] + TFS
if not TFS:
    sys.exit("no TFs given")

BASE = "."
OUT  = f"{BASE}/phase3_figs/batch_{args.tag}"
os.makedirs(OUT, exist_ok=True)

TARGET = ["Q1", "Q2", "Q3"]
LABEL, FNAME = "Q1→Q3", "Q1Q3"
NN   = 100 if args.fast else 200
FRAC = 0.5 if args.fast else 1.0

# ---- 1. Oracle + GRN (identical to the Srebf2 figures) ----
stamp(f"celloracle {co.__version__}  TFs={TFS}  nn={NN} frac={FRAC}")
oracle = co.load_hdf5(f"{BASE}/r1r9.celloracle.oracle")
links  = co.load_hdf5(f"{BASE}/r1r9.celloracle.links")
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
stamp("GRN fit done")

emb    = oracle.adata.obsm["X_diffmap2d"]
pt     = oracle.adata.obs["palantir_pseudotime"].values
phases = oracle.adata.obs["phase_q"].astype(str).values
in_target = np.isin(phases, TARGET)
phase_colors = {"Q1":"#440154","Q2":"#3b528b","Q3":"#21918c","Q4":"#5ec962","Q5":"#fde725"}
bg_colors = np.array([phase_colors.get(p, "#888") for p in phases])
legend_handles = [Patch(color=c, label=p) for p, c in phase_colors.items()]

# which requested TFs are usable as regulators?
reg_union = set()
for cl, df in links.filtered_links.items():
    reg_union |= set(df["source"].unique())
usable, skipped = [], []
for tf in TFS:
    if tf not in oracle.adata.var_names:
        skipped.append((tf, "not in adata.var_names"))
    elif tf not in reg_union:
        skipped.append((tf, "not a regulator in filtered GRN (no targets) — KO would be a no-op"))
    else:
        usable.append(tf)
for tf, why in skipped:
    stamp(f"SKIP {tf}: {why}")
if not usable:
    sys.exit("no usable TFs")

# ---- 2. Natural rejuvenation flow (KNN pseudotime gradient, per cell) ----
stamp("Natural flow...")
K = 50
tree = cKDTree(emb)
_, idx = tree.query(emb, k=K+1)
idx = idx[:, 1:]
natural_vec = np.zeros_like(emb)
for i in range(len(emb)):
    nbrs = idx[i]
    w = np.clip(pt[nbrs] - pt[i], 0, None)
    if w.sum() < 1e-8:
        continue
    natural_vec[i] = ((emb[nbrs] - emb[i]) * w[:, None]).sum(axis=0) / w.sum()
nat_norm = np.linalg.norm(natural_vec, axis=1)
nat_unit = np.divide(natural_vec, nat_norm[:, None], out=np.zeros_like(natural_vec), where=nat_norm[:, None] > 0)
stamp("Natural flow done")

# ---- 3. grid helper restricted to Q1–Q3 ----
def gridify_q(vec, n_grid=28, min_count=10):
    x_min, x_max = emb[:,0].min(), emb[:,0].max()
    y_min, y_max = emb[:,1].min(), emb[:,1].max()
    pad = 0.05 * max(x_max-x_min, y_max-y_min)
    xs = np.linspace(x_min-pad, x_max+pad, n_grid)
    ys = np.linspace(y_min-pad, y_max+pad, n_grid)
    X, Y = np.meshgrid(xs, ys)
    U = np.zeros_like(X); V = np.zeros_like(X); M = np.zeros_like(X, dtype=bool)
    dx = xs[1]-xs[0]; dy = ys[1]-ys[0]
    for i in range(n_grid):
        for j in range(n_grid):
            mask = ((emb[:,0] >= xs[j]-dx) & (emb[:,0] < xs[j]+dx) &
                    (emb[:,1] >= ys[i]-dy) & (emb[:,1] < ys[i]+dy))
            if mask.sum() < min_count:
                continue
            if np.mean(np.isin(phases[mask], TARGET)) < 0.5:
                continue
            U[i,j] = vec[mask,0].mean(); V[i,j] = vec[mask,1].mean(); M[i,j] = True
    return X, Y, U, V, M

def gridify_all(vec, n_grid=22, min_count=15):
    x_min, x_max = emb[:,0].min(), emb[:,0].max()
    y_min, y_max = emb[:,1].min(), emb[:,1].max()
    pad = 0.05 * max(x_max-x_min, y_max-y_min)
    xs = np.linspace(x_min-pad, x_max+pad, n_grid)
    ys = np.linspace(y_min-pad, y_max+pad, n_grid)
    X, Y = np.meshgrid(xs, ys)
    U = np.zeros_like(X); V = np.zeros_like(X); M = np.zeros_like(X, dtype=bool)
    dx = xs[1]-xs[0]; dy = ys[1]-ys[0]
    for i in range(n_grid):
        for j in range(n_grid):
            mask = ((emb[:,0] >= xs[j]-dx) & (emb[:,0] < xs[j]+dx) &
                    (emb[:,1] >= ys[i]-dy) & (emb[:,1] < ys[i]+dy))
            if mask.sum() >= min_count:
                U[i,j] = vec[mask,0].mean(); V[i,j] = vec[mask,1].mean(); M[i,j] = True
    return X, Y, U, V, M

def plot_panel(ax, X, Y, U, V, M, title, color, zoom=True):
    ax.scatter(emb[:,0], emb[:,1], c="#cccccc", s=4, alpha=0.5, edgecolors="none")
    ax.scatter(emb[in_target,0], emb[in_target,1], c="#888888", s=4, alpha=0.7, edgecolors="none")
    msk = M.flatten()
    ax.quiver(X.flatten()[msk], Y.flatten()[msk], U.flatten()[msk], V.flatten()[msk],
              color=color, angles="xy", scale_units="xy", width=0.006, headwidth=4.5, headlength=5)
    ax.set_aspect("equal"); ax.set_xlabel("DC1"); ax.set_ylabel("DC2")
    ax.set_title(title, fontsize=13, fontweight="bold")
    if zoom:
        q = emb[in_target]; pad = 0.001
        ax.set_xlim(q[:,0].min()-pad, q[:,0].max()+pad)
        ax.set_ylim(q[:,1].min()-pad, q[:,1].max()+pad)

def grid_cos(U1, V1, U2, V2, M):
    m = M.flatten()
    u1, v1, u2, v2 = U1.flatten()[m], V1.flatten()[m], U2.flatten()[m], V2.flatten()[m]
    n1 = np.hypot(u1, v1); n2 = np.hypot(u2, v2)
    ok = (n1 > 1e-12) & (n2 > 1e-12)
    c = (u1[ok]*u2[ok] + v1[ok]*v2[ok]) / (n1[ok]*n2[ok])
    return c

def cell_scores(ko_vec):
    ko_norm = np.linalg.norm(ko_vec, axis=1)
    ok = in_target & (ko_norm > 1e-12) & (nat_norm > 1e-12)
    cos = (ko_vec[ok] * nat_unit[ok]).sum(1) / ko_norm[ok]
    proj = (ko_vec[ok] * nat_unit[ok]).sum(1)
    return dict(n_cells=int(ok.sum()),
                cell_cos_mean=float(cos.mean()), cell_cos_median=float(np.median(cos)),
                cell_frac_opposite=float((cos < 0).mean()),
                ko_mag_mean=float(ko_norm[in_target].mean()),
                proj_on_natural_mean=float(proj.mean()))

Xn, Yn, Un, Vn, Mn = gridify_q(natural_vec)

# ---- 4. KO loop ----
rows, panels = [], [("Natural", Xn, Yn, Un, Vn, Mn, "black")]
for tf in usable:
    stamp(f"=== {tf} KO ===")
    oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=3, ignore_warning=True)
    oracle.estimate_transition_prob(n_neighbors=NN, knn_random=True, sampled_fraction=FRAC, n_jobs=args.n_jobs)
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    ko_vec = np.asarray(oracle.delta_embedding).copy()
    rnd_vec = getattr(oracle, "delta_embedding_random", None)

    sc_ = cell_scores(ko_vec); sc_["TF"] = tf
    Xk, Yk, Uk, Vk, Mk = gridify_q(ko_vec)
    gc = grid_cos(Un, Vn, Uk, Vk, Mn & Mk)
    sc_["grid_cos_mean"] = float(gc.mean()) if len(gc) else np.nan
    sc_["grid_frac_opposite"] = float((gc < 0).mean()) if len(gc) else np.nan
    sc_["grid_n"] = int(len(gc))
    if rnd_vec is not None:
        r = cell_scores(np.asarray(rnd_vec)); sc_["random_cos_mean"] = r["cell_cos_mean"]
    rows.append(sc_)
    panels.append((f"{tf} KO", Xk, Yk, Uk, Vk, Mk, "darkred"))
    stamp(f"  {tf}: cell_cos={sc_['cell_cos_mean']:+.3f}  frac_opp={sc_['cell_frac_opposite']:.2f}  "
          f"grid_cos={sc_['grid_cos_mean']:+.3f}  |KO|={sc_['ko_mag_mean']:.2e}")

    # per-TF zoomed natural vs KO
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    plot_panel(axes[0], Xn, Yn, Un, Vn, Mn, f"Natural EGM rejuvenation\n({LABEL})", "black")
    plot_panel(axes[1], Xk, Yk, Uk, Vk, Mk, f"{tf} KO simulation\n({LABEL})", "darkred")
    fig.suptitle(f"{tf} KO vs natural — {LABEL}   (cell cos = {sc_['cell_cos_mean']:+.2f}, "
                 f"{sc_['cell_frac_opposite']*100:.0f}% opposite)", fontsize=14, fontweight="bold", y=1.00)
    fig.tight_layout()
    fig.savefig(f"{OUT}/{tf}_{FNAME}_natural_vs_KO.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT}/{tf}_{FNAME}_natural_vs_KO.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # per-TF full view
    Xa, Ya, Ua, Va, Ma = gridify_all(ko_vec)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(emb[:,0], emb[:,1], c=bg_colors, s=2, alpha=0.35, edgecolors="none")
    m = Ma.flatten()
    ax.quiver(Xa.flatten()[m], Ya.flatten()[m], Ua.flatten()[m], Va.flatten()[m],
              color="darkred", angles="xy", scale_units="xy", width=0.005, headwidth=4)
    ax.set_aspect("equal"); ax.set_xlabel("DC1"); ax.set_ylabel("DC2")
    ax.set_title(f"{tf} KO simulation flow", fontsize=14, fontweight="bold")
    fig.legend(handles=legend_handles, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.02), frameon=False)
    fig.tight_layout()
    fig.savefig(f"{OUT}/{tf}_KO_fullview.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # save scores incrementally (safe if the run is interrupted)
    pd.DataFrame(rows).to_csv(f"{OUT}/ko_alignment_scores.csv", index=False)

# ---- 5. combined panel ----
ncol = 3
nrow = int(np.ceil(len(panels)/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(6*ncol, 6*nrow))
axes = np.array(axes).flatten()
for ax, (title, X, Y, U, V, M, col) in zip(axes, panels):
    plot_panel(ax, X, Y, U, V, M, title, col, zoom=True)
for ax in axes[len(panels):]:
    ax.axis("off")
fig.suptitle(f"TF KO panel — {LABEL}  (tag={args.tag})", fontsize=16, fontweight="bold", y=1.00)
fig.tight_layout()
fig.savefig(f"{OUT}/batch_{args.tag}_panel.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/batch_{args.tag}_panel.png", dpi=200, bbox_inches="tight")
plt.close(fig)

df = pd.DataFrame(rows).sort_values("cell_cos_mean")
df.to_csv(f"{OUT}/ko_alignment_scores.csv", index=False)
print("\n=== KO alignment summary (sorted: most opposite to rejuvenation first) ===")
print(df[["TF","n_cells","cell_cos_mean","cell_frac_opposite","grid_cos_mean","ko_mag_mean","proj_on_natural_mean"]
        + (["random_cos_mean"] if "random_cos_mean" in df.columns else [])].to_string(index=False))
stamp(f"DONE — outputs in {OUT}")
