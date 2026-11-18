import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings, gc, os
warnings.filterwarnings('ignore')

PREPROC = "../preprocessed"
OUT = "."

print("Loading merged data...", flush=True)
adata = sc.read_h5ad(os.path.join(PREPROC, "EGM_all_samples_filtered.h5ad"))
print(f"  {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

# Store raw counts
adata.layers['counts'] = adata.X.copy()

# Normalize + log
print("Normalizing...", flush=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG
print("Finding HVGs...", flush=True)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key='sample')
print(f"  HVGs: {adata.var['highly_variable'].sum()}", flush=True)

# Scale + PCA
print("PCA...", flush=True)
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50)
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
del adata_hvg
gc.collect()

# Harmony integration
print("Harmony integration...", flush=True)
import harmonypy
ho = harmonypy.run_harmony(adata.obsm['X_pca'], adata.obs, 'sample', max_iter_harmony=20)
# Auto-detect shape
Z = ho.Z_corr
if Z.shape[0] != adata.n_obs:
    Z = Z.T
adata.obsm['X_pca_harmony'] = Z
print("  Harmony done.", flush=True)

# Neighbors + UMAP
print("Computing neighbors...", flush=True)
sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=30, n_pcs=30)
print("Computing UMAP...", flush=True)
sc.tl.umap(adata, min_dist=0.3)

# Save
umap_path = os.path.join(OUT, "EGM_umap_integrated.h5ad")
adata.write_h5ad(umap_path)
print(f"Saved: {umap_path}", flush=True)

# ── Plot ──
print("Plotting...", flush=True)

SAMPLE_ORDER = ['Y1','A1','R0','R1','R2','R3','R4','R5','R6','R7','R8','R9']
SAMPLE_COLORS = {
    'Y1': '#2ecc71', 'A1': '#e74c3c',
    'R0': '#1a5276', 'R1': '#1f6f8b', 'R2': '#2589a0', 'R3': '#2ca3b5',
    'R4': '#45b7c9', 'R5': '#6bc5d2', 'R6': '#91d3dc', 'R7': '#add9e0',
    'R8': '#c9e4e8', 'R9': '#dff0f2',
}

# Condition color
COND_COLORS = {
    'Young': '#2ecc71', 'Aged': '#e74c3c', 'EGM_treated': '#3498db'
}

umap_coords = adata.obsm['X_umap']

fig_dir = os.path.join(OUT, 'figures')

# --- Plot 1: By sample (12 colors) ---
fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
# Plot in reverse order so Y1/A1 are on top
for sample in reversed(SAMPLE_ORDER):
    mask = adata.obs['sample'].values == sample
    ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
               c=SAMPLE_COLORS[sample], s=0.3, alpha=0.4, rasterized=True, label=sample)

# Legend
handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=SAMPLE_COLORS[s],
            markersize=6, label=s) for s in SAMPLE_ORDER]
leg = ax.legend(handles=handles, fontsize=7, loc='upper right', frameon=True,
                markerscale=1, handletextpad=0.3, labelspacing=0.25, borderpad=0.4,
                edgecolor='#cccccc')
ax.set_xlabel('UMAP 1', fontsize=12, fontweight='bold')
ax.set_ylabel('UMAP 2', fontsize=12, fontweight='bold')
ax.tick_params(labelsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
fig.savefig(os.path.join(fig_dir, 'umap_by_sample.png'), dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.close(fig)
print("  Saved: umap_by_sample.png", flush=True)

# --- Plot 2: By group (Young / Aged / EGM_treated) ---
fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
for grp in ['EGM_treated', 'Aged', 'Young']:
    mask = adata.obs['group'].values == grp
    ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
               c=COND_COLORS[grp], s=0.3, alpha=0.4, rasterized=True, label=grp)

handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=COND_COLORS[g],
            markersize=8, label=g) for g in ['Young','Aged','EGM_treated']]
leg = ax.legend(handles=handles, fontsize=9, loc='upper right', frameon=True,
                markerscale=1, handletextpad=0.4, borderpad=0.5, edgecolor='#cccccc')
ax.set_xlabel('UMAP 1', fontsize=12, fontweight='bold')
ax.set_ylabel('UMAP 2', fontsize=12, fontweight='bold')
ax.tick_params(labelsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
fig.savefig(os.path.join(fig_dir, 'umap_by_group.png'), dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.close(fig)
print("  Saved: umap_by_group.png", flush=True)

# --- Plot 3: Hub 0 SCAP score on UMAP ---
scored_path = os.path.join(OUT, "step4_scores/EGM_all_senepy_scored.h5ad")
scored_obs = sc.read_h5ad(scored_path, backed='r').obs
hub0_scores = scored_obs['senepy_hub0'].values

fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
sc_plot = ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
                      c=hub0_scores, cmap='RdYlGn_r', s=0.3, alpha=0.5,
                      vmin=np.percentile(hub0_scores, 2),
                      vmax=np.percentile(hub0_scores, 98), rasterized=True)
cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label('SCAP-0 Score', fontsize=10, fontweight='bold')
cbar.ax.tick_params(labelsize=8)
ax.set_xlabel('UMAP 1', fontsize=12, fontweight='bold')
ax.set_ylabel('UMAP 2', fontsize=12, fontweight='bold')
ax.tick_params(labelsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
fig.savefig(os.path.join(fig_dir, 'umap_scap0_score.png'), dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.close(fig)
print("  Saved: umap_scap0_score.png", flush=True)

print("\nDone!", flush=True)
