"""hdWGCNA on EGM 12-sample dataset (Y1, R0-R9, A1)."""
import os, time, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import seaborn as sns
from py_hdWGCNA import HDWGCNA

H5  = '../senepy_pipeline/responder_analysis/EGM_annotated_responder.h5ad'
OUT = '.'
FIG = os.path.join(OUT, 'figures'); os.makedirs(FIG, exist_ok=True)
WGN = 'EGM_12sample'
GROUP_ORDER = ['Y1','R0','R1','R2','R3','R4','R5','R6','R7','R8','R9','A1']

print("Loading data...", flush=True); t0 = time.time()
adata = sc.read_h5ad(H5)
print(f"  shape: {adata.shape}  ({time.time()-t0:.1f}s)", flush=True)

np.random.seed(42)
keep = []
for g in GROUP_ORDER:
    idx = np.where(adata.obs['sample'].values == g)[0]
    keep.extend(np.random.choice(idx, min(len(idx), 3000), replace=False))
adata = adata[sorted(keep)].copy()
print(f"  subsampled: {adata.n_obs} cells", flush=True)

sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
adata.obs['cell_type'] = pd.Categorical(['fibroblast']*adata.n_obs)
adata.obs['sample'] = adata.obs['sample'].astype('category').cat.reorder_categories(GROUP_ORDER)

# Need a PCA for metacell KNN
print("Computing PCA for metacell KNN...", flush=True); t = time.time()
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, use_highly_variable=True)
print(f"  done ({time.time()-t:.1f}s)", flush=True)

# After scale, .X is scaled; hdWGCNA wants log-normalised back. Reload counts? No — use layers
# Actually hdWGCNA reads from .X via layer='data'. Save log-norm in .layers and restore .X
# Simpler: redo pipeline storing log-norm in .layers['data'] BEFORE scale
print("Re-preparing layers (log-norm in .layers['data'])...", flush=True)
adata2 = sc.read_h5ad(H5)
adata2 = adata2[adata.obs_names].copy()
sc.pp.normalize_total(adata2, target_sum=1e4); sc.pp.log1p(adata2)
adata2.layers['data'] = adata2.X.copy()
adata2.obsm['X_pca'] = adata.obsm['X_pca']
adata2.obs = adata.obs.copy()
adata2.var = adata.var.copy()
adata = adata2

print("\nInitialising HDWGCNA...", flush=True)
hdw = HDWGCNA(adata)

print(f"\n[1/7] setup_for_wgcna...", flush=True); t = time.time()
hdw.setup_for_wgcna(gene_select='fraction', fraction=0.05, wgcna_name=WGN)
print(f"  done ({time.time()-t:.1f}s)", flush=True)

print(f"\n[2/7] metacells_by_groups (k=25)...", flush=True); t = time.time()
hdw.metacells_by_groups(group_by=['cell_type','sample'],
                        reduction='X_pca',
                        k=25, max_shared=10, min_cells=50,
                        target_metacells=1000, wgcna_name=WGN)
print(f"  done ({time.time()-t:.1f}s)", flush=True)

print(f"\n[3/7] normalize_metacells...", flush=True); t = time.time()
try:
    hdw.normalize_metacells(wgcna_name=WGN)
except TypeError:
    hdw.normalize_metacells()
print(f"  done ({time.time()-t:.1f}s)", flush=True)

print(f"\n[4/7] set_dat_expr...", flush=True); t = time.time()
hdw.set_dat_expr(group_name='fibroblast', group_by='cell_type', wgcna_name=WGN)
print(f"  done ({time.time()-t:.1f}s)", flush=True)

print(f"\n[5/7] test_soft_powers (signed)...", flush=True); t = time.time()
hdw.test_soft_powers(network_type='signed', wgcna_name=WGN)
ptab = hdw.adata.uns['hdWGCNA'][WGN]['power_table']
ptab.to_csv(f'{OUT}/power_table.csv', index=False)
sp = hdw.adata.uns['hdWGCNA'][WGN].get('softPower', None)
print(f"  power picked = {sp}  ({time.time()-t:.1f}s)", flush=True)

print(f"\n[6/7] construct_network...", flush=True); t = time.time()
hdw.construct_network(wgcna_name=WGN)
mods_df = hdw.adata.uns['hdWGCNA'][WGN]['modules_df']
mods_df.to_csv(f'{OUT}/modules_df.csv', index=False)
print(f"  modules: {mods_df['module'].nunique()}  ({time.time()-t:.1f}s)")
print("  module sizes:"); print(mods_df['module'].value_counts())

print(f"\n[7/7] module_eigengenes (Harmony on sample)...", flush=True); t = time.time()
hdw.module_eigengenes(group_by_vars='sample', wgcna_name=WGN)
hMEs = hdw.adata.uns['hdWGCNA'][WGN]['hMEs']
print(f"  hMEs shape: {hMEs.shape}  ({time.time()-t:.1f}s)", flush=True)
hMEs.to_csv(f'{OUT}/hMEs_cells.csv')

print("\n[+] module_connectivity...", flush=True); t = time.time()
hdw.module_connectivity(group_by='cell_type', group_name='fibroblast', wgcna_name=WGN)
print(f"  done ({time.time()-t:.1f}s)", flush=True)

# Sample × module matrix
print("\n=== Sample × module mean hME ===", flush=True)
sample_vec = adata.obs['sample'].astype(str).values
me_df = hMEs.copy(); me_df['sample'] = sample_vec
sample_x_mod = me_df.groupby('sample').mean().loc[GROUP_ORDER]
sample_x_mod.to_csv(f'{OUT}/sample_by_module.csv')
print(sample_x_mod.round(3))

mod_cols = [c for c in sample_x_mod.columns if str(c).lower() != 'grey']
mat = sample_x_mod[mod_cols].values

# Clustered heatmap
zmat = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-9)
zdf  = pd.DataFrame(zmat, index=GROUP_ORDER, columns=mod_cols)
g = sns.clustermap(zdf, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                   figsize=(max(8, 0.6*len(mod_cols)), 6),
                   row_cluster=True, col_cluster=True,
                   cbar_kws={'label':'z-scored ME'})
g.ax_heatmap.set_xlabel('Module', fontsize=12, fontweight='bold')
g.ax_heatmap.set_ylabel('Sample', fontsize=12, fontweight='bold')
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
g.savefig(f'{FIG}/sample_x_module_clustermap.png', dpi=300, bbox_inches='tight')
g.savefig(f'{FIG}/sample_x_module_clustermap.pdf', bbox_inches='tight')
plt.close()

# Sample dendrogram
D = pdist(mat, metric='correlation')
Z = linkage(D, method='average')
fig, ax = plt.subplots(figsize=(10, 4))
dendrogram(Z, labels=GROUP_ORDER, ax=ax, leaf_font_size=14,
           color_threshold=0.4, above_threshold_color='#666')
ax.set_ylabel('1 − Pearson r\n(across module ME)', fontsize=12, fontweight='bold')
ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
ax.tick_params(labelsize=12)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
fig.savefig(f'{FIG}/sample_dendrogram_by_modules.png', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(f'{FIG}/sample_dendrogram_by_modules.pdf', bbox_inches='tight', facecolor='white')
plt.close()

# Soft-power diagnostic
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
slope_col = 'slope' if 'slope' in ptab.columns else ptab.columns[2]
sft_col   = 'SFT.R.sq' if 'SFT.R.sq' in ptab.columns else ptab.columns[1]
a1.plot(ptab['Power'], -np.sign(ptab[slope_col])*ptab[sft_col], 'o-', color='steelblue')
a1.axhline(0.85, color='red', ls='--', lw=1)
a1.set_xlabel('soft power'); a1.set_ylabel('signed R²'); a1.set_title('Scale-Free Topology fit')
mk_col = 'mean.k.' if 'mean.k.' in ptab.columns else 'mean_k'
a2.plot(ptab['Power'], ptab[mk_col], 'o-', color='darkorange')
a2.set_xlabel('soft power'); a2.set_ylabel('mean connectivity'); a2.set_title('Mean connectivity')
plt.tight_layout()
fig.savefig(f'{FIG}/soft_power.png', dpi=200, bbox_inches='tight'); plt.close()

print("\nSaving augmented adata...", flush=True)
hdw.adata.write_h5ad(f'{OUT}/adata_with_hdwgcna.h5ad')
print("\n=== ALL DONE ===", flush=True)
