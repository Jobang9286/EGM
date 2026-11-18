"""Cell-level module eigengenes (Harmony on sample) for the May 2026 hdWGCNA network, with
  current   = May's final modules (merge as run in May: M4 -> M3, M12 -> M2), as a check against May's saved hMEs
  corrected = May's pre-merge modules with the merge the 0.2 cut calls for (M9 -> M5, M13 -> M6; M4, M12 kept)
  unmerged  = all 13 pre-merge modules
on the same 36,000 cells and the same function and settings as run_hdwgcna.py. Prints only small summary tables.
Run from hdwgcna_fixmerge/ after may_merge_check.py, in conda env senepy with PYTHONPATH=pyoverlay
(numpy 2.0.2, as in the May run): python corrected_hmes.py unmerged   (or: corrected)
"""
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import py_hdWGCNA.modules as M
from py_hdWGCNA.utils import get_hdWGCNA_data

WGN = 'EGM_12sample'
adata = sc.read_h5ad('../hdwgcna_pipeline/adata_with_hdwgcna.h5ad')
wd = get_hdWGCNA_data(adata, WGN)
may_hmes = wd['hMEs'].copy()
genes = wd['modules_df']['gene_name'].values
may_mod = pd.Series(wd['modules_df']['module'].values, index=genes)
pre = pd.read_csv('may_premerge_vs_final.csv').set_index('gene_name')
assert (pre.loc[genes, 'may_final'].values == may_mod.values).all()
import sys
sets = {'current': may_mod, 'corrected': pre.loc[genes, 'premerge'].replace({'M9': 'M5', 'M13': 'M6'}),
        'unmerged': pre.loc[genes, 'premerge']}   # all 13 pre-merge modules, no merging
sets = {k: sets[k] for k in (sys.argv[1:] or sets)}

samples = adata.obs['sample'].astype(str).values
order = ['A1', 'Y1'] + ['R%d' % i for i in range(10)]
for tag, mods in sets.items():
    wd['modules_df'] = pd.DataFrame({'gene_name': genes, 'module': mods.values, 'color': 'grey'})
    a = M.module_eigengenes(adata, group_by_vars='sample', harmonize=True, wgcna_name=WGN)
    h = get_hdWGCNA_data(a, WGN)['hMEs'].copy()
    h.to_csv(f'hMEs_cells_{tag}.csv')
    if tag == 'current':
        common = [c for c in h.columns if c in may_hmes.columns]
        d = np.abs(h[common].values - may_hmes.loc[h.index, common].values).max()
        print(f"[check] current vs May saved hMEs: modules {len(common)}/{may_hmes.shape[1]}, max |diff| = {d:.2e}")
    mods_sorted = sorted([c for c in h.columns if c != 'grey'], key=lambda m: int(m[1:]))
    mean = h.groupby(samples).mean().loc[order, mods_sorted]
    p = pd.DataFrame(np.nan, index=mods_sorted, columns=order)
    for m in mods_sorted:
        ref = h.loc[samples == 'A1', m].values
        for s in order[1:]:
            p.loc[m, s] = mannwhitneyu(h.loc[samples == s, m].values, ref, alternative='two-sided',
                                      use_continuity=True, method='asymptotic').pvalue
    print(f"\n### {tag} module sizes:", mods.value_counts().drop('grey').reindex(mods_sorted).to_dict())
    print(f"### {tag} mean hME (samples x modules)"); print(mean.to_csv(float_format='%.6g'))
    print(f"### {tag} Mann-Whitney p vs A1 (modules x samples, unadjusted)"); print(p.to_csv(float_format='%.6g'))
for m in (() if 'unmerged' in sets else ('M4', 'M12')):
    print(f"\n### premerge {m} genes:", ','.join(sorted(pre.index[pre['premerge'] == m])))
print("\n=== DONE")
