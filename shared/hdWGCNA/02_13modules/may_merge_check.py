"""Why are M4 and M12 missing from the May 2026 hdWGCNA (hdwgcna_pipeline)? Uses only what that run saved.

Re-cuts May's stored TOM with May's parameters to get the 13 pre-merge modules, then applies py_hdWGCNA's
merge exactly as installed (and, for comparison, with the ME-row lookup that skips grey). The installed merge must
reproduce May's final module_labels gene for gene.
"""
import inspect, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, h5py
import py_hdWGCNA.network as net
from py_hdWGCNA.utils import compute_module_eigengenes, hclust_and_cut

g = h5py.File('../hdwgcna_pipeline/adata_with_hdwgcna.h5ad', 'r')['uns/hdWGCNA/EGM_12sample']
tom   = g['tom_dissim'][()]
dat   = g['dat_expr'][()].astype(np.float64)
genes = np.array([x.decode() for x in g['dat_expr_genes'][()]])
may   = g['module_labels'][()]

pre = hclust_and_cut(tom, deepSplit=4, pamRespectsDendro=True, pamStage=False, minClusterSize=50,
                         cutHeight=0.995, method='average')
print("pre-merge modules:", len(set(pre)) - (0 in set(pre)), "| grey genes:", int((pre == 0).sum()))
MEs, _, valid = compute_module_eigengenes(dat, pre)
print("ME rows (first-appearance order):", ['M%d' % v for v in valid])

installed = net._merge_via_python
src = inspect.getsource(installed)
old = "label_to_me_index = {lbl: i for i, lbl in enumerate(_fa_labels) if lbl != 0}"
exec(compile(src.replace(old, "label_to_me_index = {lbl: i for i, lbl in enumerate([m for m in _fa_labels if m != 0])}"),
             net.__file__, 'exec'), net.__dict__)
fixed = net._merge_via_python

lab_inst  = installed(MEs, pre, 0.2, dat)[1]
lab_fixed = fixed(MEs, pre, 0.2, dat)[1]
print("\ninstalled merge == May module_labels:", bool(np.array_equal(lab_inst, may)), int((lab_inst == may).sum()), "/", len(may))

me_r = np.corrcoef(MEs)
name = {v: i for i, v in enumerate(valid)}
def show(lab, title):
    print(f"\n--- {title}")
    for old_m in sorted(set(pre) - {0}):
        new_m = np.unique(lab[pre == old_m])[0]
        if new_m != old_m:
            print(f"  M{old_m} ({int((pre == old_m).sum())} genes) -> M{new_m} ({int((pre == new_m).sum())} genes)"
                  f" | ME r = {me_r[name[old_m], name[new_m]]:.3f}")
show(lab_inst, "merges applied in May (installed py_hdWGCNA)")
show(lab_fixed, "merges the 0.2 cut actually calls for (ME r > 0.8)")
print("\nME r > 0.8 pairs among pre-merge modules:")
for i in range(len(valid)):
    for j in range(i + 1, len(valid)):
        if me_r[i, j] > 0.8: print(f"  M{valid[i]} - M{valid[j]}: r = {me_r[i, j]:.3f}")
pd.DataFrame({'gene_name': genes, 'premerge': ['M%d' % l if l else 'grey' for l in pre],
              'may_final': ['M%d' % l if l else 'grey' for l in may]}).to_csv('may_premerge_vs_final.csv', index=False)
