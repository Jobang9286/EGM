"""
CellOracle Phase 1 — maximized resource usage
- 128-core Xeon, 187GB RAM
- n_jobs = 48 (leave ~80 cores for OS / other)
- KNN k auto, but cap n_jobs for bagging across clusters
"""
import sys, time, os
os.environ["OPENBLAS_NUM_THREADS"] = "2"   # avoid oversubscription within each worker
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["OMP_NUM_THREADS"] = "2"

import numpy as np
import scanpy as sc
import celloracle as co

t0 = time.time()
def stamp(msg):
    print(f"[{time.time()-t0:6.1f}s] {msg}", flush=True)

BASE = "."
IN   = f"{BASE}/oracle_input_hvg_plus_TFs.h5ad"
OUT_ORACLE = f"{BASE}/r1r9.celloracle.oracle"
OUT_LINKS  = f"{BASE}/r1r9.celloracle.links"

N_JOBS = 48   # 48 workers × 2 BLAS threads = 96 cores (out of 128)

stamp(f"celloracle {co.__version__}, scanpy {sc.__version__}, N_JOBS={N_JOBS}")

# ---- 1. Load + verify ------------------------------------------
adata = sc.read_h5ad(IN)
stamp(f"Loaded: {adata.shape}")

adata.X = adata.layers["raw_count"].copy()
stamp(f"X set to raw_count (max={adata.X.max()}, dtype={adata.X.dtype})")

adata.obs["phase_q"] = adata.obs["phase_q"].astype(str)
print("phase_q value counts:")
print(adata.obs["phase_q"].value_counts().sort_index())
stamp(f"X_diffmap2d shape: {adata.obsm['X_diffmap2d'].shape}")


# ---- 2. Oracle 객체 ------------------------------------------
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name="phase_q",
    embedding_name="X_diffmap2d",
)
stamp("Oracle.import_anndata_as_raw_count done")

# ---- 3. PCA + KNN imputation ---------------------------------
oracle.perform_PCA()
ev = np.cumsum(oracle.pca.explained_variance_ratio_)
n_pcs = int(np.searchsorted(ev, 0.5)) + 1
n_pcs = min(max(n_pcs, 20), 50)
stamp(f"PCA done, n_pcs={n_pcs} (cum_var={ev[n_pcs-1]:.3f})")

n_cell = adata.n_obs
k_knn  = int(0.025 * n_cell)
stamp(f"KNN imputation: k={k_knn}, n_pca_dims={n_pcs}, n_jobs={N_JOBS}")
oracle.knn_imputation(
    n_pca_dims=n_pcs,
    k=k_knn,
    balanced=True,
    b_sight=k_knn * 8,
    b_maxl=k_knn * 4,
    n_jobs=N_JOBS,
)
stamp("KNN imputation done")


# ---- 4. Base GRN 결합 ---------------------------------------
stamp("Loading mouse base GRN...")
base_GRN = co.data.load_mouse_scATAC_atlas_base_GRN()
oracle.import_TF_data(TF_info_matrix=base_GRN)
stamp("Base GRN imported")


# ---- 5. GRN 추론 ---------------------------------------------
stamp(f"Starting get_links with n_jobs={N_JOBS}...")
links = oracle.get_links(
    cluster_name_for_GRN_unit="phase_q",
    alpha=10,
    model_method="bagging_ridge",
    verbose_level=1,
    n_jobs=N_JOBS,
)
stamp(f"get_links DONE. clusters: {list(links.links_dict.keys())}")
for cl, df in links.links_dict.items():
    print(f"  cluster {cl}: {len(df)} edges")


# ---- 6. 저장 ------------------------------------------------
oracle.to_hdf5(OUT_ORACLE)
links.to_hdf5(file_path=OUT_LINKS)
stamp(f"Saved:\n  {OUT_ORACLE}\n  {OUT_LINKS}")
stamp(f"TOTAL TIME: {time.time()-t0:.1f}s")
