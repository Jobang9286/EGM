PY_SCANPY <- path.expand("~/.conda/envs/scanpy_env/bin/python")   # conda env with anndata

# oracle_input_hvg.h5ad → Seurat .rds (V3 assay 강제 버전)
# 모든 obs, obsm, var, layers 보존
# ============================================================
suppressPackageStartupMessages({
library(Seurat); library(Matrix); library(reticulate)
})
# ---- 0. V3 assay 강제 (핵심) ------------------------------------
options(Seurat.object.assay.version = "v3")
# ---- 1. Python 연결 + h5ad 로드 ---------------------------------
Sys.setenv(RETICULATE_PYTHON = PY_SCANPY)
use_python(PY_SCANPY, required = TRUE)
ad    <- import("anndata",      convert = FALSE)
scipy <- import("scipy.sparse", convert = FALSE)
np    <- import("numpy",        convert = FALSE)
h5ad_path <- "EGM_multiflex/oracle/oracle_input_hvg.h5ad"
adata <- ad$read_h5ad(h5ad_path)
# ---- 2. helper: scipy CSR → R dgCMatrix --------------------------
to_dgC <- function(spX) {
if (py_to_r(scipy$issparse(spX))) {
csc <- spX$tocsc()
Matrix::sparseMatrix(
i      = as.integer(py_to_r(csc$indices)) + 1L,
p      = as.integer(py_to_r(csc$indptr)),
x      = as.numeric(py_to_r(csc$data)),
dims   = as.integer(unlist(py_to_r(csc$shape))),
index1 = TRUE
)
} else {
as(py_to_r(spX), "CsparseMatrix")
}
}
# ---- 3. counts / logcounts (genes × cells) ----------------------
counts_py <- adata$layers["raw_count"]
log_py    <- adata$layers["logcounts"]
counts <- t(to_dgC(counts_py))
logX   <- t(to_dgC(log_py))
cell_ids <- py_to_r(adata$obs_names$values)
gene_ids <- py_to_r(adata$var_names$values)
rownames(counts) <- rownames(logX) <- gene_ids
colnames(counts) <- colnames(logX) <- cell_ids
cat("counts:", dim(counts), "max:", max(counts), "\n")
cat("logX  :", dim(logX),   "max:", round(max(logX), 3), "\n")
# ---- 4. obs / var -----------------------------------------------
obs <- py_to_r(adata$obs)
rownames(obs) <- cell_ids
# 적절한 컬럼 factor 복원
factor_cols <- c("sample","batch","group","orig.ident",
"seurat_clusters","RNA_snn_res.0.5",
paste0("dm_res", c(0.1,0.2,0.3,0.5,0.8)),
paste0("hb_res", c(0.1,0.2,0.3,0.45,0.5,0.55,0.6,0.7,0.8)),
"local_sample","phase","cluster6","time_group")
for (cc in intersect(factor_cols, colnames(obs))) {
obs[[cc]] <- factor(obs[[cc]])
}
if ("sample" %in% colnames(obs))
obs$sample <- factor(obs$sample, levels = paste0("R", 0:9))
# ---- 5. Seurat 객체 + data layer --------------------------------
seu <- CreateSeuratObject(counts = counts, meta.data = obs)
seu <- SetAssayData(seu, layer = "data", new.data = logX)
# var info → feature metadata
var_df <- py_to_r(adata$var)
rownames(var_df) <- gene_ids
for (cc in colnames(var_df)) {
seu[["RNA"]][[cc]] <- var_df[[cc]][match(rownames(seu), rownames(var_df))]
}
# ---- 6. obsm → Reductions ---------------------------------------
obsm_keys <- py_to_r(adata$obsm_keys())
cat("\nLoading obsm:\n"); print(obsm_keys)

add_reduction <- function(seu, obsm_key, red_name, key_prefix) {
emb <- as.matrix(py_to_r(adata$obsm[obsm_key]))
storage.mode(emb) <- "double"        # ← bool/int 모두 numeric으로 강제
rownames(emb) <- cell_ids
colnames(emb) <- paste0(key_prefix, seq_len(ncol(emb)))
seu[[red_name]] <- CreateDimReducObject(
embeddings = emb, key = key_prefix, assay = "RNA"
)
seu
}
red_map <- list(
list("X_pca",                       "pca",            "PC_"),
list("X_pca_sample",                "pca_sample",     "PCs_"),
list("X_harmony",                   "harmony",        "harm_"),
list("X_harmony_batch",             "harmony_batch",  "harmb_"),
list("X_umap_harmony",              "umap_harmony",   "UMAP_"),
list("X_diffmap",                   "diffmap",        "DC_"),
list("X_diffmap2d",                 "diffmap2d",      "DC2d_"),
list("DM_EigenVectors",             "dm_eigen",       "DCe_"),
list("DM_EigenVectors_multiscaled", "dm_eigen_ms",    "DCms_"),
list("palantir_fate_probabilities", "palantir_fate",  "FP_"),
list("branch_masks",                "branch_masks",   "BR_")
)
for (rm in red_map) {
if (rm[[1]] %in% obsm_keys) {
seu <- add_reduction(seu, rm[[1]], rm[[2]], rm[[3]])
cat("  +", rm[[1]], "→", rm[[2]], "\n")
}
}
# ---- 7. 확인 -----------------------------------------------------
print(seu)
cat("\nReductions:\n"); print(Reductions(seu))
cat("\nMeta cols:\n"); print(colnames(seu@meta.data))
cat("\nSamples:\n"); print(table(seu$sample))
cat("\nPalantir pt range:\n"); print(range(seu$palantir_pseudotime))
# ---- 8. 저장 -----------------------------------------------------
out_path <- "EGM_multiflex/oracle/oracle_input_hvg.rds"
saveRDS(seu, out_path)

# ============================================================
# Pseudotime-based phase clustering on diffusion map
# ============================================================
suppressPackageStartupMessages({
library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
library(patchwork); library(RColorBrewer); library(viridisLite); library(mclust)
})
options(Seurat.object.assay.version = "v3")
set.seed(42)
base_dir <- "EGM_multiflex/oracle"
seu <- readRDS(file.path(base_dir, "oracle_input_hvg.rds"))
pt   <- seu$palantir_pseudotime
mods <- c("M1","M3","M5","M6","M9","M2","M11")
# ---- A. Equal-quantile binning -----------------------------------
n_phase <- 5     # 3, 5, 7 비교하고 싶으면 바꿔서 다시
seu$phase_q <- cut(pt,
breaks  = quantile(pt, seq(0, 1, length = n_phase + 1)),
include.lowest = TRUE,
labels  = paste0("Q", 1:n_phase))

# ---- phase = phase_q, save ------------------------------------------------------
seu$phase <- seu$phase_q
saveRDS(seu, file.path(base_dir, "oracle_input_hvg.rds"))

# ============================================================
# Export updated Seurat .rds → h5ad for CellOracle
# 모든 obs (phase_q/gmm/cp 포함) + obsm + layers 보존
# ============================================================
suppressPackageStartupMessages({
library(Seurat); library(Matrix); library(reticulate)
})
options(Seurat.object.assay.version = "v3")
base_dir <- "EGM_multiflex/oracle"
Sys.setenv(RETICULATE_PYTHON = PY_SCANPY)
use_python(PY_SCANPY, required = TRUE)
ad    <- import("anndata",      convert = FALSE)
scipy <- import("scipy.sparse", convert = FALSE)
np    <- import("numpy",        convert = FALSE)
# ---- 1. rds 로드 ------------------------------------------------
seu <- readRDS(file.path(base_dir, "oracle_input_hvg.rds"))
cat("Loaded:", ncol(seu), "cells ×", nrow(seu), "genes\n")
cat("Reductions:\n"); print(Reductions(seu))
# ---- 2. helper: dgCMatrix → scipy CSR (cells × genes) -----------
to_csr <- function(m) {
m <- as(m, "CsparseMatrix")
csc <- scipy$csc_matrix(
reticulate::tuple(m@x, as.integer(m@i), as.integer(m@p)),
shape = reticulate::tuple(nrow(m), ncol(m))
)
csc$T$tocsr()
}
# ---- 3. X / layers ----------------------------------------------
counts <- GetAssayData(seu, layer = "counts")
data_  <- GetAssayData(seu, layer = "data")
X_py    <- to_csr(counts)
logX_py <- to_csr(data_)
# ---- 4. obs (factor → character) --------------------------------
obs <- seu@meta.data
obs[] <- lapply(obs, function(x) if (is.factor(x)) as.character(x) else x)
# ---- 5. var -----------------------------------------------------
feat_meta <- seu[["RNA"]][[]]                # var info
var_df    <- if (nrow(feat_meta) == nrow(seu)) feat_meta else data.frame(row.names = rownames(seu))
# ---- 6. obsm ----------------------------------------------------
obsm_list <- list()
for (rn in Reductions(seu)) {
obsm_list[[paste0("X_", rn)]] <- Embeddings(seu, rn)
}
# ---- 7. AnnData & 저장 ------------------------------------------
adata <- ad$AnnData(
X      = X_py,
obs    = obs,
var    = var_df,
layers = reticulate::dict(raw_count = X_py,
logcounts = logX_py),
obsm   = do.call(reticulate::dict, obsm_list)
)
out_path <- file.path(base_dir, "oracle_input_hvg_with_phase.h5ad")
adata$write_h5ad(out_path, compression = "gzip")
print(adata)
cat(sprintf("\nSaved: %s (%.2f GB)\n",
out_path, file.info(out_path)$size / 1e9))
gc()
