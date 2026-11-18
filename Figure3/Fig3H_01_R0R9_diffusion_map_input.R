PY_SCANPY <- path.expand("~/.conda/envs/scanpy_env/bin/python")   # conda env with anndata

# ---- R0-R9 diffusion map --------------------------------------------------------
suppressPackageStartupMessages({
library(Seurat); library(harmony); library(destiny)
library(dplyr); library(tidyr); library(tibble); library(ggplot2)
library(patchwork); library(ggrepel); library(viridisLite)
library(RColorBrewer); library(scales); library(colorspace); library(mgcv)
})
set.seed(42)
base_dir <- "EGM_multiflex"
mod_dir  <- file.path(base_dir, "module_gene_lists")
out_dir  <- file.path(base_dir, "diffmap_figs")
dir.create(out_dir, showWarnings = FALSE)
# ---- 1. Load + re-add module scores ----------------------------------
seu <- readRDS(file.path(base_dir, "seu_ds.rds"))
mods       <- c("M1","M3","M5","M6","M9","M2","M11")
mod_colors <- c(M1="#C0392B", M3="#E67E22", M5="#2A9D8F",
M6="#52B788", M9="#377EB8", M2="#7B5BA6", M11="#C39BD3")
mod_gene_list <- setNames(lapply(mods, function(m)
intersect(read.csv(file.path(mod_dir, paste0(m, "_genes.csv")))$gene_name,
rownames(seu))), mods)
if (!all(mods %in% colnames(seu@meta.data))) {
seu <- AddModuleScore(seu, features = mod_gene_list,
name = "modscore_", nbin = 24, ctrl = 100)
score_cols <- paste0("modscore_", seq_along(mods))
names(seu@meta.data)[match(score_cols, names(seu@meta.data))] <- mods
}
# ---- 2. R0-R9 subset + figure-grade downsample (3000/sample) ---------
r_samples <- paste0("R", 0:9)
seu_r <- subset(seu, subset = sample %in% r_samples)
seu_r$sample <- droplevels(factor(seu_r$sample, levels = r_samples))
n_per <- 3000          # 메모리 부족하면 2500 또는 2000으로
keep <- unlist(lapply(levels(seu_r$sample), function(s) {
cells <- colnames(seu_r)[seu_r$sample == s]
sample(cells, size = min(n_per, length(cells)))
}))
seu_r <- subset(seu_r, cells = keep)
cat("Cells / sample:\n"); print(table(seu_r$sample))
# ---- 3. Re-normalize → PCA → Harmony(batch only) ---------------------
DefaultAssay(seu_r) <- "RNA"
seu_r <- NormalizeData(seu_r) |>
FindVariableFeatures(nfeatures = 3000) |>
ScaleData(vars.to.regress = "percent.mt", verbose = FALSE) |>
RunPCA(npcs = 50, verbose = FALSE)
seu_r <- RunHarmony(seu_r,
group.by.vars  = "batch",     # sample이 아닌 batch!
reduction.use  = "pca",
dims.use       = 1:50,
reduction.save = "harmony_batch",
plot_convergence = FALSE)
# ---- 4. F-stat gene selection ---------------------------------------
expr_hvg  <- as.matrix(GetAssayData(seu_r, layer = "data")[VariableFeatures(seu_r), ])
sample_f  <- seu_r$sample
fstat <- apply(expr_hvg, 1, function(x) {
m <- tapply(x, sample_f, mean)
v <- tapply(x, sample_f, var)
(var(m) * (length(m) - 1)) / (sum(v) / length(v))
})
top_genes <- names(sort(fstat, decreasing = TRUE))[1:2000]
saveRDS(top_genes, file.path(out_dir, "fstat_top2000.rds"))
# ---- 5. PCA on F-stat genes → DiffusionMap --------------------------
seu_r <- ScaleData(seu_r, features = top_genes, verbose = FALSE)
seu_r <- RunPCA(seu_r, features = top_genes, npcs = 50,
reduction.name = "pca_sample", verbose = FALSE)
embed_s <- Embeddings(seu_r, "pca_sample")[, 1:50]
t0 <- Sys.time()
dm <- DiffusionMap(embed_s, k = 30, sigma = "local",
density_norm = TRUE, n_eigs = 20, verbose = TRUE)
cat("DM done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
dc_mat <- as.matrix(dm@eigenvectors[, 1:5])
colnames(dc_mat) <- paste0("DC", 1:5)
rownames(dc_mat) <- colnames(seu_r)
seu_r[["diffmap"]] <- CreateDimReducObject(
embeddings = dc_mat, key = "DC_", assay = "RNA"
)
# ---- 6. DPT ---------------------------------------------------------
dc1      <- dc_mat[, 1]
r0_cells <- which(seu_r$sample == "R0")
tip      <- r0_cells[which.max(abs(dc1[r0_cells] - median(dc1)))]
dpt      <- DPT(dm, tips = tip)
seu_r$dpt <- dpt$dpt
if (mean(seu_r$dpt[seu_r$sample == "R0"]) >
mean(seu_r$dpt[seu_r$sample == "R9"])) {
seu_r$dpt <- max(seu_r$dpt) - seu_r$dpt
}

# ---- export seu_r -> h5ad -------------------------------------------------------
# ============================================================
# Export seu_r (R0-R9 only) → seu_R0R9_for_UNAGI.h5ad
# Raw counts preserved in both X and layers["counts"]
# ============================================================
suppressPackageStartupMessages({
library(Seurat); library(Matrix); library(reticulate)
})
base_dir <- "EGM_multiflex"
# Python (scanpy_env) 연결
Sys.setenv(RETICULATE_PYTHON = PY_SCANPY)
use_python(PY_SCANPY, required = TRUE)
ad    <- import("anndata",      convert = FALSE)
scipy <- import("scipy.sparse", convert = FALSE)
np    <- import("numpy",        convert = FALSE)
# ---- 1. seu_r 준비 (v5 → v3 다운캐스트 안전장치) -----------------
DefaultAssay(seu_r) <- "RNA"
if (inherits(seu_r[["RNA"]], "Assay5")) {
seu_r <- JoinLayers(seu_r)
seu_r[["RNA"]] <- as(seu_r[["RNA"]], "Assay")
}
# ---- 2. sparse(genes × cells) → scipy CSR(cells × genes) --------
to_csr <- function(m) {
m <- as(m, "CsparseMatrix")
csc <- scipy$csc_matrix(
reticulate::tuple(m@x, as.integer(m@i), as.integer(m@p)),
shape = reticulate::tuple(nrow(m), ncol(m))
)
csc$T$tocsr()
}
counts <- GetAssayData(seu_r, assay = "RNA", layer = "counts")
data_  <- GetAssayData(seu_r, assay = "RNA", layer = "data")
# 확인 — counts는 정수 max, data는 소수 max
cat("counts max:", max(counts), " (정수여야 OK)\n")
cat("data   max:", round(max(data_), 3), " (소수)\n")
X    <- to_csr(counts)        # ← raw counts
logX <- to_csr(data_)         # ← log-normalized
# ---- 3. obs / var -----------------------------------------------
obs <- seu_r@meta.data
# factor → character (anndata 호환)
obs[] <- lapply(obs, function(x) if (is.factor(x)) as.character(x) else x)
# UNAGI에서 쓰기 좋게 stage(int) 추가
obs$stage <- as.integer(gsub("R", "", obs$sample))
var <- data.frame(row.names = rownames(seu_r))
# ---- 4. embeddings → obsm ---------------------------------------
obsm_list <- list()
for (rn in c("pca", "harmony", "harmony_batch", "pca_sample",
"diffmap", "umap_harmony")) {
if (rn %in% Reductions(seu_r)) {
obsm_list[[paste0("X_", rn)]] <- Embeddings(seu_r, rn)
}
}
# ---- 5. AnnData 생성 & 저장 -------------------------------------
adata <- ad$AnnData(
X      = X,
obs    = obs,
var    = var,
layers = reticulate::dict(counts    = X,         # ← raw counts (X 미러)
logcounts = logX),      # ← log-norm
obsm   = do.call(reticulate::dict, obsm_list)
)
out_path <- file.path(base_dir, "seu_R0R9_for_oracle_with_raw.h5ad")
adata$write_h5ad(out_path, compression = "gzip")
