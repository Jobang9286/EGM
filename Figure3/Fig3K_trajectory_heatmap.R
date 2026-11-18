PY_CELLORACLE <- path.expand("~/.conda/envs/celloracle_env/bin/python")   # env with anndata

# ---- setup ---------------------------------------------------------------------
suppressPackageStartupMessages({
library(Seurat); library(dplyr); library(tidyr); library(reticulate)
library(ComplexHeatmap); library(circlize); library(RColorBrewer)
})
options(Seurat.object.assay.version = "v3")
base_dir <- "EGM_multiflex/oracle"
mod_dir  <- "EGM_multiflex/1_module_heatmap_trajectory_script/module_gene_lists"
seu <- readRDS(file.path(base_dir, "oracle_input_hvg.rds"))
# ----- 1) Python + full-gene h5ad -----
Sys.setenv(RETICULATE_PYTHON = PY_CELLORACLE)
use_python(PY_CELLORACLE, required = TRUE)
ad <- import("anndata", convert = FALSE)
adata_full <- ad$read_h5ad(file.path(base_dir, "seu_R0R9_for_oracle_with_raw.h5ad"))
all_genes_full <- py_to_r(adata_full$var_names$values)
cat("Full gene count:", length(all_genes_full), "\n")

# ---- module gene sets (metadata stored in trajectory_env) ----------------------
mods_use <- c("M1","M6","M5","M9","M2","M11")    # M3 빼고 M6, M9 추가
TOP_N    <- 80                              # 동일 크기 권장
must_include <- list(
M5 = c("Srebf2","Hmgcr","Insig1","Akt1","Acly","Fasn",
"Fdps","Dhcr24","Msmo1","Ldlr")
)
mod_colors <- c(M1="#C0392B", M6="#52B788", M5="#2A9D8F",
M9="#377EB8", M2="#7B5BA6", M11="#C39BD3")
mod_bio <- c(
M1  = "Aged / ECM senescence",
M6  = "Proteostasis & vesicle trafficking",
M5  = "Cholesterol & lipid biosynthesis",
M9  = "Cytoskeleton / actin",
M2  = "Young identity (mid)",
M11 = "Cell cycle / mitosis"
)
# ----- 2. mod_genes -----
mod_genes <- lapply(mods_use, function(m) {
d <- read.csv(file.path(mod_dir, paste0(m, "_genes.csv")))
gene_col <- intersect(c("gene_name","gene","symbol","Gene"), colnames(d))[1]
if (is.na(gene_col)) gene_col <- colnames(d)[1]
score_col <- grep("kME|score|weight", colnames(d), ignore.case=TRUE, value=TRUE)[1]
if (!is.na(score_col)) d <- d[order(-d[[score_col]]), ]
top_genes <- head(d[[gene_col]], TOP_N)
if (!is.null(must_include[[m]])) {
forced <- intersect(must_include[[m]], d[[gene_col]])
top_genes <- unique(c(top_genes, forced))
}
intersect(top_genes, all_genes_full)
})
names(mod_genes) <- mods_use
print(sapply(mod_genes, length))
# ----- 3. Expression + binning -----
all_genes <- unique(unlist(mod_genes))
adata_sub <- adata_full[, all_genes]
logX <- py_to_r(adata_sub$layers["logcounts"]$toarray())
rownames(logX) <- py_to_r(adata_sub$obs_names$values)
colnames(logX) <- all_genes
expr <- t(logX[colnames(seu), ])
N_BIN <- 80
cell_bin <- cut(seu$palantir_pseudotime,
breaks = quantile(seu$palantir_pseudotime,
probs = seq(0,1,length.out=N_BIN+1)),
include.lowest = TRUE, labels = FALSE)
bin_mean <- sapply(seq_len(N_BIN), function(b)
rowMeans(expr[, cell_bin == b, drop = FALSE]))
rownames(bin_mean) <- rownames(expr)
z_mat <- t(scale(t(bin_mean)))
z_mat[z_mat >  1.5] <-  1.5
z_mat[z_mat < -1.5] <- -1.5
# ----- 4. Gene ordering -----
gene_module <- rep(names(mod_genes), times=sapply(mod_genes, length))
names(gene_module) <- unlist(mod_genes)
gene_module <- gene_module[rownames(z_mat)]
peak_bin <- apply(z_mat, 1, which.max)
gene_order <- data.frame(gene=rownames(z_mat), module=gene_module, peak=peak_bin) %>%
mutate(module=factor(module, levels=mods_use)) %>%
arrange(module, peak) %>% pull(gene)
z_mat <- z_mat[gene_order, ]

# ---- pseudotime-dynamic genes ---------------------------------------------------
suppressPackageStartupMessages({
library(Seurat); library(dplyr); library(reticulate)
library(ComplexHeatmap); library(circlize); library(RColorBrewer)
})
# ----- 1. Full-gene expression matrix (R로 가져오기) -----
adata_sub_all <- adata_full   # 16k gene
logX_all <- py_to_r(adata_sub_all$layers["logcounts"]$toarray())
rownames(logX_all) <- py_to_r(adata_sub_all$obs_names$values)
colnames(logX_all) <- all_genes_full
# HVG로 후보 좁히기 (속도 + noise 줄임)
hvg <- VariableFeatures(seu)              # 3000 HVG
if (length(hvg) == 0) {
# 없으면 variance top 3000 직접
v <- apply(logX_all[colnames(seu), ], 2, var)
hvg <- names(sort(v, decreasing = TRUE))[1:3000]
}
hvg_in <- intersect(hvg, all_genes_full)
expr <- t(logX_all[colnames(seu), hvg_in])     # genes × cells
cat("Candidate genes:", nrow(expr), "\n")
# ----- 2. Quantile bin + 평균 -----
N_BIN <- 80
cell_bin <- cut(seu$palantir_pseudotime,
breaks = quantile(seu$palantir_pseudotime,
probs = seq(0,1,length.out=N_BIN+1)),
include.lowest = TRUE, labels = FALSE)
bin_mean <- sapply(seq_len(N_BIN), function(b)
rowMeans(expr[, cell_bin == b, drop = FALSE]))
rownames(bin_mean) <- rownames(expr)
# ----- 3. Dynamic gene 선별 (amplitude 기반) -----
# 각 gene의 smooth (rolling mean) → max - min
smooth_window <- 5
bin_smooth <- t(apply(bin_mean, 1, function(x) {
zoo::rollmean(x, k = smooth_window, fill = NA, align = "center")
}))
# NA는 원래 값으로
bin_smooth[is.na(bin_smooth)] <- bin_mean[is.na(bin_smooth)]
amp  <- apply(bin_smooth, 1, function(x) max(x) - min(x))
peak <- apply(bin_smooth, 1, which.max)
mean_expr <- rowMeans(expr)
# 충분히 발현 + 충분히 dynamic만 keep
keep <- amp > quantile(amp, 0.6) & mean_expr > 0.1
dyn_genes <- rownames(bin_mean)[keep]
cat("Dynamic genes selected:", length(dyn_genes), "\n")

# ---- k-means phases + heatmap (Fig 3K left) -------------------------------------
suppressPackageStartupMessages({
library(ComplexHeatmap); library(circlize); library(RColorBrewer)
library(clusterProfiler); library(org.Mm.eg.db)
library(dplyr); library(ggplot2); library(patchwork)
})
# ----- 1. Dynamic gene을 trajectory-shape으로 k-means 클러스터링 -----
# bin_smooth (이전 단계의 smoothed bin matrix) 사용 — dynamic gene만 추림
dyn_smooth <- bin_smooth[dyn_genes, ]      # genes × bins
dyn_z <- t(scale(t(dyn_smooth)))           # row-wise z
# k-means: 4 phase
set.seed(42)
K <- 4
km <- kmeans(dyn_z, centers = K, nstart = 25, iter.max = 100)
gene_phase <- factor(paste0("P", km$cluster), levels = paste0("P", 1:K))
# Cluster를 peak 시간 순으로 재정렬 (P1=early → P4=late)
phase_peak <- sapply(seq_len(K), function(k) {
mean_z <- colMeans(dyn_z[km$cluster == k, , drop = FALSE])
which.max(mean_z)
})
phase_order <- order(phase_peak)
new_label <- paste0("P", seq_len(K))[match(seq_len(K), phase_order)]
gene_phase <- factor(new_label[km$cluster], levels = paste0("P", 1:K))
# 각 phase의 gene 수
print(table(gene_phase))
# ----- 2. Heatmap (middle-focused, phase-colored sidebar) -----
gene_order <- names(sort(setNames(km$cluster, rownames(dyn_z)),
decreasing = FALSE))
# 정렬: phase → 그 안에서 peak 위치
peaks <- apply(dyn_z, 1, which.max)
gene_order <- data.frame(gene = rownames(dyn_z),
phase = gene_phase[rownames(dyn_z)],
peak  = peaks) %>%
arrange(phase, peak) %>% pull(gene)
z_mat <- dyn_z[gene_order, ]
z_mat[z_mat >  1.5] <-  1.5
z_mat[z_mat < -1.5] <- -1.5
# 끝 trim
trim_n <- 2
z_mat <- z_mat[, seq(trim_n+1, N_BIN-trim_n)]
phase_colors <- setNames(brewer.pal(K, "Set2"), paste0("P", 1:K))
col_fun <- colorRamp2(seq(-1.5, 1.5, length.out=11),
rev(brewer.pal(11, "Spectral")))
ht <- Heatmap(
z_mat, name="z-score", col=col_fun,
cluster_rows=FALSE, cluster_columns=FALSE,
show_row_names=FALSE, show_column_names=FALSE,
left_annotation = rowAnnotation(
Phase = factor(gene_phase[rownames(z_mat)], levels=paste0("P",1:K)),
col = list(Phase = phase_colors),
width = unit(0.4, "cm")
),
row_split = factor(gene_phase[rownames(z_mat)], levels=paste0("P",1:K)),
row_gap = unit(0, "mm"),
row_title_rot = 0,
border = FALSE, use_raster=TRUE, raster_quality=5,
column_title = "Palantir pseudotime",
column_title_side = "bottom",
heatmap_legend_param = list(title="z-score")
)
draw(ht, column_title=sprintf("Pseudotime-dynamic genes — %d phase clusters", K),
column_title_gp=gpar(fontsize=14, fontface="bold"))

# ---- save trajectory_env.rds ----------------------------------------------------
# ----- Save 환경 (RStudio Desktop으로 옮길 핵심만) -----
trajectory_env <- list(
# Heatmap 재현용
z_mat        = z_mat,
gene_module  = gene_module,
gene_order   = gene_order,
mods_use     = mods_use,
mod_colors   = mod_colors,
mod_bio      = mod_bio,
# Module별 gene list
mod_genes    = mod_genes,
# GSEA용 union (M5 + M6)
M5M6_union   = unique(c(mod_genes$M5, mod_genes$M6)),
# 개별
M5_genes     = mod_genes$M5,
M6_genes     = mod_genes$M6,
# 메타 정보 (figure 캡션용)
N_BIN        = N_BIN,
driver_genes = must_include$M5,
# 전체 dynamic gene set (옵션 — GO universe로)
dyn_genes_all = if (exists("dyn_genes")) dyn_genes else NULL
)
out_path <- "EGM_multiflex/oracle/trajectory_env.rds"
saveRDS(trajectory_env, out_path)
cat("Saved:", out_path, "(size:", round(file.info(out_path)$size/1e6, 1), "MB)\n")
cat("Contents:\n"); str(trajectory_env, max.level = 1)
