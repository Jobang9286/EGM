# ---------------------------------------------------------------------------
# A. Load, QC, normalise, Harmony, cluster; remove outlier cluster 12
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
library(Seurat)
library(harmony)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(future)
})
plan("multicore", workers = 128)
options(future.globals.maxSize = 192 * 1024^3)  # 50 GB
# ---- 2. 샘플 메타 ----------------------------------------------
#  batch : CellRanger multi 다중화 배치 (경로 해석 + 기술적 batch)
#  group : 생물학적 그룹
base_dir <- "EGM_multiflex/20260226_Multiflex_DA250141/CellRanger_decode"
sample_info <- tribble(
~sample, ~batch,
"A1",    "Set3",
"Y1",    "Set1",
"R0",    "Set1",
"R1",    "Set1",
"R2",    "Set1",
"R3",    "Set2",
"R4",    "Set2",
"R5",    "Set2",
"R6",    "Set2",
"R7",    "Set3",
"R8",    "Set3",
"R9",    "Set3"
) %>%
mutate(
path  = file.path(base_dir, batch, "per_sample_outs",
sample, "count", "sample_filtered_feature_bc_matrix"),
group = case_when(
sample == "A1"      ~ "Aged",
sample == "Y1"      ~ "Young",
grepl("^R", sample) ~ "Treat",
TRUE                ~ "Other"
)
)
stopifnot(all(dir.exists(sample_info$path)))
# ---- 3. 샘플별 Seurat 객체 생성 --------------------------------
seu_list <- lapply(seq_len(nrow(sample_info)), function(i) {
s <- sample_info$sample[i]
message("Loading: ", s)
mtx <- Read10X(sample_info$path[i])
obj <- CreateSeuratObject(
counts       = mtx,
project      = s,
min.cells    = 3,
min.features = 200
)
obj$sample <- s
obj$batch  <- sample_info$batch[i]
obj$group  <- sample_info$group[i]
obj <- RenameCells(obj, add.cell.id = s)     # 바코드 충돌 방지
obj
})
names(seu_list) <- sample_info$sample
# ---- 4. Merge ---------------------------------------------------
seu <- merge(seu_list[[1]], y = seu_list[-1], project = "EGM_multiflex")
# Seurat v5: count layer 통합 (있다면)
if (inherits(seu[["RNA"]], "Assay5")) {
seu <- JoinLayers(seu)
}
# species 확인 후 mito pattern 조정:  human "^MT-"   mouse "^mt-"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
# ---- 6. QC 필터링 (임계값은 분포 보고 조정) --------------------
seu <- subset(seu,
subset = nFeature_RNA >= 500 &
nFeature_RNA <= 7500 &
percent.mt   <  10
)
# ---- 7. Normalize → HVG → Scale → PCA ---------------------------
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu, vars.to.regress = c("percent.mt"))
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
# ---- 8. Harmony 통합 -------------------------------------------
seu <- RunHarmony(
seu,
group.by.vars    = c("sample"),   # 샘플 + 다중화 batch 동시 보정
reduction.use    = "pca",
dims.use         = 1:20,
reduction.save   = "harmony",
plot_convergence = FALSE
)
# ---- 9. Harmony embedding 기반 downstream ----------------------
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:15,
reduction.name = "umap_harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:15)
seu <- FindClusters(seu, resolution = 0.5)
# ---- 10. 저장 (선택) -------------------------------------------
saveRDS(seu, "EGM_multiflex/seu_harmony.rds")
seu2 <- subset(seu, idents = "12", invert = TRUE)
# 기존 harmony reduction을 그대로 쓰고 그래프/UMAP만 다시
seu2 <- RunUMAP(seu2, reduction = "harmony", dims = 1:20,
reduction.name = "umap_harmony")
seu2 <- FindNeighbors(seu2, reduction = "harmony", dims = 1:20)
seu2 <- FindClusters(seu2, resolution = 0.5)
saveRDS(seu2, "EGM_multiflex/seu_harmony_rmov12.rds")

# ---------------------------------------------------------------------------
# B1. First 5,000-cell-per-sample draw -> seu_ds.rds
#     The v5 -> v3 assay downcast ran before the save, so seu_ds.rds holds a v3 "Assay".
# ---------------------------------------------------------------------------
seu2 <- readRDS("EGM_multiflex/seu_harmony_rmov12.rds")
library(Seurat)
set.seed(42)
n_per <- 5000
# 샘플별 셀 바코드에서 최대 n_per개씩 추출
keep_cells <- unlist(lapply(unique(seu2$sample), function(s) {
cells <- colnames(seu2)[seu2$sample == s]
sample(cells, size = min(n_per, length(cells)))
}))
seu2_ds <- subset(seu2, cells = keep_cells)
# ---- 2. seu2_ds v5 → v3 다운캐스트 ---------------------------
DefaultAssay(seu2_ds) <- "RNA"
if (inherits(seu2_ds[["RNA"]], "Assay5")) {
seu2_ds <- JoinLayers(seu2_ds)
seu2_ds[["RNA"]] <- as(seu2_ds[["RNA"]], "Assay")
}
saveRDS(seu2_ds, "EGM_multiflex/seu_ds.rds")

# ---------------------------------------------------------------------------
# B2. Re-cluster at res 0.2, then the same draw again
# ---------------------------------------------------------------------------
seu2 <- readRDS("EGM_multiflex/seu_harmony_rmov12.rds")
seu2 <- FindClusters(seu2, resolution = 0.2)
library(Seurat)
set.seed(42)
n_per <- 5000
# 샘플별 셀 바코드에서 최대 n_per개씩 추출
keep_cells <- unlist(lapply(unique(seu2$sample), function(s) {
cells <- colnames(seu2)[seu2$sample == s]
sample(cells, size = min(n_per, length(cells)))
}))
seu2_ds <- subset(seu2, cells = keep_cells)

# ---------------------------------------------------------------------------
# C. Cluster re-annotation (res 0.2: 3+4 merged, renumbered 1-6), sample labels A1->A, Y1->Y,
#    5 time groups
# ---------------------------------------------------------------------------
# ----- 0. 클러스터/샘플 reannotation --------------------------
# 0~6 → 0~5 (3,4 합치고 5→4, 6→5)
remap <- c("0"="0","1"="1","2"="2","3"="3","4"="3","5"="4","6"="5")
seu2_ds$seurat_clusters <- factor(
unname(remap[as.character(seu2_ds$RNA_snn_res.0.2)]),
levels = c("0","1","2","3","4","5")
)
# 0~5 → 1~6
seu2_ds$seurat_clusters <- factor(
as.integer(as.character(seu2_ds$seurat_clusters)) + 1,
levels = 1:6
)
Idents(seu2_ds) <- "seurat_clusters"
# sample 라벨 (A1→A, Y1→Y)
seu2_ds$sample <- factor(
dplyr::recode(as.character(seu2_ds$sample), "A1"="A", "Y1"="Y"),
levels = c("A", paste0("R", 0:9), "Y")
)
# ----- 0. sample_group 메타 추가 ------------------------------
seu2_ds$sample_group <- dplyr::case_when(
seu2_ds$sample %in% c("A","R0")        ~ "Aged+R0",
seu2_ds$sample %in% paste0("R", 1:3)   ~ "R1-3",
seu2_ds$sample %in% paste0("R", 4:6)   ~ "R4-6",
seu2_ds$sample %in% paste0("R", 7:9)   ~ "R7-9",
seu2_ds$sample == "Y"                  ~ "Young"
)
seu2_ds$sample_group <- factor(seu2_ds$sample_group,
levels = c("Aged+R0","R1-3","R4-6","R7-9","Young"))

# ---------------------------------------------------------------------------
# D1. Module scores (May 2026 11-module hdWGCNA lists) + 6 time groups
# ---------------------------------------------------------------------------
library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
mod_dir  <- "EGM_multiflex/module_gene_lists"
mods_use <- c("M1","M3","M5","M6","M9","M2","M11")
mod_gene_list <- setNames(lapply(mods_use, function(m)
intersect(read.csv(file.path(mod_dir, paste0(m, "_genes.csv")))$gene_name,
rownames(seu2_ds))), mods_use)
seu2_ds <- AddModuleScore(seu2_ds, features = mod_gene_list,
name = "modscore_", nbin = 24, ctrl = 100)
# Seurat가 modscore_1, _2, ... 식으로 붙임 → 이름 변경
score_cols <- paste0("modscore_", seq_along(mods_use))
names(seu2_ds@meta.data)[match(score_cols, names(seu2_ds@meta.data))] <- mods_use
# group_split: A, R0, R1-3, R4-6, R7-9, Young
seu2_ds$sample_group6 <- dplyr::case_when(
seu2_ds$sample == "A"                 ~ "Aged",
seu2_ds$sample == "R0"                ~ "R0",
seu2_ds$sample %in% paste0("R", 1:3)  ~ "R1-3",
seu2_ds$sample %in% paste0("R", 4:6)  ~ "R4-6",
seu2_ds$sample %in% paste0("R", 7:9)  ~ "R7-9",
seu2_ds$sample == "Y"                 ~ "Young"
)
seu2_ds$sample_group6 <- factor(seu2_ds$sample_group6,
levels = c("Aged","R0","R1-3","R4-6","R7-9","Young"))

# ---------------------------------------------------------------------------
# D2. Supp Fig 8H: UMAP FeaturePlot per module
# ---------------------------------------------------------------------------
library(Seurat); library(ggplot2); library(patchwork); library(dplyr); library(tibble)
# ----- 1. 모듈별 cutoff + 색 설정 ------------------------------
mod_cfg <- tribble(
~module, ~color,     ~min_cut, ~max_cut,
# Aged/Somatic
"M1",    "#C0392B",  "q05",    "q95",   # brick red
"M3",    "#E67E22",  "q05",    "q95",   # orange
# Young/Early
"M5",    "#2A9D8F",  "q05",    "q95",   # teal
"M6",    "#52B788",  "q05",    "q95",   # green
"M9",    "#377EB8",  "q05",    "q95",   # blue
# Young/Mid
"M2",    "#7B5BA6",  "q05",    "q95",   # purple
"M11",   "#C39BD3",  "q05",    "q95"    # lavender
)
# 카테고리 매핑 (행 단위)
mod_cat <- list(
"Aged/Somatic" = c("M1","M3"),
"Young/Early"  = c("M5","M6","M9"),
"Young/Mid"    = c("M2","M11")
)
# ----- 2. panel 생성 함수 -------------------------------------
make_panel <- function(mod) {
cfg <- mod_cfg %>% filter(module == mod)
FeaturePlot(seu2_ds,
features    = mod,
reduction   = "umap_harmony",
order       = TRUE,
min.cutoff  = cfg$min_cut,
max.cutoff  = cfg$max_cut,
pt.size     = 0.4) +
scale_color_gradient(low = "grey92", high = cfg$color) +
labs(title = mod, x = NULL, y = NULL) +
theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
axis.text  = element_blank(),
axis.ticks = element_blank(),
axis.title = element_blank())
}
# ----- 3. 행 단위 조립 (3-col 폭 고정) -------------------------
build_row <- function(mods, n_cols = 3) {
plots <- lapply(mods, make_panel)
while (length(plots) < n_cols) plots <- c(plots, list(plot_spacer()))
Reduce(`|`, plots)
}
p_aged  <- build_row(mod_cat[["Aged/Somatic"]])
p_early <- build_row(mod_cat[["Young/Early"]])
p_mid   <- build_row(mod_cat[["Young/Mid"]])
# ----- 4. 최종 plot ------------------------------------------
final <- p_aged / p_early / p_mid +
plot_layout(heights = c(1, 1, 1))
final
# 특정 module cutoff 조정
mod_cfg$max_cut[mod_cfg$module == "M5"]  <- "q99"   # 진한 값 살리기
mod_cfg$min_cut[mod_cfg$module == "M11"] <- "q10"   # 약한 노이즈 컷
mod_cfg$max_cut[mod_cfg$module == "M3"]  <- 0.3     # 절대값도 가능
mod_cfg$max_cut[mod_cfg$module %in% c("M5","M6","M9")] <- "q99"
# 다시 그리기
final <- build_row(mod_cat[["Aged/Somatic"]]) /
build_row(mod_cat[["Young/Early"]])  /
build_row(mod_cat[["Young/Mid"]])
final

# ---------------------------------------------------------------------------
# D3. Supp Fig 8I: mean module score per cluster across time groups
# ---------------------------------------------------------------------------
library(ggplot2); library(dplyr); library(tidyr); library(RColorBrewer)
mods_use <- c("M1","M3","M5","M6","M9","M2","M11")
# ----- 1. cluster × group × module 평균 + SEM -----------------
bar_df <- seu2_ds@meta.data %>%
dplyr::select(cluster = seurat_clusters,
group   = sample_group6,
all_of(mods_use)) %>%
pivot_longer(cols = all_of(mods_use),
names_to = "module", values_to = "score") %>%
group_by(cluster, group, module) %>%
summarise(mean_score = mean(score),
sem        = sd(score)/sqrt(n()),
.groups = "drop") %>%
mutate(module  = factor(module,  levels = mods_use),
cluster = factor(cluster, levels = 1:6),
group   = factor(group,
levels = c("Aged","R0","R1-3","R4-6","R7-9","Young")))
# 모듈 하나 그리는 함수
make_panel <- function(mod_name) {
ggplot(bar_df %>% filter(module == mod_name),
aes(x = as.numeric(group), y = mean_score,
color = cluster, group = cluster)) +
geom_hline(yintercept = 0, color = "grey60",
linewidth = 0.3, linetype = "dashed") +
geom_smooth(method = "loess", span = 1.2, se = FALSE, linewidth = 1) +
scale_x_continuous(breaks = 1:6, labels = levels(bar_df$group)) +
scale_color_manual(values = brewer.pal(6, "Set2"), name = "Cluster") +
labs(x = NULL, y = "Mean module score", title = mod_name) +
theme_bw(base_size = 10) +
theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
panel.grid.minor = element_blank())
}
# 행별 3-col 폭에 맞춰 spacer로 정렬
row1 <- (make_panel("M1")  | make_panel("M3")  | plot_spacer())
row2 <- (make_panel("M5")  | make_panel("M6")  | make_panel("M9"))
row3 <- (make_panel("M2")  | make_panel("M11") | plot_spacer())
final <- (row1 / row2 / row3) + plot_layout(guides = "collect")
final
