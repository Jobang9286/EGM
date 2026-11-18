# ============================================================
# 01_processing.R
# EGM_multiflex scRNA-seq: Load → QC → Harmony → Cluster → Reannotation
# Output : seu2 (post-filter, reannotated)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(dplyr); library(tibble)
  library(ggplot2); library(patchwork); library(future)
  library(RColorBrewer); library(viridisLite)
})

plan("multicore", workers = 16)
options(future.globals.maxSize = 192 * 1024^3)


# ---- 1. 샘플 메타 + CellRanger 경로 ---------------------------
base_dir <- "EGM_multiflex/20260226_Multiflex_DA250141/CellRanger_decode"

# 디스크 폴더명은 "A1","Y1" 이지만 객체에는 처음부터 "A","Y" 로 저장
sample_info <- tribble(
  ~sample, ~dir_name, ~batch,
  "A",  "A1", "Set3",
  "Y",  "Y1", "Set1",
  "R0", "R0", "Set1",
  "R1", "R1", "Set1",
  "R2", "R2", "Set1",
  "R3", "R3", "Set2",
  "R4", "R4", "Set2",
  "R5", "R5", "Set2",
  "R6", "R6", "Set2",
  "R7", "R7", "Set3",
  "R8", "R8", "Set3",
  "R9", "R9", "Set3"
) %>%
  mutate(
    path  = file.path(base_dir, batch, "per_sample_outs",
                      dir_name, "count", "sample_filtered_feature_bc_matrix"),
    group = case_when(
      sample == "A"        ~ "Aged",
      sample == "Y"        ~ "Young",
      grepl("^R", sample)  ~ "Treat"
    )
  )
stopifnot(all(dir.exists(sample_info$path)))


# ---- 2. Per-sample Seurat 객체 → merge ------------------------
seu_list <- lapply(seq_len(nrow(sample_info)), function(i) {
  s <- sample_info$sample[i]
  message("Loading: ", s)
  mtx <- Read10X(sample_info$path[i])
  obj <- CreateSeuratObject(mtx, project = s,
                            min.cells = 3, min.features = 200)
  obj$sample <- s
  obj$batch  <- sample_info$batch[i]
  obj$group  <- sample_info$group[i]
  RenameCells(obj, add.cell.id = s)
})
names(seu_list) <- sample_info$sample

seu <- merge(seu_list[[1]], y = seu_list[-1], project = "EGM_multiflex")
if (inherits(seu[["RNA"]], "Assay5")) seu <- JoinLayers(seu)


# ---- 3. QC metric + VlnPlot with thresholds -------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")   # mouse
seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

p <- VlnPlot(seu, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
             group.by = "sample", pt.size = 0, combine = FALSE)
p[[1]] <- p[[1]] + geom_hline(yintercept = c(500, 7500),
                              color = "red", linetype = "dashed", linewidth = 0.6)
p[[3]] <- p[[3]] + geom_hline(yintercept = 10,
                              color = "red", linetype = "dashed", linewidth = 0.6)
wrap_plots(p, ncol = 1)


# ---- 4. QC 필터 -----------------------------------------------
seu <- subset(seu, subset = nFeature_RNA >= 500 &
                            nFeature_RNA <= 7500 &
                            percent.mt   <  10)


# ---- 5. Normalize → HVG → Scale → PCA → Harmony ---------------
seu <- NormalizeData(seu) |>
       FindVariableFeatures(selection.method = "vst", nfeatures = 3000) |>
       ScaleData(vars.to.regress = "percent.mt") |>
       RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(seu, ndims = 50)

seu <- RunHarmony(seu,
                  group.by.vars    = "sample",
                  reduction.use    = "pca",
                  dims.use         = 1:20,
                  reduction.save   = "harmony",
                  plot_convergence = FALSE)


# ---- 6. 1차 cluster (res 0.5) → outlier cluster 제거 ----------
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:15,
               reduction.name = "umap_harmony") |>
       FindNeighbors(reduction = "harmony", dims = 1:15) |>
       FindClusters(resolution = 0.5)

seu2 <- subset(seu, idents = "12", invert = TRUE)


# ---- 7. 2차 cluster (res 0.2, seu2 기준) ----------------------
seu2 <- RunUMAP(seu2, reduction = "harmony", dims = 1:20,
                reduction.name = "umap_harmony") |>
        FindNeighbors(reduction = "harmony", dims = 1:20) |>
        FindClusters(resolution = 0.2)


# ---- 8. Reannotation ------------------------------------------
# (a) 0~6 → 0~5 (3,4 합치고 5→4, 6→5)
remap <- c("0"="0","1"="1","2"="2","3"="3","4"="3","5"="4","6"="5")
seu2$seurat_clusters <- factor(
  unname(remap[as.character(seu2$RNA_snn_res.0.2)]),
  levels = c("0","1","2","3","4","5")
)
# (b) 0~5 → 1~6
seu2$seurat_clusters <- factor(
  as.integer(as.character(seu2$seurat_clusters)) + 1,
  levels = 1:6
)
Idents(seu2) <- "seurat_clusters"

# (c) sample factor 순서 고정 (A, R0~R9, Y)
seu2$sample <- factor(seu2$sample,
                      levels = c("A", paste0("R", 0:9), "Y"))

# (d) 그룹 메타 컬럼들
seu2$sample_group <- dplyr::case_when(
  seu2$sample %in% c("A","R0")        ~ "Aged+R0",
  seu2$sample %in% paste0("R", 1:3)   ~ "R1-3",
  seu2$sample %in% paste0("R", 4:6)   ~ "R4-6",
  seu2$sample %in% paste0("R", 7:9)   ~ "R7-9",
  seu2$sample == "Y"                  ~ "Young"
) |>
  factor(levels = c("Aged+R0","R1-3","R4-6","R7-9","Young"))

seu2$sample_group6 <- dplyr::case_when(
  seu2$sample == "A"                  ~ "Aged",
  seu2$sample == "R0"                 ~ "R0",
  seu2$sample %in% paste0("R", 1:3)   ~ "R1-3",
  seu2$sample %in% paste0("R", 4:6)   ~ "R4-6",
  seu2$sample %in% paste0("R", 7:9)   ~ "R7-9",
  seu2$sample == "Y"                  ~ "Young"
) |>
  factor(levels = c("Aged","R0","R1-3","R4-6","R7-9","Young"))


# ---- 9. 기본 시각화 -------------------------------------------
pal_clust <- brewer.pal(6, "Set2")

DimPlot(seu2, reduction = "umap_harmony", label = FALSE,
        cols = pal_clust, pt.size = 1.3)

DimPlot(seu2, reduction = "umap_harmony", label = FALSE,
        cols = pal_clust, pt.size = 1.2, split.by = "sample_group6", ncol = 3)

# Cluster × sample stacked bar (A/Y 흰 점선)
sample_order <- c("A", paste0("R",0:9), "Y")
df <- as.data.frame(table(Cluster = seu2$seurat_clusters,
                          Sample  = seu2$sample)) %>%
  group_by(Cluster) %>% mutate(Fraction = Freq/sum(Freq)*100) %>% ungroup()
df$Sample <- factor(df$Sample, levels = sample_order)

pal_sample <- c(A = "gray",
                setNames(viridis(10), paste0("R", 0:9)),
                Y = "#dee")
border_col <- setNames(ifelse(sample_order %in% c("A","Y"), "white", NA), sample_order)
border_lty <- setNames(ifelse(sample_order %in% c("A","Y"), "dashed","solid"), sample_order)

ggplot(df, aes(x = Cluster, y = Fraction,
               fill = Sample, color = Sample, linetype = Sample)) +
  geom_col(width = 0.85, linewidth = 0.5) +
  scale_fill_manual(values = pal_sample, breaks = sample_order) +
  scale_color_manual(values = border_col, breaks = sample_order, guide = "none") +
  scale_linetype_manual(values = border_lty, breaks = sample_order, guide = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
  labs(x = "Clusters", y = "Sample Fraction (%)") +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"))


# ---- 10. 저장 -------------------------------------------------
out_dir <- "EGM_multiflex"
saveRDS(seu2, file.path(out_dir, "seu2.rds"))
qs::qsave(seu2, file.path(out_dir, "seu2.qs"), nthreads = 8)

cat("=== DONE: 01_processing ===\n")
