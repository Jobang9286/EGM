# ============================================================
# 02_downsample_modulescores.R
# Input  : seu2 (from 01_processing.R)
# Steps  :
#   - Downsample to 5000 cells / sample → seu2_ds
#   - AddModuleScore for hdWGCNA modules M1,M3,M5,M6,M9,M2,M11
#   - Visualizations: UMAP FeaturePlot per module (module-specific colors),
#                     line plot of module dynamics per cluster
# Output : seu2_ds (rds + qs)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(RColorBrewer)
})

set.seed(42)

base_dir <- "EGM_multiflex"
mod_dir  <- file.path(base_dir, "module_gene_lists")

# ---- 0. seu2 로드 ---------------------------------------------
seu2 <- readRDS(file.path(base_dir, "seu2.rds"))


# ---- 1. 샘플당 5000 셀 downsample ----------------------------
n_per <- 5000
keep_cells <- unlist(lapply(levels(seu2$sample), function(s) {
  cells <- colnames(seu2)[seu2$sample == s]
  sample(cells, size = min(n_per, length(cells)))
}))
seu2_ds <- subset(seu2, cells = keep_cells)
cat("After downsample:\n"); print(table(seu2_ds$sample))


# ---- 2. 모듈 gene set 로드 + AddModuleScore -------------------
mods_use <- c("M1","M3","M5","M6","M9","M2","M11")

mod_gene_list <- setNames(lapply(mods_use, function(m)
  intersect(read.csv(file.path(mod_dir, paste0(m, "_genes.csv")))$gene_name,
            rownames(seu2_ds))), mods_use)

seu2_ds <- AddModuleScore(seu2_ds,
                          features = mod_gene_list,
                          name     = "modscore_",
                          nbin     = 24,
                          ctrl     = 100)

# AddModuleScore가 modscore_1, _2... 식으로 붙임 → M? 이름으로 변경
score_cols <- paste0("modscore_", seq_along(mods_use))
names(seu2_ds@meta.data)[match(score_cols, names(seu2_ds@meta.data))] <- mods_use


# ---- 3. 모듈별 색 + cutoff 설정 -------------------------------
mod_cfg <- tribble(
  ~module, ~color,     ~min_cut, ~max_cut,
  "M1",    "#C0392B",  "q05",    "q95",   # brick red    (Aged/Somatic)
  "M3",    "#E67E22",  "q05",    "q95",   # orange       (Aged/Somatic)
  "M5",    "#2A9D8F",  "q05",    "q99",   # teal         (Young/Early)
  "M6",    "#52B788",  "q05",    "q99",   # green        (Young/Early)
  "M9",    "#377EB8",  "q05",    "q99",   # blue         (Young/Early)
  "M2",    "#7B5BA6",  "q05",    "q95",   # purple       (Young/Mid)
  "M11",   "#C39BD3",  "q05",    "q95"    # lavender     (Young/Mid)
)

mod_cat <- list(
  "Aged/Somatic" = c("M1","M3"),
  "Young/Early"  = c("M5","M6","M9"),
  "Young/Mid"    = c("M2","M11")
)


# ---- 4. UMAP FeaturePlot per module --------------------------
make_panel <- function(mod) {
  cfg <- mod_cfg %>% filter(module == mod)
  FeaturePlot(seu2_ds, features = mod, reduction = "umap_harmony",
              order = TRUE,
              min.cutoff = cfg$min_cut,
              max.cutoff = cfg$max_cut,
              pt.size = 0.4) +
    scale_color_gradient(low = "grey92", high = cfg$color) +
    labs(title = mod, x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
}

build_row <- function(mods, n_cols = 3) {
  plots <- lapply(mods, make_panel)
  while (length(plots) < n_cols) plots <- c(plots, list(plot_spacer()))
  Reduce(`|`, plots)
}

p_module_umap <- build_row(mod_cat[["Aged/Somatic"]]) /
                 build_row(mod_cat[["Young/Early"]])  /
                 build_row(mod_cat[["Young/Mid"]])
print(p_module_umap)


# ---- 5. Line plot: module × cluster × group -------------------
bar_df <- seu2_ds@meta.data %>%
  dplyr::select(cluster = seurat_clusters,
                group   = sample_group6,
                all_of(mods_use)) %>%
  pivot_longer(cols = all_of(mods_use),
               names_to = "module", values_to = "score") %>%
  group_by(cluster, group, module) %>%
  summarise(mean_score = mean(score), .groups = "drop") %>%
  mutate(module  = factor(module,  levels = mods_use),
         cluster = factor(cluster, levels = 1:6),
         group   = factor(group,
                          levels = c("Aged","R0","R1-3","R4-6","R7-9","Young")))

make_line_panel <- function(mod_name) {
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
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          panel.grid.minor = element_blank())
}

row1 <- (make_line_panel("M1") | make_line_panel("M3") | plot_spacer())
row2 <- (make_line_panel("M5") | make_line_panel("M6") | make_line_panel("M9"))
row3 <- (make_line_panel("M2") | make_line_panel("M11") | plot_spacer())

p_module_dyn <- (row1 / row2 / row3) + plot_layout(guides = "collect")
print(p_module_dyn)


# ---- 6. 저장 --------------------------------------------------
saveRDS(seu2_ds, file.path(base_dir, "seu2_ds.rds"))
qs::qsave(seu2_ds, file.path(base_dir, "seu2_ds.qs"), nthreads = 8)

cat("=== DONE: 02_downsample_modulescores ===\n")
cat("Saved:\n",
    file.path(base_dir, "seu2_ds.rds"), "\n",
    file.path(base_dir, "seu2_ds.qs"),  "\n")
