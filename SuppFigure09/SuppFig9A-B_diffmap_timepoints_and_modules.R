## ============================ Section 1 =====================================
## R0-R9 diffusion map.
## Not duplicated here: the same code is the section "R0-R9 diffusion map" of
## ../Figure3/Fig3H_01_R0R9_diffusion_map_input.R (Fig 3H upstream). The lines below run only
## that section, in this R session (the h5ad export that follows it
## in Fig3H_01 is not needed here). It creates seu_r ("diffmap" reduction DC1-DC5, dpt), dm,
## dc_mat, mods and mod_colors, and loads Seurat, harmony, destiny, dplyr, tidyr, tibble,
## ggplot2, patchwork, ggrepel, viridisLite, RColorBrewer, scales, colorspace and mgcv.
FIG3H_01 <- "../Figure3/Fig3H_01_R0R9_diffusion_map_input.R"
fig3h_src  <- readLines(FIG3H_01, encoding = "UTF-8")
fig3h_from <- grep("^# ---- R0-R9 diffusion map", fig3h_src)
fig3h_to   <- grep("^# ---- export seu_r -> h5ad", fig3h_src) - 1
stopifnot(length(fig3h_from) == 1, length(fig3h_to) == 1, fig3h_to > fig3h_from)
eval(parse(text = fig3h_src[fig3h_from:fig3h_to], encoding = "UTF-8"))
# ---- 7. Visualization df --------------------------------------------
pal_R <- viridis(10); names(pal_R) <- paste0("R", 0:9)
df <- as.data.frame(dc_mat)
df$sample <- seu_r$sample
df$batch  <- seu_r$batch
df$dpt    <- seu_r$dpt
df[mods]  <- seu_r@meta.data[, mods]

## ============================ Section 2 (9A) ================================
# (b) per-sample facet
p_facet <- ggplot(df, aes(DC1, DC2)) +
  geom_point(data = df %>% select(-sample),
             color = "grey88", size = 0.25) +
  geom_point(aes(color = sample), size = 0.45) +
  scale_color_manual(values = pal_R, guide = "none") +
  facet_wrap(~ sample, ncol = 5) +
  theme_classic(base_size = 10) + coord_fixed() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(title = "Per-sample distribution")
ggsave(file.path(out_dir, "02_DC1DC2_facet.pdf"),
       p_facet, width = 12, height = 5)

## ============================ Section 3 (9B) ================================
plot_mod <- function(m) {
  v   <- df[[m]]
  lim <- quantile(v, c(0.05, 0.95), na.rm = TRUE)
  hi  <- mod_colors[m]
  pal <- c("grey92", lighten(hi, 0.3), hi, darken(hi, 0.3))
  ggplot(df[order(v), ], aes(DC1, DC2, color = .data[[m]])) +
    geom_point(size = 0.4, alpha = 0.95) +
    scale_color_gradientn(colors = pal, limits = lim, oob = squish,
                          guide = "none") +
    theme_classic(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 11, face = "bold")) +
    labs(title = m, x = NULL, y = NULL)
}
p_mod <- wrap_plots(lapply(mods, plot_mod), ncol = 4) +
  plot_annotation(title = "Modules on DC1-DC2")
p_mod

## ============================ Section 4 (9B variant) ========================
## Same overlay on the CellOracle input object (DC1-DC2 of its "diffmap" reduction).
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(RColorBrewer); library(viridisLite)
})
options(Seurat.object.assay.version = "v3")
base_dir <- "EGM_multiflex/oracle"
seu <- readRDS(file.path(base_dir, "oracle_input_hvg.rds"))
## Leiden res 0.1 on the diffusion map -> seu$dm_cluster.
## dm_cluster is only carried in df_dm below; it does not change the plot.
# DM 위에서 Leiden res 0.1 재계산
seu <- FindNeighbors(seu, reduction = "diffmap", dims = 1:5,
                     graph.name = c("dm_nn", "dm_snn"), verbose = FALSE)
seu <- FindClusters(seu, graph.name = "dm_snn",
                    resolution = 0.1, algorithm = 4,
                    leiden_method = "igraph", random.seed = 42,
                    cluster.name = "dm_leiden_0.1", verbose = FALSE)
seu$dm_cluster <- factor(seu$dm_leiden_0.1)
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(RColorBrewer)
})
mods <- c("M1","M3","M5","M6","M9","M2","M11")
mod_colors <- c(M1="#C0392B", M3="#E67E22", M5="#2A9D8F",
                M6="#52B788", M9="#377EB8", M2="#7B5BA6", M11="#C39BD3")
# ---- (c) DM 위에 module overlay (cluster 외곽 보기) ------------
library(scales); library(colorspace)
dc <- Embeddings(seu, "diffmap")
df_dm <- data.frame(DC1 = dc[,1], DC2 = dc[,2],
                    cluster = seu$dm_cluster)
df_dm[mods] <- seu@meta.data[, mods]
plot_mod <- function(m) {
  v <- df_dm[[m]]
  lim <- quantile(v, c(0.05, 0.95), na.rm = TRUE)
  hi <- mod_colors[m]
  pal <- c("grey92", lighten(hi, 0.3), hi, darken(hi, 0.3))
  ggplot(df_dm[order(v), ], aes(DC1, DC2, color = .data[[m]])) +
    geom_point(size = 0.35, alpha = 0.95) +
    scale_color_gradientn(colors = pal, limits = lim, oob = squish,
                          guide = "none") +
    theme_classic(base_size = 9) + coord_fixed() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(face = "bold", size = 10)) +
    labs(title = m, x = NULL, y = NULL)
}
p_mods_dm <- wrap_plots(lapply(mods, plot_mod), ncol = 4) +
  plot_annotation(title = "Module score on DM")
print(p_mods_dm)
