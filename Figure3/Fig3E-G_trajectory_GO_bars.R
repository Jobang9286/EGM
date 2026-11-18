# =============================================================================
# Trajectory plots (3 panels) + module-specific GO BP bar charts
#   Group 1: Aging feature, decreasing      (M1, M3)        A1=100, Y1=0
#   Group 2: Young feature, early           (M5, M6, M9)    Y1=100, A1=0
#   Group 3: Young feature, mid-stage       (M2, M11)       Y1=100, A1=0
# Smoothing:
#   - line   : loess(family="symmetric", degree=2, span=0.75)   ("LOWESS-like")
#   - ribbon : loess-smoothed (mean +/- sample SEM)             (no endpoint blow-up)
# Right side of each panel: top-2 GO BP terms per module, bar = -log10(FDR)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})

# ---- paths --------------------------------------------------------------
DIR_IN  <- "egm_final/for_rstudio"
OUT_PDF <- "egm_final/trajectory_with_go.pdf"
OUT_PNG <- "egm_final/trajectory_with_go.png"

# ---- input + normalization ---------------------------------------------
df <- read.csv(file.path(DIR_IN, "trajectory_long.csv"))

sample_lvl <- c("A1","R0","R1","R2","R3","R4","R5","R6","R7","R8","R9","Y1")
df$sample  <- factor(df$sample, levels = sample_lvl)
df$x_num   <- as.integer(df$sample)

aging_mods <- c("M1","M3")
young_mods <- c("M5","M6","M9","M2","M11")

df_norm <- df %>%
  filter(module %in% c(aging_mods, young_mods)) %>%
  group_by(module) %>%
  mutate(
    x_A1 = mean_hME[sample == "A1"],
    x_Y1 = mean_hME[sample == "Y1"],
    is_age = first(module) %in% aging_mods,
    norm = ifelse(is_age,
                  100 * (mean_hME - x_Y1) / (x_A1 - x_Y1),
                  100 * (mean_hME - x_A1) / (x_Y1 - x_A1)),
    norm_sem = 100 * sem_hME / abs(x_A1 - x_Y1)
  ) %>% ungroup()

# ---- LOWESS-style smoothing on fine grid --------------------------------
xgrid <- seq(1, 12, length.out = 400)

smooth_line_sem <- function(d, span = 0.75) {
  d <- d %>% arrange(x_num)
  f_mid <- loess(norm ~ x_num,                    data = d, span = span,
                 family = "symmetric", degree = 2)
  f_lo  <- loess(I(norm - norm_sem) ~ x_num,      data = d, span = span,
                 family = "symmetric", degree = 2)
  f_hi  <- loess(I(norm + norm_sem) ~ x_num,      data = d, span = span,
                 family = "symmetric", degree = 2)
  tibble(
    x_num = xgrid,
    norm  = predict(f_mid, data.frame(x_num = xgrid)),
    lo    = predict(f_lo , data.frame(x_num = xgrid)),
    hi    = predict(f_hi , data.frame(x_num = xgrid))
  )
}
df_smooth <- df_norm %>% group_by(module) %>%
  group_modify(~ smooth_line_sem(.x)) %>% ungroup()

# ---- curated top-2 GO BP terms per module + FDR ------------------------
# Values from clusterProfiler-style Fisher enrichment vs module gene lists
# (background = all WGCNA universe genes; MSigDB C5 GO:BP, mouse)
go_df <- tribble(
  ~module, ~term,                              ~padj,
  "M1",  "ECM / collagen organization",        1.89e-10,
  "M1",  "Inflammatory response",              1.40e-9,
  "M3",  "Interferon-mediated signaling",      1.69e-7,
  "M3",  "Cytokine production",                7.41e-8,
  "M5",  "Sterol biosynthesis",                1.15e-12,
  "M5",  "Lipid biosynthesis",                 3.28e-6,
  "M6",  "Vesicle-mediated transport",         8.84e-3,
  "M6",  "Autophagy",                          1.02e-2,
  "M9",  "Actin filament organization",        1.06e-7,
  "M9",  "Cytoskeleton organization",          4.20e-5,
  "M2",  "RNA processing",                     6.18e-36,
  "M2",  "Ribosome biogenesis",                3.33e-25,
  "M11", "Cell cycle / mitosis",               1.26e-32,
  "M11", "Chromosome segregation",             7.24e-31
)
go_df$nlog <- pmin(-log10(go_df$padj), 30)   # cap for visualization

# ---- palettes / sizes ---------------------------------------------------
pal1 <- c(M1 = "#c0392b", M3 = "#e67e22")
pal2 <- c(M5 = "#16a085", M6 = "#27ae60", M9 = "#2980b9")
pal3 <- c(M2 = "#8e44ad", M11 = "#d4a5dc")

S_CANVAS  <- 0.448   # canvas scale relative to base 8 x 9 inch
BASE_FONT <- 8

# ---- panel constructors -------------------------------------------------
make_panel <- function(mods, pal, ttl, ylab, yr) {
  d <- df_smooth %>% filter(module %in% mods)
  d$module <- factor(d$module, levels = mods)
  ggplot(d, aes(x_num, norm, color = module, fill = module, group = module)) +
    geom_hline(yintercept = c(0, 100), linetype = "dashed",
               color = "grey70", linewidth = 0.25) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(breaks = 1:12, labels = sample_lvl,
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = c(0, 100)) +
    coord_cartesian(ylim = yr) +
    labs(title = ttl, x = NULL, y = ylab) +
    theme_classic(base_size = BASE_FONT) +
    theme(
      plot.title       = element_text(size = BASE_FONT, face = "bold"),
      legend.position  = "right", legend.title = element_blank(),
      legend.key.size  = unit(BASE_FONT, "pt")
    )
}

make_go_bar <- function(mods, pal) {
  d <- go_df %>% filter(module %in% mods)
  d$module <- factor(d$module, levels = mods)
  d$term   <- factor(d$term,   levels = rev(d$term))
  ggplot(d, aes(nlog, term, fill = module)) +
    geom_vline(xintercept = -log10(0.05), color = "black", linewidth = 0.3) +
    geom_col(width = 0.65) +
    scale_fill_manual(values = pal) +
    facet_grid(module ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = expression(-log[10]*"[FDR]"), y = NULL) +
    theme_classic(base_size = BASE_FONT) +
    theme(
      axis.text.y      = element_text(size = BASE_FONT - 1, hjust = 1),
      axis.text.x      = element_text(size = BASE_FONT - 1),
      axis.title.x     = element_text(size = BASE_FONT),
      strip.placement  = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1,
                                       size = BASE_FONT, face = "bold"),
      panel.spacing.y  = unit(2, "pt"),
      legend.position  = "none"
    )
}

# ---- assemble panels ----------------------------------------------------
p1 <- make_panel(aging_mods,        pal1,
                 "Aging feature, decreasing (A1=100, Y1=0)",
                 "Aging-similarity (%)",   c(-25, 120))
p2 <- make_panel(c("M5","M6","M9"), pal2,
                 "Young feature, early acquisition (Y1=100, A1=0)",
                 "Young-similarity (%)",   c(-25, 140))
p3 <- make_panel(c("M2","M11"),     pal3,
                 "Young feature, mid-stage acquisition (Y1=100, A1=0)",
                 "Young-similarity (%)",   c(-25, 120))

g1 <- make_go_bar(aging_mods,        pal1)
g2 <- make_go_bar(c("M5","M6","M9"), pal2)
g3 <- make_go_bar(c("M2","M11"),     pal3)

g <- (p1 + g1 + plot_layout(widths = c(2, 1.2))) /
     (p2 + g2 + plot_layout(widths = c(2, 1.2))) /
     (p3 + g3 + plot_layout(widths = c(2, 1.2)))

# ---- save --------------------------------------------------------------
W <- (8 * S_CANVAS) * 1.8
H <-  9 * S_CANVAS
ggsave(OUT_PDF, g, width = W, height = H)
ggsave(OUT_PNG, g, width = W, height = H, dpi = 300)

cat("saved:\n  ", OUT_PDF, "\n  ", OUT_PNG, "\n", sep = "")
