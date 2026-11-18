# Legends of the module heatmap (Z-score, -log(adjP)) on their own, titles not bold, Arial, 1200 dpi.
# Same colour functions and breaks as Fig3D_compare_module_heatmaps.R.
# Output: heatmap_legend_plain.pdf / _1200dpi.png
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(RColorBrewer); library(grid) })
ROOT <- "4FP_plasma/EGM_hdWGCNA_code"
DIR <- file.path(ROOT, "03_figures_R/output")
FONT <- "Arial"
col_z     <- colorRamp2(c(-2.4, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.4), rev(brewer.pal(11, "RdBu")))
col_aging <- colorRamp2(c(0, 1.3, 3, 6, 10), c("white", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"))
tgp <- gpar(fontsize = 7.5, fontface = "plain", fontfamily = FONT)
lgp <- gpar(fontsize = 6.5, fontfamily = FONT)
lg <- packLegend(
  Legend(col_fun = col_z, title = "Z-score", at = c(-2, -1, 0, 1, 2), labels = c("-2", "-1", "0", "1", "2"),
         direction = "horizontal", legend_width = unit(2, "cm"), title_gp = tgp, labels_gp = lgp),
  Legend(col_fun = col_aging, title = "-log(adjP)", at = c(0, 1.3, 3, 6, 10), labels = c("0", "1.3", "3", "6", ">=10"),
         direction = "horizontal", legend_width = unit(2, "cm"), title_gp = tgp, labels_gp = lgp),
  direction = "horizontal", gap = unit(5, "mm"))
pad <- unit(1, "mm")
w <- convertWidth(ComplexHeatmap:::width(lg) + 2 * pad, "in", valueOnly = TRUE) + 0.12   # room for ">=10" past the bar end
h <- convertHeight(ComplexHeatmap:::height(lg) + 2 * pad, "in", valueOnly = TRUE)
f <- file.path(DIR, "heatmap_legend_plain")
quartz(type = "pdf", file = paste0(f, ".pdf"), width = w, height = h, family = FONT); draw(lg); invisible(dev.off())
ragg::agg_png(paste0(f, "_1200dpi.png"), width = w, height = h, units = "in", res = 1200, background = "white")
draw(lg); invisible(dev.off())
cat(sprintf("%.2f x %.2f in\n", w, h))
