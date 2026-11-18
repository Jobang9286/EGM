# =============================================================================
# Module heatmap with all 13 modules (M1-M13, no merging):
# z-score block (A, Y, R0-R9) + 12 aging hallmarks, 5-mm cells, no identity block, legends left out.
# Values come from Fig3D_compare_module_heatmaps.R ("unmerged"), same cells and settings as May.
# Output: heatmap_13modules_labeled.pdf / _1200dpi.png (with labels), heatmap_13modules_tiles.pdf / _1200dpi.png
# =============================================================================
ROOT <- "4FP_plasma/EGM_hdWGCNA_code"
CODE_DIR <- "."   # this repository's Figure3/ folder (holds Fig3D_compare_module_heatmaps.R)
DIR <- file.path(ROOT, "03_figures_R/output")
ex  <- parse(file.path(CODE_DIR, "Fig3D_compare_module_heatmaps.R"))
skip <- vapply(ex, function(e) grepl("^(draw_one\\(|compare <-|for \\(o in)", deparse(e)[1]), logical(1))
for (e in ex[!skip]) eval(e)                    # builds res$unmerged, col_z, col_aging, CELL
# the "SASP" column is Reactome + WP SASP plus Hallmark TNFA_SIGNALING_VIA_NFKB and INFLAMMATORY_RESPONSE: label it as such
for (k in names(res)) colnames(res[[k]]$mat_a)[colnames(res[[k]]$mat_a) == "SASP"] <-
  colnames(res[[k]]$sig_a)[colnames(res[[k]]$sig_a) == "SASP"] <- "SASP/Inflammation"
FONT <- "Arial"

paper_ht <- function(r, labels) {
  z <- r$z; colnames(z) <- c("A", "Y", paste0("R", 0:9)); sz <- r$sig_z
  txt <- function(s) function(j, i, x, y, w, h, fill) grid.text(s[i, j], x, y, gp = gpar(fontsize = 5, fontfamily = FONT))
  Heatmap(z, name = "z", col = col_z, cluster_rows = FALSE, cluster_columns = FALSE,
          show_row_names = labels, row_names_side = "left", row_names_gp = gpar(fontsize = 7, fontfamily = FONT),
          show_column_names = labels, column_names_gp = gpar(fontsize = 7, fontfamily = FONT),
          column_names_rot = 0, column_names_centered = TRUE,
          column_split = factor(c("A", "Y", rep("R", 10)), levels = c("A", "Y", "R")), column_title = NULL,
          column_gap = unit(1, "mm"), width = unit(CELL * 12, "mm"), height = unit(CELL * nrow(z), "mm"),
          rect_gp = gpar(col = "white", lwd = 0.6), cell_fun = txt(sz), show_heatmap_legend = FALSE) +
  Heatmap(r$mat_a, name = "h", col = col_aging, cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
          show_column_names = labels, column_names_gp = gpar(fontsize = 7, fontfamily = FONT), column_names_rot = 45,
          width = unit(CELL * ncol(r$mat_a), "mm"), height = unit(CELL * nrow(z), "mm"),
          rect_gp = gpar(col = "white", lwd = 0.6), cell_fun = txt(r$sig_a), show_heatmap_legend = FALSE)
}

for (v in list(list(TRUE, "heatmap_13modules_labeled"), list(FALSE, "heatmap_13modules_tiles"))) {
  ht <- paper_ht(res$unmerged, v[[1]])
  pdf(NULL); htd <- draw(ht, padding = unit(c(1, 1, 1, 1), "mm"))    # measure the drawn size
  w <- convertWidth(ComplexHeatmap:::width(htd), "in", valueOnly = TRUE)
  h <- convertHeight(ComplexHeatmap:::height(htd), "in", valueOnly = TRUE); invisible(dev.off())
  f <- file.path(DIR, v[[2]])
  quartz(type = "pdf", file = paste0(f, ".pdf"), width = w, height = h, family = FONT)
  draw(ht, padding = unit(c(1, 1, 1, 1), "mm")); invisible(dev.off())
  ragg::agg_png(paste0(f, "_1200dpi.png"), width = w, height = h, units = "in", res = 1200, background = "white")
  draw(ht, padding = unit(c(1, 1, 1, 1), "mm")); invisible(dev.off())
  cat(sprintf("%s: %.2f x %.2f in\n", v[[2]], w, h))
}

# labeled PNG: cut the empty strip left under the rotated hallmark labels, keeping 50 px below the last inked row
# (the PDF is left as drawn)
f <- file.path(DIR, "heatmap_13modules_labeled_1200dpi.png"); im <- png::readPNG(f)
g <- (0.299 * im[, , 1] + 0.587 * im[, , 2] + 0.114 * im[, , 3]) * 255
last <- max(which(apply(g < 245, 1, any)))
png::writePNG(im[seq_len(min(last + 50, nrow(g))), , , drop = FALSE], f, dpi = 1200)
cat(sprintf("labeled PNG cropped to %d x %d px\n", dim(im)[2], min(last + 50, nrow(g))))
