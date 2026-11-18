# ------------------------------------------------------------------
# Liver DEL heatmap (lipidomics, HGPS EGM vs control, n = 2 vs 2)
#
# Data and settings:
#   input  n_liver_final2.csv (67 lipids x 4 samples)
#   z-score per lipid (pheatmap scale = "row"), lipids clustered on the z-scores
#   (euclidean, complete), blue-white-red over +-max|z|, cell border grey60,
#   group colours grey75 (Control) / cornflowerblue (EGM)
# Drawn horizontally (3.07 x 1.69 in) with the lipid order of the vertical heatmap rotated clockwise
# (its bottom lipid on the left, so the dendrogram is mirrored), lipid names reading bottom-to-top,
# the group bar on the LEFT, Arial.
# Outputs:
#   Liver_DEL_heatmap.pdf / _1200dpi.png         heatmap + dendrogram + lipid names + group bar
#   Liver_DEL_heatmap_legend.pdf / _1200dpi.png  Group + Z-score legends
# ------------------------------------------------------------------
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid); library(pheatmap) })

IN  <- "LCMS/n_liver_final2.csv"
OUT <- "4FP_plasma/EGM_lipid_heatmap"
FONT <- "Arial"

# ---- sizes ------------------------------------------------------------
FIG_W <- 3.07; FIG_H <- 1.69          # in
LIPID_PT <- 3.2                       # lipid names (pt)
DEND_PT  <- 22                        # dendrogram height (bp)
LWD_CELL <- 0.3 / 0.75                # 0.3 pt cell borders (lwd units are 1/96 in)
LWD_DEND <- 0.45 / 0.75

# ---- data + clustering (as pheatmap does it) ------------------------------
df  <- read.csv(IN, check.names = FALSE)
mat <- as.matrix(df[, c("P_Con_Liver_1", "P_Con_Liver_2", "P_EGM_Liver_1", "P_EGM_Liver_2")])
rownames(mat) <- df$name
z  <- t(scale(t(mat)))                                    # z-score per lipid
hc <- hclust(dist(z), method = "complete")
ph <- pheatmap(mat, scale = "row", cluster_cols = FALSE, clustering_distance_rows = "euclidean",
               clustering_method = "complete", silent = TRUE)
stopifnot(identical(hc$order, ph$tree_row$order))        # same order as pheatmap
m_abs <- max(abs(z))
col_fun <- colorRamp2(c(-m_abs, 0, m_abs), c("blue", "white", "red"), space = "RGB")
grp_col <- c(Control = "grey75", EGM = "cornflowerblue")

# horizontal: samples as rows, EGM on top and Control below; lipids as columns
row_ord <- c("P_EGM_Liver_1", "P_EGM_Liver_2", "P_Con_Liver_1", "P_Con_Liver_2")
m_h <- t(z)[row_ord, ]
grp <- ifelse(grepl("EGM", row_ord), "EGM", "Control")
dend <- rev(as.dendrogram(hc))                            # mirrored = vertical figure rotated clockwise
cat("left -> right:", head(labels(dend), 3), "...", tail(labels(dend), 3), "\n")

# PDF device: quartz (CoreText) keeps exact Arial glyph advances; cairo_pdf rounds them at these
# small sizes (letters of a 3.2 pt name end up unevenly spaced, e.g. "L PC").
open_pdf <- function(f, w, h) quartz(type = "pdf", file = f, width = w, height = h, family = FONT)
# open the device first: ComplexHeatmap measures the lipid names when the object is built
f_ht <- file.path(OUT, "Liver_DEL_heatmap.pdf")
open_pdf(f_ht, FIG_W, FIG_H)
ht <- Heatmap(m_h, name = "Z-score", col = col_fun,
              cluster_rows = FALSE, cluster_columns = dend, column_dend_reorder = FALSE,
              column_dend_height = unit(DEND_PT, "bigpts"), column_dend_gp = gpar(lwd = LWD_DEND),
              show_row_names = FALSE, column_names_side = "bottom", column_names_rot = 90,
              column_names_gp = gpar(fontsize = LIPID_PT, fontfamily = FONT),
              rect_gp = gpar(col = "grey60", lwd = LWD_CELL),
              left_annotation = rowAnnotation(Group = grp, col = list(Group = grp_col),
                                               show_annotation_name = FALSE, show_legend = FALSE,
                                               simple_anno_size = unit(5, "bigpts"),
                                               gp = gpar(col = "grey60", lwd = LWD_CELL)),
              show_heatmap_legend = FALSE,
              heatmap_width = unit(1, "npc"), heatmap_height = unit(1, "npc"))

ht_opt$HEATMAP_LEGEND_PADDING <- unit(0, "mm"); ht_opt$ANNOTATION_LEGEND_PADDING <- unit(0, "mm")
draw(ht, padding = unit(c(1, 0.5, 0.5, 0.5), "bigpts"), gap = unit(1.5, "bigpts"))   # bottom, left, top, right
dev.off()

# ---- legends: Group (title 5 pt, labels 4 pt) and Z-score bar (labels 4 pt, title 6 pt below) ----
LEG_W <- 0.42; LEG_H <- 1.30          # in
f_lg <- file.path(OUT, "Liver_DEL_heatmap_legend.pdf")
open_pdf(f_lg, LEG_W, LEG_H)
grid.newpage()
pushViewport(viewport(xscale = c(0, LEG_W * 72), yscale = c(LEG_H * 72, 0)))   # bp, y downwards
txt <- function(lab, x, y, size, just = c("left", "center"))
  grid.text(lab, unit(x, "native"), unit(y, "native"), just = just, gp = gpar(fontsize = size, fontfamily = FONT))
txt("Group", 1, 3.5, 5)
sq <- 5.7
for (i in 1:2) {
  y0 <- 8 + (i - 1) * (sq + 1.6)
  grid.rect(unit(1, "native"), unit(y0, "native"), unit(sq, "native"), unit(-sq, "native"),
            just = c("left", "top"), gp = gpar(fill = grp_col[i], col = NA))
  txt(names(grp_col)[i], 1 + sq + 2.2, y0 + sq / 2, 4)
}
bar_x <- 3; bar_w <- 3.9; bar_top <- 28.5; bar_h <- 57          # bar from +m_abs (top) to -m_abs (bottom)
ys <- seq(m_abs, -m_abs, length.out = 201)
for (k in seq_len(length(ys) - 1)) {
  grid.rect(unit(bar_x, "native"), unit(bar_top + (k - 1) * bar_h / 200, "native"), unit(bar_w, "native"),
            unit(-(bar_h / 200 + 0.05), "native"), just = c("left", "top"),
            gp = gpar(fill = col_fun(mean(ys[k:(k + 1)])), col = NA))
}
for (v in c(1, 0.5, 0, -0.5, -1)) {
  yv <- bar_top + (m_abs - v) / (2 * m_abs) * bar_h
  grid.segments(unit(bar_x + bar_w - 0.9, "native"), unit(yv, "native"), unit(bar_x + bar_w, "native"), unit(yv, "native"),
                gp = gpar(col = "white", lwd = 0.25 / 0.75))
  txt(sub("-", "−", format(v)), bar_x + bar_w + 2.2, yv, 4)
}
txt("Z-score", 5.0, bar_top + bar_h + 6, 6)
popViewport()
dev.off()

for (f in c(f_ht, f_lg))
  system2("pdftoppm", c("-png", "-r", "1200", "-singlefile", shQuote(f), shQuote(sub("\\.pdf$", "_1200dpi", f))))
cat("max |z| =", round(m_abs, 3), "\n")
