# ------------------------------------------------------------------
# Added boxplot rows for the deconvolution aging figure (HGPS_EGM vs HGPS_control, n = 3 vs 3)
#   Kidney : distal convoluted tubule (DCT)  -> Klf6, Atf4, Ccnb1, Igf1
#   Liver  : midzonal hepatocytes            -> Lcn2, Cp, Ppara, Cyp8b1
#
# Panel style is the one of the existing rows (kidney PT / Endo and liver zone boxplots):
#   - kidney: y ticks from breaks_half() (<= 3 ticks, 0.5 steps); p shown in the figure = t-test on log2(Z + 1)
#   - liver : y ticks from y_scale_half() (axis widened to the ticks); p shown = DESeq2 padj
#   - points jittered horizontally only (height = 0): vertical jitter, with the liver axis limits set to the
#     data range, pushes edge points out of the scale
# Layout: plots only, no text and no white margin. Panel frames are placed with the geometry of the
# existing rows, so the new row lines up when it is scaled to the same width. Panel height = 156.2 pt as in
# the existing R figures, so frame, tick and point sizes are the same. The page is cut at the frame edge
# (top/right) and at the tick ends (left/bottom).
# Outputs:
#   <tissue>_<celltype>_boxplot_notext.pdf / _notext_1200dpi.png  (PNG rendered from the PDF with pdftoppm)
#   <tissue>_<celltype>_boxplot_ylabels.pdf / _ylabels_1200dpi.png (same, plus y tick values in Arial, sized to
#        3.5 pt once the row is scaled to the existing rows; Arial embedded)
#   added_rows_stats.csv  (means, LFC, padj, t-test p, y tick values bottom -> top)
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(grid); library(gtable); library(InstaPrism)
})

ROOT <- "EGM_project_for_deconvolution/integrated_analysis_for_lizard_factors"
OUT  <- "4FP_plasma/EGM_deconv_added_rows"

# ---- shared style (identical to the existing figures) -------------
cond_colors <- c("HGPS_control" = "#D9D9D9", "HGPS_EGM" = "#7FA6FF")
point_size  <- 1.6 * 0.7
TICK        <- 3.25                                   # pt, theme_bw(base_size = 13) tick length
BORDER_LW   <- 1.3                                    # panel frame line
base_theme <- function(ylab_size = NULL) theme_bw(base_size = 13) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = BORDER_LW),
        axis.text = element_blank(), axis.title = element_blank(),            # no text at all
        plot.title = element_blank(), plot.caption = element_blank(),
        axis.ticks.length = unit(TICK, "pt"), plot.margin = margin(0, 0, 0, 0)) +
  (if (is.null(ylab_size)) theme() else theme(axis.text.y = element_text(size = ylab_size, colour = "black", family = "Arial")))
gene_df <- function(Z, gene, meta) {
  e <- Z[gene, ]
  data.frame(Sample = names(e), Expression = e, Group = meta$Group[match(names(e), meta$Sample)]) %>%
    filter(Group %in% c("HGPS_control", "HGPS_EGM")) %>% mutate(Group = factor(Group, levels = c("HGPS_control", "HGPS_EGM")))
}
tt_one <- function(Z, g, meta) {                       # Z = raw counts (genes x samples)
  ctl <- meta$Sample[meta$Group == "HGPS_control"]; egm <- meta$Sample[meta$Group == "HGPS_EGM"]
  e <- Z[g, ]; l <- log2(e + 1)
  data.frame(gene = g, Mean_Control = mean(e[ctl]), Mean_EGM = mean(e[egm]),
             P_ttest_raw = t.test(e[egm], e[ctl])$p.value, P_ttest_log2 = t.test(l[egm], l[ctl])$p.value)
}

# ---- y-axis rules ---------------------------------------------------
breaks_half <- function(lims) {                       # kidney rows
  steps <- c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 15, 20)
  for (s in steps) { b <- seq(ceiling(lims[1] / s) * s, floor(lims[2] / s) * s, by = s); if (length(b) <= 3) return(b) }
  b
}
SEARCH <- 0.10
y_scale_half <- function(v) {                         # liver rows
  r <- range(v); span <- diff(r); ex <- r + c(-1, 1) * span * SEARCH
  b <- NULL
  for (s in c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 15, 20)) {
    bb <- seq(ceiling(ex[1] / s) * s, floor(ex[2] / s) * s, by = s)
    if (length(bb) >= 2 && length(bb) <= 3) { b <- bb; break }
    if (length(bb) < 2) break
  }
  if (is.null(b)) b <- seq(floor(r[1] / 0.5) * 0.5, ceiling(r[2] / 0.5) * 0.5, by = 0.5)
  lim <- c(min(r[1], min(b)), max(r[2], max(b)))
  scale_y_continuous(breaks = b, limits = lim, expand = expansion(mult = 0.06))
}
box_panel <- function(Z, gene, meta, y_rule, ylab_size = NULL) {
  d <- gene_df(Z, gene, meta)
  ggplot(d, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.6, colour = "black", linewidth = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, size = point_size, alpha = 0.8) +
    scale_fill_manual(values = cond_colors) +
    (if (y_rule == "kidney") scale_y_continuous(breaks = breaks_half) else y_scale_half(d$Expression)) +
    base_theme(ylab_size)
}
y_ticks <- function(Z_log, g, meta, y_rule) {
  bt <- ggplot_build(box_panel(Z_log, g, meta, y_rule))$layout$panel_params[[1]]$y$get_labels()
  paste(bt[!is.na(bt)], collapse = ", ")
}

# ---- row geometry ----------------------------------------------------
# frame-centre coordinates (px) measured in the existing row images, divided by the panel height
GEOM <- list(liver  = list(w = 238.4  / 306.75, left = c(0, 319.5, 622, 941.5) / 306.75),
             kidney = list(w = 225.25 / 309.5,  left = c(0, 281.5, 571, 860.5) / 309.5))
PANEL_H <- 156.2                                     # bp (1/72 in), measured in the original PDFs
BHALF   <- BORDER_LW * .pt * 72 / 96 / 2             # half of the panel.border line (bp)
TICK_BP <- TICK * 72 / 72.27                         # tick length in bp (grid "pt" is 1/72.27 in)
# span from the first to the last frame line of the existing rows (in): liver 1180 of 1210 px at
# 2.07 in, kidney 1086 of 1126 px at 2.09 in -> scale of the new row once it is placed next to them
PPT_SPAN_IN <- c(liver = 1180 / 1210 * 2.07, kidney = 1086 / 1126 * 2.09)
PPT_TICK_PT <- 3.5                                   # tick label size (pt) at the final scale

draw_row <- function(Z_log, genes, meta, y_rule, file, ylabels = FALSE) {
  g <- GEOM[[y_rule]]; pw <- g$w * PANEL_H; lefts <- g$left * PANEL_H
  ppt_scale <- PPT_SPAN_IN[[y_rule]] * 72 / (max(lefts) + pw)          # bp at the final scale per bp here
  ysz <- if (ylabels) PPT_TICK_PT / ppt_scale else NULL
  open_dev <- function(f, w = 7, h = 7) if (ylabels) cairo_pdf(f, width = w, height = h, family = "Arial") else pdf(f, width = w, height = h)
  make_gt <- function(i) {
    set.seed(42 + i)                                   # reproducible jitter
    gt <- ggplotGrob(box_panel(Z_log, genes[i], meta, y_rule, ysz))
    pc <- gt$layout[gt$layout$name == "panel", ]
    # make every column/row absolute (the axis column carries 'null' units that would otherwise share
    # space with the panel), then fix the panel itself to pw x PANEL_H
    aw <- convertWidth(gt$widths, "bigpts", valueOnly = TRUE);   aw[pc$l] <- pw;      gt$widths  <- unit(aw, "bigpts")
    ah <- convertHeight(gt$heights, "bigpts", valueOnly = TRUE); ah[pc$t] <- PANEL_H; gt$heights <- unit(ah, "bigpts")
    list(gt = gt, lw = sum(aw[seq_len(pc$l - 1)]), bh = sum(ah[-seq_len(pc$b)]), W = sum(aw), H = sum(ah))
  }
  tmp <- tempfile(fileext = ".pdf"); open_dev(tmp); lw1 <- make_gt(1)$lw; dev.off(); unlink(tmp)   # left axis width of panel 1
  # tick labels at the top/bottom tick stick out of the frame by up to half a digit height
  over <- c(t = 0, b = 0)
  if (ylabels) for (i in seq_along(genes)) {
    yp <- ggplot_build(box_panel(Z_log, genes[i], meta, y_rule))$layout$panel_params[[1]]$y
    br <- yp$breaks[!is.na(yp$breaks)]; yb <- (br - yp$continuous_range[1]) / diff(yp$continuous_range) * PANEL_H
    half <- 0.38 * ysz                                 # Arial digit height is ~0.72 em
    over <- pmax(over, c(t = max(yb) + half - PANEL_H, b = half - min(yb)))
  }
  pad <- c(l = max(TICK_BP, lw1), r = BHALF, t = max(BHALF, over[["t"]]), b = max(TICK_BP, over[["b"]]))  # all in bp
  W <- pad[["l"]] + max(lefts) + pw + pad[["r"]]; H <- pad[["b"]] + PANEL_H + pad[["t"]]
  open_dev(file, W / 72, H / 72)
  grid.newpage()
  for (i in seq_along(genes)) {
    m <- make_gt(i)
    pushViewport(viewport(x = unit(pad[["l"]] + lefts[i] - m$lw, "bigpts"), y = unit(pad[["b"]] - m$bh, "bigpts"),
                          width = unit(m$W, "bigpts"), height = unit(m$H, "bigpts"), just = c("left", "bottom")))
    grid.draw(m$gt); popViewport()
  }
  dev.off()
  system2("pdftoppm", c("-png", "-r", "1200", "-singlefile", shQuote(file), shQuote(sub("\\.pdf$", "_1200dpi", file))))
  cat(sprintf("%s: %.1f x %.1f bp -> width in the PPT %.3f in%s\n", basename(file), W, H, W * ppt_scale / 72,
              if (ylabels) sprintf(" (tick labels %.1f pt here = %g pt in the PPT)", ysz, PPT_TICK_PT) else ""))
}

# ---- kidney: DCT ---------------------------------------------------
load(file.path(ROOT, "kidney/deg_results_2/kidney_deconvolution_main_objects.RData"))   # Z_array, Z_PT_full, sample_metadata
Z_DCT <- t(Z_array[, , "DCT"])
stopifnot(isTRUE(all.equal(t(Z_array[, , "PT"]), Z_PT_full, check.attributes = FALSE)))  # same layout as the PT row source
kid_meta  <- sample_metadata
kid_genes <- c("Klf6", "Atf4", "Ccnb1", "Igf1")
kid_deg   <- read.csv(file.path(ROOT, "kidney/deg_results/cell_type_specific/DCT_DEG_HGPS_EGM_vs_control.csv"))
kid_tt    <- bind_rows(lapply(kid_genes, function(g) tt_one(Z_DCT, g, kid_meta))) %>%
  left_join(kid_deg %>% select(gene, log2FoldChange, padj), by = "gene")
Z_DCT_log <- log2(Z_DCT + 1)

# ---- liver: midzonal -----------------------------------------------
Z_liver   <- readRDS(file.path(ROOT, "liver/liver_deg_results/Z_celltype_expression.rds"))
liver_deg <- readRDS(file.path(ROOT, "liver/liver_deg_results/deg_results_list.rds"))
Z_mid     <- Z_liver$Midzonal
liv_meta  <- data.frame(Sample = colnames(Z_mid), Group = rep(c("HGPS_control", "HGPS_EGM"), each = 3))
stopifnot(all(grepl("control", liv_meta$Sample[1:3])), all(grepl("EGM", liv_meta$Sample[4:6])))
liv_genes <- c("Lcn2", "Cp", "Ppara", "Cyp8b1")
liv_tt    <- bind_rows(lapply(liv_genes, function(g) tt_one(Z_mid, g, liv_meta))) %>%
  left_join(liver_deg$Midzonal %>% select(gene, log2FoldChange, padj), by = "gene")
Z_mid_log <- log2(Z_mid + 1)

# ---- draw + save ---------------------------------------------------
draw_row(Z_DCT_log, kid_genes, kid_meta, "kidney", file.path(OUT, "Kidney_DCT_boxplot_notext.pdf"))
draw_row(Z_mid_log, liv_genes, liv_meta, "liver",  file.path(OUT, "Liver_Midzonal_boxplot_notext.pdf"))
draw_row(Z_DCT_log, kid_genes, kid_meta, "kidney", file.path(OUT, "Kidney_DCT_boxplot_ylabels.pdf"), ylabels = TRUE)
draw_row(Z_mid_log, liv_genes, liv_meta, "liver",  file.path(OUT, "Liver_Midzonal_boxplot_ylabels.pdf"), ylabels = TRUE)

stats <- bind_rows(
  kid_tt %>% mutate(Tissue = "Kidney", CellType = "DCT", p_shown = "t-test log2",
                    y_ticks = sapply(gene, function(g) y_ticks(Z_DCT_log, g, kid_meta, "kidney"))),
  liv_tt %>% mutate(Tissue = "Liver", CellType = "Midzonal", p_shown = "DESeq2 padj",
                    y_ticks = sapply(gene, function(g) y_ticks(Z_mid_log, g, liv_meta, "liver")))) %>%
  select(Tissue, CellType, gene, Mean_Control, Mean_EGM, log2FoldChange, padj_DESeq2 = padj,
         P_ttest_raw, P_ttest_log2, p_shown, y_ticks)
write.csv(stats, file.path(OUT, "added_rows_stats.csv"), row.names = FALSE)
print(as.data.frame(stats %>% mutate(across(where(is.numeric), ~signif(.x, 3)))), row.names = FALSE)
