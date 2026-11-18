# ------------------------------------------------------------------
# NES of selected pathways across pseudotime quintiles (Q1 -> Q5) of the rejuvenation trajectory
# Input : quintile_gsea_Q1..Q5.csv from Fig3L_quintile_window_GSEA.R (window GSEA, each quintile vs the rest)
# Curves: shape-preserving cubic interpolation (PCHIP, Fritsch-Carlson) through the five NES values,
#         so a curve never overshoots its points; area between curve and 0 is filled.
# Points: filled = FDR < 0.05 in that quintile, open = n.s.
# Output: quintile_NES_curves.pdf / .png
# ------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(grid) })

qd_dir <- "4FP_plasma/EGM_Fig3K_pathways"
qd_q   <- paste0("Q", 1:5)

qd_sel <- tibble::tribble(
  ~Description,                              ~group,
  "lipid metabolic process",                 "Q2-Q3",   # the two lipid terms with FDR < 0.05 in Q2
  "monocarboxylic acid metabolic process",   "Q2-Q3",
  "cell redox homeostasis",                  "Q2-Q3",   # sharp, Q2-specific rise (FDR 0.032)
  "response to fibroblast growth factor",    "Q4",
  "cell junction organization",              "Q4",
  "morphogenesis of an epithelial sheet",    "Q4",
  "telomere maintenance",                    "Q5",
  "RNA splicing",                            "Q5",
  "chromatin remodeling",                    "Q5")
qd_pal   <- list("Q2-Q3" = c("#D62839", "#F18F01", "#7A1E48"),
                 "Q4"    = c("#2A9D8F", "#57CC99", "#1B6F66"),
                 "Q5"    = c("#6A4C93", "#3A86FF", "#B07CC6"))
qd_band  <- list("Q2-Q3" = c(1.5, 3.5), "Q4" = c(3.5, 4.5), "Q5" = c(4.5, 5.2))   # x axis runs 0.8-5.2
qd_title <- c("Q2-Q3" = "Early transition (Q2-Q3)", "Q4" = "Mid-late (Q4)", "Q5" = "Late (Q5)")

qd_pts <- bind_rows(lapply(qd_q, function(w)
            read.csv(file.path(qd_dir, sprintf("quintile_gsea_%s.csv", w))) %>% mutate(window = w))) %>%
  filter(Description %in% qd_sel$Description) %>% left_join(qd_sel, by = "Description") %>%
  mutate(x = as.integer(sub("Q", "", window)), sig = p.adjust < 0.05,
         label = paste0(toupper(substr(Description, 1, 1)), substring(Description, 2)))
stopifnot(nrow(qd_pts) == nrow(qd_sel) * 5)

qd_pchip <- function(x, y, n = 200) {      # slopes are 0 at local extrema -> no overshoot between points
  h <- diff(x); d <- diff(y) / h; k <- length(x); m <- numeric(k)
  m[1] <- d[1]; m[k] <- d[k - 1]
  for (i in 2:(k - 1)) m[i] <- if (d[i - 1] * d[i] <= 0) 0 else {
    w1 <- 2 * h[i] + h[i - 1]; w2 <- h[i] + 2 * h[i - 1]; (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]) }
  xs <- seq(min(x), max(x), length.out = n); j <- pmin(findInterval(xs, x, rightmost.closed = TRUE), k - 1)
  t <- (xs - x[j]) / h[j]
  data.frame(x = xs, y = (2*t^3 - 3*t^2 + 1) * y[j] + (t^3 - 2*t^2 + t) * h[j] * m[j] +
                          (-2*t^3 + 3*t^2) * y[j + 1] + (t^3 - t^2) * h[j] * m[j + 1])
}
qd_curve <- qd_pts %>% group_by(group, label) %>% arrange(x) %>% reframe(qd_pchip(x, NES))
qd_ylim  <- range(qd_pts$NES) + c(-0.25, 0.3)

# k scales line widths, point sizes and legend keys, kt the text and theme spacing, kl the legend text
# (k = kt = kl = 1: the 10.07 in figure); xgap = space between the x tick labels and "Pseudotime quintile",
# in pt (theme default 2.75 * kt); legend names longer than wrap characters go on two lines, except those
# given in lab_break (name -> label with its own "\n")
qd_panel <- function(gp, show_y = TRUE, show_shape = FALSE, show_title = TRUE, k = 1, kt = k, kl = kt,
                     xgap = 2.75 * kt, wrap = Inf, lab_break = NULL, family = "") {
  cv <- qd_curve %>% filter(group == gp); pt <- qd_pts %>% filter(group == gp)
  labs_lv <- qd_sel$Description[qd_sel$group == gp]
  labs_lv <- paste0(toupper(substr(labs_lv, 1, 1)), substring(labs_lv, 2))
  cols <- setNames(qd_pal[[gp]][seq_along(labs_lv)], labs_lv)
  ggplot() +
    annotate("rect", xmin = qd_band[[gp]][1], xmax = qd_band[[gp]][2], ymin = -Inf, ymax = Inf, fill = "grey94") +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35 * k) +
    geom_ribbon(data = cv, aes(x, ymin = pmin(y, 0), ymax = pmax(y, 0), fill = label), alpha = 0.13) +
    geom_line(data = cv, aes(x, y, colour = label), linewidth = 0.95 * k) +
    geom_point(data = pt, aes(x, NES, colour = label, shape = sig), size = 2.3 * k, fill = "white", stroke = 0.85 * k) +
    scale_colour_manual(values = cols, breaks = labs_lv, name = NULL, labels = function(x) vapply(x, function(l)
      if (l %in% names(lab_break)) lab_break[[l]] else if (nchar(l) > wrap) paste(strwrap(l, wrap + 1), collapse = "\n")
      else l, "")) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."), name = NULL) +
    scale_x_continuous(breaks = 1:5, labels = qd_q, expand = expansion(add = 0.2)) +
    scale_y_continuous(breaks = seq(-2, 4, 2)) +
    coord_cartesian(ylim = qd_ylim) +
    labs(title = if (show_title) qd_title[[gp]] else NULL, x = "Pseudotime quintile", y = if (show_y) "NES" else NULL) +
    guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 2 * k)),
           shape = if (show_shape) guide_legend(order = 2, override.aes = list(colour = "grey25")) else "none") +
    theme_classic(base_size = 11 * kt, base_family = family) +
    theme(plot.title = element_text(face = "bold", size = 11 * kt), axis.text = element_text(colour = "black"),
          axis.title.x = element_text(margin = margin(t = xgap)),
          legend.key.size = unit(1.2 * kt, "lines"),   # theme default, but "lines" does not follow base_size
          legend.key.width = unit(0.7 * k, "cm"), legend.text = element_text(size = 9.5 * kl),
          legend.spacing.y = unit(1 * k, "pt"),
          # one row of panels: each legend sits under its panel, aligned to the left edge of the plot and to the top
          legend.position = "bottom", legend.direction = "vertical", legend.box = "horizontal",
          legend.justification = c(0, 1), legend.location = "plot", legend.margin = margin(2 * k, 0, 0, 0))
}
qd_fig <- qd_panel("Q2-Q3") | qd_panel("Q4", show_y = FALSE) | qd_panel("Q5", show_y = FALSE, show_shape = TRUE)

ggsave(file.path(qd_dir, "quintile_NES_curves.pdf"), qd_fig, width = 10.07, height = 2.98, device = cairo_pdf)
ggsave(file.path(qd_dir, "quintile_NES_curves.png"), qd_fig, width = 10.07, height = 2.98, dpi = 1200, bg = "white")
print(qd_fig)

# Fig 3L: no panel titles, drawn at the final figure size (1200 dpi).
# Lines and points use k = 3.62 / 9.775 (3.62 in of content over the content width of
# the 10.07 in version). Text is 15% larger than that (kt = 1.15 k), "Pseudotime quintile" sits 2.2 pt below
# the tick labels, the legend text is 20% larger again (kl = 1.2 kt) with names over 26 characters on two
# lines and 2 pt between terms, the Late legend has its n.s. / FDR keys in one row under the terms, and each
# panel is fixed at 0.978 x 0.656 in (10% narrower than when the panels filled the 3.62 in width); the figure is
# then exactly as large as its content needs.
qd_k   <- 3.62 / 9.775
qd_kt  <- 1.15 * qd_k
qd_pw  <- 0.978; qd_ph <- 0.656
qd_gap <- 2   # pt between legend terms, and above the n.s. / FDR row
qd_brk <- c("Morphogenesis of an epithelial sheet" = "Morphogenesis of\nan epithelial sheet")   # break after "of"
qd_p3L <- lapply(list(list("Q2-Q3", TRUE, FALSE), list("Q4", FALSE, FALSE), list("Q5", FALSE, TRUE)), function(a)
  qd_panel(a[[1]], show_y = a[[2]], show_shape = a[[3]], show_title = FALSE, k = qd_k, kt = qd_kt, kl = 1.2 * qd_kt,
           xgap = 2.2, wrap = 26, lab_break = qd_brk, family = "Arial") + theme(legend.key.spacing.y = unit(qd_gap, "pt")))
qd_p3L[[3]] <- qd_p3L[[3]] + guides(shape = guide_legend(order = 2, nrow = 1, override.aes = list(colour = "grey25"))) +
  theme(legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(qd_gap, "pt"))
qd_f <- file.path(qd_dir, "Fig3L_quintile_NES_curves.pdf")
# measure on the PNG device: cairo_pdf rounds glyph widths at 4 pt and under-measures text by ~9%,
# so the PDF is written with the quartz device, whose Arial widths match ragg's
ragg::agg_png(tempfile(fileext = ".png"), 10, 10, units = "in", res = 1200)
in_w <- function(u) convertWidth(u, "in", valueOnly = TRUE)
# leave room on the right in case the Late legend is wider than its plot
qd_g3   <- ggplotGrob(qd_p3L[[3]])
qd_over <- max(0, in_w(sum(qd_g3$grobs[[which(qd_g3$layout$name == "guide-box-bottom")]]$widths)) -
                  in_w(qd_g3$widths[qd_g3$layout$l[qd_g3$layout$name == "axis-l"]]) - qd_pw -
                  in_w(unit(5.5 * qd_kt, "pt")))
qd_fig_3L <- (qd_p3L[[1]] | qd_p3L[[2]] | qd_p3L[[3]]) +
  plot_layout(widths = unit(qd_pw, "in"), heights = unit(qd_ph, "in")) +
  plot_annotation(theme = theme(plot.margin = margin(0, qd_over * 72 + 1.5, 0, 0)))   # no other outer margin
qd_g <- patchwork::patchworkGrob(qd_fig_3L)
qd_w <- in_w(sum(qd_g$widths))
qd_h <- convertHeight(sum(qd_g$heights), "in", valueOnly = TRUE)
dev.off()
ggsave(qd_f, qd_fig_3L, width = qd_w, height = qd_h,
       device = function(filename, ...) quartz(type = "pdf", file = filename, family = "Arial", ...))
ggsave(sub("\\.pdf$", "_1200dpi.png", qd_f), qd_fig_3L, width = qd_w, height = qd_h, dpi = 1200, bg = "white",
       device = ragg::agg_png)
