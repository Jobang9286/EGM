# ------------------------------------------------------------------------------
# Fig 5E — monkey blood chemistry (CHEM) and complete blood count (CBC), Young / Saline / EGM
# Data : CHEM, CBC data.xlsx (sheets CHEM, CBC; one value per group)
# Colour = z-score of each parameter across the 3 groups (scale(), so with n = 3 the largest |z| is 1.155)
# Drawn as one figure (2.35 x 1.34 in) with titles, group and parameter names and the Z-score bar
# Output: Fig5E_CHEM_CBC.pdf, _1200dpi.png, _1200dpi_transparent.png
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(readxl); library(ggplot2); library(patchwork) })
OUT  <- "4FP_plasma/EGM_Fig5E_CHEM_CBC"
XLSX <- "CHEM, CBC data.xlsx"

rd <- function(s) {
  x <- suppressMessages(read_excel(XLSX, sheet = s))
  names(x)[1:4] <- c("no", "parameter", "abbr", "unit")
  x <- x[!is.na(x$parameter) & !is.na(x$abbr), ]
  for (g in c("Young", "Saline", "EGM")) x[[g]] <- suppressWarnings(as.numeric(x[[g]]))
  x$sheet <- s
  x[, c("sheet", "parameter", "abbr", "unit", "Young", "Saline", "EGM")]
}
d <- rbind(rd("CHEM"), rd("CBC"))
stopifnot(sum(d$sheet == "CHEM") == 19, sum(d$sheet == "CBC") == 21, !anyNA(d[, c("Young", "Saline", "EGM")]))
zs <- t(apply(as.matrix(d[, c("Young", "Saline", "EGM")]), 1, scale))
colnames(zs) <- c("Young", "Saline", "EGM")
long <- data.frame(sheet = rep(d$sheet, 3), abbr = rep(d$abbr, 3),
                   grp = factor(rep(colnames(zs), each = nrow(d)), rev(colnames(zs))), z = as.vector(zs))
ZM  <- 1.155
LOW <- "#4A7DB5"; HIGH <- "#C4443C"

# ---- sizes (in / pt) ---------------------------------------------------------
W <- 2.35; H <- 1.34
FONT <- "Arial"
FS_TITLE <- 7; FS_GRP <- 5.5; FS_PAR <- 4.5; FS_LEG <- 4.5
LW_TILE <- 0.15
LEG_TOP <- 9        # pt above the CHEM title, so the Z-score title over the bar fits

blk <- function(sh, legend = FALSE) {
  L <- long[long$sheet == sh, ]
  L$abbr <- factor(L$abbr, unique(d$abbr[d$sheet == sh]))
  ggplot(L, aes(abbr, grp, fill = z)) +
    geom_tile(colour = "black", linewidth = LW_TILE) +
    scale_fill_gradient2(low = LOW, mid = "white", high = HIGH, midpoint = 0, limits = c(-ZM, ZM),
                         breaks = c(-1, 0, 1), name = "Z-score") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(title = sh, x = NULL, y = NULL) +
    theme_minimal(base_size = FS_GRP, base_family = FONT) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = FS_PAR, margin = margin(t = 1)),
          axis.text.y = element_text(size = FS_GRP, margin = margin(r = 1)),
          plot.title = element_text(face = "bold", size = FS_TITLE, hjust = 0.5, margin = margin(b = 1.5)),
          plot.margin = margin(if (legend) LEG_TOP else 1, 1, 2, 1),   # room above CHEM for the taller legend
          # the Z-score bar sits in the title band, right of "CHEM" (bottom-right corner on the panel's top-right)
          legend.position = if (legend) "inside" else "none", legend.position.inside = c(1, 1.13),
          legend.justification.inside = c(1, 0), legend.direction = "horizontal",
          legend.title = element_text(size = FS_LEG + 0.5, hjust = 0.5, margin = margin(b = 1)),
          legend.title.position = "top",
          legend.text = element_text(size = FS_LEG, margin = margin(t = 0.5)),
          legend.key.width = unit(0.3 / 5, "in"),    # a colourbar is 5 keys long: 0.3 in bar
          legend.key.height = unit(0.045, "in"), legend.margin = margin(0, 0, 0, 0),
          legend.background = element_blank())
}

fig <- blk("CHEM", legend = TRUE) / blk("CBC") +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))   # no patchwork outer margin

f <- file.path(OUT, "Fig5E_CHEM_CBC")
quartz(type = "pdf", file = paste0(f, ".pdf"), width = W, height = H, family = FONT); print(fig); invisible(dev.off())
ragg::agg_png(paste0(f, "_1200dpi.png"), width = W, height = H, units = "in", res = 1200, background = "white")
print(fig); invisible(dev.off())
none <- element_rect(fill = "transparent", colour = NA)
fig_t <- (fig & theme(plot.background = none, panel.background = none, legend.background = none)) +
  plot_annotation(theme = theme(plot.background = none))
ragg::agg_png(paste0(f, "_1200dpi_transparent.png"), width = W, height = H, units = "in", res = 1200,
              background = "transparent")
print(fig_t); invisible(dev.off())
