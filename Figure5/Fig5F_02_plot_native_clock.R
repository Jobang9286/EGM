# Fig 5F — native macaque transcriptAge clock, Brain / Liver / Kidney.
# Annotations: MAE / R (LOO, from native_clock_ymctrl_summary.csv) and ΔAgeDiff = clock(EGM) - clock(Saline).
# Output: outputs/Fig5F_native_clock_3tis.pdf / _1200dpi.png / _1200dpi_transparent.png
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
base <- "monkey/reproduce data/repro_native_clock_egm"
out  <- file.path(base, "outputs")
tis  <- c("Brain", "Liver", "Kidney")

pp <- fread(file.path(out, "native_clock_ymctrl_pubpred.csv"))    # Y / M / Saline LOO clock
tp <- fread(file.path(out, "native_clock_ymctrl_testpred.csv"))   # EGM projected clock
sm <- fread(file.path(out, "native_clock_ymctrl_summary.csv"))
pred  <- rbind(pp[, .(tissue, group, clock)], tp[, .(tissue, group, clock)])[tissue %in% tis]
labs4 <- c(Y = "Y-Ctrl", M = "M-Ctrl", CTRL = "Saline", EGM = "EGM")
xpos  <- c(Y = 5, M = 11, CTRL = 16, EGM = 16)    # x compressed: 22 yr drawn at 16
chron <- c(Y = 5, M = 11, CTRL = 22, EGM = 22)
gcol  <- c("Y-Ctrl" = "#C9C9C9", "M-Ctrl" = "#4472C4", "Saline" = "#C00000", "EGM" = "#70AD47")
pred[, `:=`(grp = factor(labs4[group], levels = names(gcol)), x = xpos[group], dAge = clock - chron[group],
            tissue = factor(tissue, levels = tis))]

# ---- sizes (pt / in) ----------------------------------------------------------
W <- 4.62; H <- 2.35
FONT <- "Arial"
FS_TITLE <- 7; FS_SUB <- 6.5; FS_AX <- 6; FS_TICK <- 5.5; FS_ANN <- 5.5; FS_LEG <- 4.5   # FS_LEG: legend inside the panels
LW <- 0.25; LW_BOX <- 0.35; LW_FIT <- 0.4; PT <- 0.9; PT_B <- 0.7
ROW_GAP <- 7                                      # extra pt between the scatter row and the boxplot row

fillbrk <- function(x) { R <- diff(x); raw <- R / 3; cand <- c(1, 2, 3, 5, 10, 20); ok <- cand[cand <= raw]
  step <- if (length(ok)) max(ok) else 1
  sort(unique(c(seq(0, ceiling(x[2] / step) * step, by = step), seq(0, floor(x[1] / step) * step, by = -step)))) }

thm <- theme_classic(base_size = FS_AX, base_family = FONT, base_line_size = LW) +
  theme(axis.text = element_text(colour = "black", size = FS_TICK), axis.title = element_text(size = FS_AX),
        axis.ticks.length = unit(1.5, "pt"), plot.title = element_text(size = FS_TITLE, hjust = 0.5, margin = margin(b = 2)),
        plot.margin = margin(2, 5, 2, 2))

# legend key = box + dot, the look of the boxplot row it also stands for
key_box_point <- function(data, params, size)
  grid::grobTree(draw_key_boxplot(transform(data, fill = NA, linewidth = LW_BOX * 0.7), params, size),
                 draw_key_point(transform(data, size = 0.6), params, size))

mktop <- function(t) {
  d <- pred[tissue == t]; s <- sm[tissue == t]
  ggplot(d, aes(x, clock)) +
    geom_smooth(data = d[grp != "EGM"], method = "lm", formula = y ~ x, se = TRUE, colour = "black",
                fill = "grey85", fullrange = TRUE, linewidth = LW_FIT) +
    geom_point(aes(colour = grp), size = PT, key_glyph = key_box_point) +
    annotate("text", x = 0.6, y = Inf, hjust = 0, vjust = 1.2, lineheight = 0.95, size = FS_ANN / .pt, family = FONT,
             label = sprintf("MAE = %.2f\nR = %.3f", s$LOO_MAE, s$LOO_R)) +
    scale_colour_manual(values = gcol, name = NULL) +
    scale_x_continuous(limits = c(0, 17.5), expand = expansion(mult = c(0, 0.02)), breaks = c(0, 5, 10, 16),
                       labels = c("0", "5", "10", "22")) +
    scale_y_continuous(breaks = function(x) seq(0, floor(x[2]), by = 10)) +
    labs(x = "Chronological age (year)", y = "transcriptAge (year)", title = t) + thm +
    theme(plot.margin = margin(2, 5, 2 + ROW_GAP, 2),
          # legend inside each scatter, bottom right (empty there: the 22 yr points sit above 14)
          legend.position = "inside", legend.position.inside = c(1, 0.02), legend.justification = c(1, 0),
          legend.key.size = unit(5, "pt"), legend.key.spacing.y = unit(0.4, "pt"), legend.text = element_text(size = FS_LEG),
          legend.background = element_blank(), legend.margin = margin(0, 0, 0, 0))
}

mkbx <- function(t) {
  d <- pred[tissue == t]; s <- sm[tissue == t]
  ggplot(d, aes(grp, dAge)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey55", linewidth = LW) +
    geom_boxplot(aes(colour = grp), fill = NA, outlier.shape = NA, width = 0.6, linewidth = LW_BOX) +
    geom_point(aes(colour = grp), position = position_jitter(width = 0.1, height = 0, seed = 1), size = PT_B) +
    scale_colour_manual(values = gcol, guide = "none") +
    scale_y_continuous(breaks = fillbrk, expand = expansion(mult = c(0.06, 0.1))) +
    labs(x = NULL, y = "\u0394Age (year)", title = sprintf("\u0394AgeDiff = %.2f yr", s$EGM_minus_CTRL)) + thm +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), plot.title = element_text(size = FS_SUB))
}

fig <- wrap_plots(c(lapply(tis, mktop), lapply(tis, mkbx)), ncol = 3, heights = c(1.12, 1))

draw_fig <- function() print(fig)

f <- file.path(out, "Fig5F_native_clock_3tis")
quartz(type = "pdf", file = paste0(f, ".pdf"), width = W, height = H, family = FONT); draw_fig(); invisible(dev.off())
ragg::agg_png(paste0(f, "_1200dpi.png"), width = W, height = H, units = "in", res = 1200, background = "white")
draw_fig(); invisible(dev.off())

# the same figure on a transparent background (every background rect removed, patchwork's own included)
none <- element_rect(fill = "transparent", colour = NA)
fig_t <- (fig & theme(plot.background = none, panel.background = none, legend.background = none,
                      legend.box.background = none, legend.key = none)) +
  plot_annotation(theme = theme(plot.background = none))
ragg::agg_png(paste0(f, "_1200dpi_transparent.png"), width = W, height = H, units = "in", res = 1200,
              background = "transparent")
print(fig_t); invisible(dev.off())
