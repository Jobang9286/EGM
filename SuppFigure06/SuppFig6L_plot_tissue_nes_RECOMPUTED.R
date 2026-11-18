# ------------------------------------------------------------------------------
# NES bars for Liver / Kidney / Brain
#
# Input : nes_match_table.csv + tissue_gsea_<tissue>_<ranking>.csv from SuppFig6L_recompute_tissue_gsea_RECOMPUTED.R
# Ranking: lfc_all (DESeq2 on all samples of the tissue, HGPS_EGM vs HGPS_control, ranked by log2FC).
#   Liver "Cell development" has a negative NES with this ranking,
#   so it is replaced by the closest development term that goes up:
#   "positive regulation of cell development" (GO:0010720, NES 1.26, P = 0.0014, FDR 0.098).
# Layout: three equal panels (0.45 x 0.63 in each) side by side,
#   tissue strip on top, NES title on the left, term names under the bars reading bottom-to-top,
#   4 pt tick/term text, 5 pt strips. Bars run from the highest NES (left, red) to the lowest (right, pink)
#   within each tissue, so after rotating the panel 90 degrees clockwise the top bar is the
#   highest and reddest.
# Outputs:
#   Supp7I_NES.pdf / _1200dpi.png         full panel with all text
#   Supp7I_NES_notext_{Liver,Kidney,Brain}.pdf / _1200dpi.png
#                                          bars + frame + ticks only
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr); library(grid) })
OUT     <- "4FP_plasma/EGM_SuppFig7_NES"
RANKING <- Sys.getenv("NES_RANKING", "lfc_all")
FONT    <- "Arial"

tab <- read.csv(file.path(OUT, "nes_match_table.csv")) %>% filter(ranking == RANKING)
stopifnot(nrow(tab) == 13)
SWAP <- c(tissue = "Liver", old = "Cell development", term = "positive regulation of cell development",
          label = "Cell development")   # shown as "Cell development"
new <- read.csv(file.path(OUT, sprintf("tissue_gsea_Liver_%s.csv", RANKING))) %>% filter(Description == SWAP[["term"]])
stopifnot(nrow(new) == 1, new$NES > 0)
i <- which(tab$tissue == SWAP[["tissue"]] & tab$label == SWAP[["old"]])
tab[i, c("term", "label", "setSize", "NES", "p.adjust", "NES_ppt")] <-
  list(SWAP[["term"]], SWAP[["label"]], new$setSize, new$NES, new$p.adjust, NA)   # no NES_ppt value for the new term
pv  <- bind_rows(lapply(unique(tab$tissue), function(t)
         read.csv(file.path(OUT, sprintf("tissue_gsea_%s_%s.csv", t, RANKING))) %>% transmute(tissue = t, term = Description, pvalue)))
dropped <- tab %>% filter(!is.na(NES_ppt), sign(NES) != sign(NES_ppt))
tab <- tab %>% left_join(pv, by = c("tissue", "term")) %>%
  group_by(tissue) %>% arrange(desc(NES), .by_group = TRUE) %>%
  mutate(fill_col = colorRampPalette(c("#FF0000", "#FFC0CB"))(n())) %>% ungroup() %>%
  mutate(label = factor(label, levels = unique(label)))
write.csv(tab %>% select(tissue, label, term, setSize, NES_figure = NES_ppt, NES, pvalue, p.adjust),
          file.path(OUT, sprintf("Supp7I_NES_values_%s.csv", RANKING)), row.names = FALSE)
y_lim <- list(Liver = c(0, 2), Kidney = c(0, 1.5), Brain = c(-2, 2))

panel <- function(tis, text = TRUE, show_title = FALSE) {
  d <- tab %>% filter(tissue == tis)
  ggplot(d, aes(label, NES, fill = fill_col)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.25) +
    geom_col(width = 0.8, colour = "black", linewidth = 0.3) +
    scale_fill_identity() +
    scale_y_continuous(limits = y_lim[[tis]], breaks = seq(-2, 2, 0.5), expand = expansion(mult = 0.02)) +
    facet_grid(. ~ tissue) +
    labs(x = NULL, y = if (show_title) "Normalized Enrichment Score (NES)" else NULL) +
    theme_bw(base_size = 4, base_family = FONT) +
    theme(panel.border = element_rect(colour = "black", linewidth = 0.5),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
          panel.grid.minor = element_line(colour = "grey92", linewidth = 0.25),
          axis.ticks = element_line(colour = "black", linewidth = 0.3), axis.ticks.length = unit(1, "pt"),
          axis.text  = element_text(colour = "black", size = 4),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.y = element_text(size = 4.5, margin = margin(r = 1)),
          strip.background = element_rect(fill = "grey85", colour = NA),
          strip.text = element_text(face = "bold", size = 5, margin = margin(1.5, 0, 1.5, 0)),
          plot.margin = margin(0.5, 1.5, 0.5, 0.5)) +
    if (!text) theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.x = element_blank(),
                     strip.background = element_blank(), strip.text = element_blank(),
                     axis.title = element_blank(), plot.margin = margin(0.4, 0.4, 0.4, 0.4))   # keeps the frame unclipped
}

open_pdf <- function(f, w, h) quartz(type = "pdf", file = f, width = w, height = h, family = FONT)
png1200  <- function(f) system2("pdftoppm", c("-png", "-r", "1200", "-singlefile", shQuote(f), shQuote(sub("\\.pdf$", "_1200dpi", f))))

# full panel: 2.31 x 1.72 in
f <- file.path(OUT, "Supp7I_NES.pdf")
open_pdf(f, 2.31, 1.72)
print(panel("Liver", show_title = TRUE) | panel("Kidney") | panel("Brain"))
dev.off(); png1200(f)

# the same panel already turned (terms on the left, highest NES at the top, strips on the right,
# NES axis at the bottom), rotated size: 1.72 x 2.31 in
panel_h <- function(tis, show_title = FALSE, show_terms = TRUE, narrow = FALSE, gap_below = 0, xlab = NULL) {
  d <- tab %>% filter(tissue == tis) %>% mutate(label = factor(label, levels = rev(levels(droplevels(label)))))
  lw <- if (narrow) c(bar = 0.2, frame = 0.35, tick = 0.2, zero = 0.2) else c(bar = 0.3, frame = 0.5, tick = 0.3, zero = 0.25)
  step <- if (narrow) c(Liver = 1, Kidney = 0.5, Brain = 2)[[tis]] else if (tis == "Brain") 1 else 0.5   # labels must not touch
  ggplot(d, aes(NES, label, fill = fill_col)) +
    geom_vline(xintercept = 0, colour = "black", linewidth = lw[["zero"]]) +
    geom_col(width = 0.8, colour = "black", linewidth = lw[["bar"]], orientation = "y") +
    scale_fill_identity() +
    scale_x_continuous(limits = y_lim[[tis]], breaks = seq(-2, 2, step),
                       labels = function(b) sprintf("%.1f", b), expand = expansion(mult = 0.02)) +
    facet_grid(tissue ~ .) +
    labs(x = xlab, y = NULL, title = if (show_title) "Normalized Enrichment Score (NES)" else NULL) +
    theme_bw(base_size = 4, base_family = FONT) +
    theme(panel.border = element_rect(colour = "black", linewidth = lw[["frame"]]),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
          panel.grid.minor = element_line(colour = "grey92", linewidth = 0.25),
          axis.ticks = element_line(colour = "black", linewidth = lw[["tick"]]), axis.ticks.length = unit(1, "pt"),
          axis.text = element_text(colour = "black", size = 4),
          plot.title = element_text(size = 4.5, hjust = 0.5, margin = margin(b = 1.5)),
          axis.title.x = element_text(size = 4, margin = margin(t = 1)),
          strip.background = element_rect(fill = "grey85", colour = NA),
          strip.text = element_text(face = "bold", size = if (narrow) 4.3 else 5, margin = margin(0, 1.5, 0, 1.5)),
          plot.margin = margin(0.5, 0.5, 1.5 + gap_below, if (narrow) 4.5 else 0.5)) +   # left: room for half of "-2.0"
    if (!show_terms) theme(axis.text.y = element_blank())
}
f <- file.path(OUT, "Supp7I_NES_horizontal.pdf")
open_pdf(f, 1.72, 2.31)
print(panel_h("Liver", show_title = TRUE) / panel_h("Kidney") / panel_h("Brain"))
dev.off(); png1200(f)

# narrow version: 4.5 pt extra space under Liver and Kidney, panels 0.59 in tall,
# thinner lines, "NES" (4 pt) under the Brain axis.
#   _noterms : 0.79 x 2.27 in, no term names
#   _terms   : same panel size, term names (4 pt) on the left; the figure is wider by the longest name
narrow_fig <- function(show_terms) {
  (panel_h("Liver",  show_terms = show_terms, narrow = TRUE, gap_below = 4.5) /
   panel_h("Kidney", show_terms = show_terms, narrow = TRUE, gap_below = 4.5) /
   panel_h("Brain",  show_terms = show_terms, narrow = TRUE, xlab = "NES")) +
    plot_annotation(theme = theme(plot.margin = margin(0, 0, 0.5, 0)))   # no outer patchwork margin
}
f <- file.path(OUT, "Supp7I_NES_horizontal_noterms.pdf")
NARROW_H <- 2.27   # panels 0.59 in + the wider gaps
NARROW_W <- 0.79
open_pdf(f, NARROW_W, NARROW_H); print(narrow_fig(FALSE)); dev.off(); png1200(f)

f <- file.path(OUT, "Supp7I_NES_horizontal_terms.pdf")
open_pdf(f, 1, 1)   # measure the longest term name in the device font first
lab_w <- max(sapply(as.character(tab$label), function(l) convertWidth(grobWidth(textGrob(l,
           gp = gpar(fontsize = 4, fontfamily = FONT))), "in", valueOnly = TRUE)))
dev.off()
# 0.85 pt: measured so both versions have the same panel width
open_pdf(f, NARROW_W + lab_w + 0.85 / 72, NARROW_H); print(narrow_fig(TRUE)); dev.off(); png1200(f)

# bars only, one file per tissue
for (tis in names(y_lim)) {
  f <- file.path(OUT, sprintf("Supp7I_NES_notext_%s.pdf", tis))
  open_pdf(f, 0.449, 0.630)
  print(panel(tis, text = FALSE))
  dev.off(); png1200(f)
}
cat("ranking:", RANKING, "| sign differs from the figure:", paste(dropped$tissue, dropped$label, collapse = ", "), "\n")
