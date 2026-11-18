# H3K27me3 coverage at Cdkn2a and the HoxA cluster, P-BP over CP.
#
# Both arms share a y scale within a locus, so bar heights are directly
# comparable. Signal is arm-mean CPM in 200 bp bins (see 01_upstream).
suppressMessages({library(tidyverse); library(patchwork); library(ggrepel)})
# project root: the directory holding data/ and figures/. Works when the script
# is run from the project root, from scripts/, or via Rscript with a full path.
find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  cands <- c(if (length(f)) dirname(dirname(normalizePath(f[1]))), getwd(), dirname(getwd()))
  hit <- cands[dir.exists(file.path(cands, "data"))]
  if (!length(hit)) stop("cannot locate the data/ directory; run from the project root")
  hit[1]
}
here <- find_root()
DATA <- file.path(here, "data"); OUT <- file.path(here, "figures")
dir.create(OUT, showWarnings = FALSE)

PAL  <- c(`P-BP` = "#5FA8CE", CP = "#ADADAD")
OUTL <- c(`P-BP` = "#2E7CA6", CP = "#7A7A7A")
# round the top tick up so a low-signal panel still gets a labelled maximum
nice_top <- function(x) { p <- 10^floor(log10(x)); ceiling(x / p * 2) / 2 * p }

# GENCODE vM38 coordinates
gm <- tribble(
  ~locus, ~gene, ~start, ~end,
  "Cdkn2a_locus", "Cdkn2a", 89192707, 89212890,
  "Cdkn2a_locus", "Cdkn2b", 89224535, 89229276,
  "HoxA_cluster", "Hoxa1",  52132569, 52135297,
  "HoxA_cluster", "Hoxa5",  52178733, 52181567,
  "HoxA_cluster", "Hoxa9",  52200079, 52208069,
  "HoxA_cluster", "Hoxa10", 52208176, 52217834,
  "HoxA_cluster", "Hoxa13", 52234673, 52237788)

d <- read_tsv(file.path(DATA, "locus_signal_H3K27me3.tsv"), show_col_types = FALSE) |>
  pivot_longer(c(sig_4FP, sig_CP), names_to = "arm", values_to = "v") |>
  mutate(arm = factor(if_else(arm == "sig_4FP", "P-BP", "CP"), c("P-BP", "CP")),
         mid = (start + end) / 2) |>
  group_by(locus) |> mutate(ymax = max(v)) |> ungroup()

pan <- function(loc, ttl) {
  dd  <- filter(d, locus == loc); rng <- range(c(dd$start, dd$end)) / 1e6
  top <- nice_top(max(dd$v)); gg <- filter(gm, locus == loc)
  trk <- ggplot(dd, aes(mid / 1e6, v, fill = arm, colour = arm)) +
    geom_area(linewidth = .2) +
    facet_grid(arm ~ ., switch = "y") +
    scale_fill_manual(values = PAL, guide = "none") +
    scale_colour_manual(values = OUTL, guide = "none") +
    scale_x_continuous(limits = rng, expand = expansion(mult = .012)) +
    scale_y_continuous(limits = c(0, top * 1.05), breaks = c(0, top), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = ttl) +
    theme_classic(base_size = 12) +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 10, face = "bold",
                                           hjust = 1, margin = margin(r = 5)),
          axis.line.y = element_line(colour = "black", linewidth = .4),
          axis.line.x = element_blank(),
          axis.ticks.y = element_line(colour = "black", linewidth = .4),
          axis.ticks.length = unit(3, "pt"),
          axis.text.y = element_text(size = 8.5, colour = "black"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = 13, face = "bold", hjust = .5, margin = margin(b = 3)),
          panel.spacing.y = unit(1.6, "mm"), plot.margin = margin(3, 18, 1, 4))
  ann <- ggplot(gg) +
    geom_segment(aes(x = start / 1e6, xend = end / 1e6, y = 0, yend = 0),
                 linewidth = 3.2, colour = "grey25") +
    geom_text_repel(aes(x = (start + end) / 2e6, y = 0, label = gene), size = 3.1,
                    colour = "grey10", nudge_y = -1, direction = "x", segment.size = .28,
                    segment.colour = "grey60", min.segment.length = 0,
                    box.padding = .13, max.overlaps = Inf, seed = 1) +
    scale_x_continuous(limits = rng, expand = expansion(mult = .012)) +
    scale_y_continuous(limits = c(-2, .5), breaks = NULL, expand = c(0, 0)) +
    labs(x = paste0(dd$chrom[1], " (Mb)"), y = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(colour = "black", linewidth = .4),
          axis.ticks.x = element_line(colour = "black", linewidth = .4),
          axis.ticks.length = unit(3, "pt"),
          axis.text.x = element_text(size = 8.5, colour = "black"),
          axis.title.x = element_text(size = 9.5, margin = margin(t = 2)),
          plot.margin = margin(1, 18, 3, 4))
  wrap_elements(trk / ann + plot_layout(heights = c(2.6, 1.6)))
}

p <- pan("Cdkn2a_locus", "Cdkn2a") / pan("HoxA_cluster", "HoxA cluster")
ggsave(file.path(OUT, "locus_tracks_H3K27me3_Cdkn2a_HoxA.pdf"), p, width = 6.0, height = 4.2)
ggsave(file.path(OUT, "locus_tracks_H3K27me3_Cdkn2a_HoxA.png"), p, width = 6.0, height = 4.2, dpi = 300)
message("wrote locus_tracks_H3K27me3_Cdkn2a_HoxA.{pdf,png}")
