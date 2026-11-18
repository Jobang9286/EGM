# ============================================================
# "Two-phase EGM rejuvenation cascade" from pseudotime WINDOWS
#
# Idea: ignore module labels; take the 1200 pseudotime-dynamic
# heatmap genes and rank them by how strongly they are induced inside a pseudotime window:
#       window score = mean z (bins inside window) - mean z (bins outside)
# then run preranked GSEA (clusterProfiler gseGO, GO BP) on that ranking.
#
#   left  : transition window  -> Lipid metabolic process
#   right : late window        -> Telomere maintenance
#
# Lollipop y = running enrichment score (weighted GSEA, p = 1) at each pathway gene, sorted
#           descending; x = pathway genes along the window ranking, rainbow-coloured.
# Q labels  = recomputed FDR from gseGO on the same ranking (Q_LABEL = "computed");
#             set Q_LABEL = "manuscript" to print 1.8e-02 / 6.4e-03 instead.
# Input   : data/trajectory_env.rds (z_mat: 1200 genes x 76 pseudotime bins)
# Output  : figs/Fig4l_window_lollipops.{pdf,png}, data/window_gsea_<panel>.csv
# ============================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(GO.db); library(org.Mm.eg.db); library(AnnotationDbi)
  library(ggplot2); library(patchwork)
})
ROOT <- "EGM_oracle"; set.seed(42)
Q_LABEL   <- "computed"
WIN_EARLY <- c(0.20, 0.60)     # transition window (fraction of pseudotime)
WIN_LATE  <- c(0.60, 1.00)     # late window

z  <- readRDS(file.path(ROOT, "data/trajectory_env.rds"))$z_mat; nb <- ncol(z)
win_rank <- function(w) {
  cols <- max(1, round(w[1] * nb)):round(w[2] * nb)
  sort(rowMeans(z[, cols]) - rowMeans(z[, -cols]), decreasing = TRUE)
}
go_genes <- function(term) {
  ids <- AnnotationDbi::select(GO.db, keys = term, keytype = "TERM", columns = "GOID")$GOID
  unique(AnnotationDbi::select(org.Mm.eg.db, keys = ids, keytype = "GOALL", columns = "SYMBOL")$SYMBOL)
}
es_profile <- function(rk_vec, set) {
  g <- names(rk_vec); N <- length(rk_vec); idx <- sort(match(intersect(set, g), g))
  w <- abs(as.numeric(rk_vec)); Nh <- sum(w[idx])
  hit <- numeric(N); hit[idx] <- w[idx] / Nh
  miss <- rep(1 / (N - length(idx)), N); miss[idx] <- 0
  run <- cumsum(hit) - cumsum(miss); ES <- run[which.max(abs(run))]
  if (ES < 0) run <- -run
  v <- run[idx]; ord <- order(v, decreasing = TRUE)
  data.frame(i = seq_along(idx), es = v[ord], gene = g[idx][ord], ES = abs(ES), flipped = ES < 0)
}

panels <- list(
  # windows are recorded in the header / console output; titles stay short to avoid overlap
  list(key = "lipid", term = "lipid metabolic process", title = "Lipid metabolic process",
       sub = "(transition)", win = WIN_EARLY, q_ms = 1.8e-2),
  list(key = "telo",  term = "telomere maintenance",    title = "Telomere maintenance",
       sub = "(late)",       win = WIN_LATE,  q_ms = 6.4e-3))

build <- function(p) {
  r <- win_rank(p$win)
  g <- clusterProfiler::gseGO(geneList = r, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                              minGSSize = 10, maxGSSize = 2000, pvalueCutoff = 1, eps = 0, verbose = FALSE)@result
  write.csv(g[, c("ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalue")],
            file.path(ROOT, sprintf("data/window_gsea_%s.csv", p$key)), row.names = FALSE)
  h <- g[g$Description == p$term, ]
  d <- es_profile(r, go_genes(p$term))
  cat(sprintf("%-24s window %s | pathways tested=%d | NES=%.2f Q=%.1e | genes=%d ES=%.3f%s\n",
              p$title, paste(p$win, collapse = "-"), nrow(g), h$NES, h$p.adjust, nrow(d), d$ES[1],
              ifelse(d$flipped[1], " (flipped)", "")))
  list(d = d, n_path = nrow(g), q = if (Q_LABEL == "computed") h$p.adjust else p$q_ms, meta = p)
}

PAL <- c("#D32F2F","#F57C00","#FBC02D","#9CCC65","#26A69A","#42A5F5","#3949AB")
lolli <- function(b) {
  YMAX <- c(lipid = 0.4)                                   # per-panel y-axis cap; others auto
  d <- b$d
  ymax <- if (b$meta$key %in% names(YMAX)) YMAX[[b$meta$key]] else max(0.5, ceiling(max(d$es) * 10) / 10)
  ggplot(d, aes(i, es, colour = i)) +
    geom_segment(aes(xend = i, yend = 0), linewidth = 1.05) +
    geom_point(size = 2.2) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.9) +
    scale_colour_gradientn(colours = PAL, guide = "none") +
    scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, 0.1),
                       labels = sprintf("%.1f", seq(0, ymax, 0.1)), expand = expansion(c(0, 0.02))) +
    scale_x_continuous(expand = expansion(c(0.03, 0.05))) +
    annotate("text", x = Inf, y = ymax * 0.93, size = 6, hjust = 1.15,
             label = sprintf("Q = %s", formatC(b$q, format = "e", digits = 1))) +
    labs(title = paste0(b$meta$title, "\n", b$meta$sub), x = NULL, y = "Enrichment score") +
    theme_classic(base_size = 15) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16, lineheight = 1.1),
          axis.title.y = element_text(size = 17), axis.text.y = element_text(size = 14, colour = "black"),
          axis.line.y = element_line(linewidth = 0.9), axis.ticks.y = element_line(linewidth = 0.9),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          plot.margin = margin(6, 12, 6, 6))
}

B <- lapply(panels, build)
fig <- (lolli(B[[1]]) | lolli(B[[2]])) +
  plot_annotation(title = "Two-phase EGM rejuvenation cascade",
                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20)))
ggsave(file.path(ROOT, "figs/Fig4l_window_lollipops.pdf"), fig, width = 11, height = 4.4, device = "pdf")
ggsave(file.path(ROOT, "figs/Fig4l_window_lollipops.png"), fig, width = 11, height = 4.4, dpi = 300, bg = "white")
cat("saved: figs/Fig4l_window_lollipops.{pdf,png}\n")

# ============================================================
# NES curve panels (NES on x, pathways on y, top 10 in red)
# Pathways shown = significant positively enriched GO BP terms (NES > 0, FDR < 0.05)
# of the same window ranking used for the lollipop.
# ============================================================
nes_curve <- function(b) {
  g <- read.csv(file.path(ROOT, sprintf("data/window_gsea_%s.csv", b$meta$key)))
  g <- g[g$NES > 0 & g$p.adjust < 0.05, ]
  g <- g[order(g$NES), ]; g$rank <- seq_len(nrow(g))
  tr <- which(g$Description == b$meta$term)
  cat(sprintf("%-24s NES curve: %d significant positive pathways | %s rank from top = %s\n",
              b$meta$title, nrow(g), b$meta$term, if (length(tr)) nrow(g) - tr + 1 else "not significant"))
  ticks3 <- function(lims) {                         # at most 3 NES ticks
    br <- pretty(lims, n = 2); br <- br[br >= lims[1] & br <= lims[2]]
    if (length(br) > 3) br <- br[round(seq(1, length(br), length.out = 3))]
    br
  }
  ggplot(g, aes(rank, NES)) +
    geom_point(colour = "black", size = 0.9) +        # no highlighted points
    scale_y_continuous(breaks = ticks3) +
    coord_flip() +
    labs(x = sprintf("Pathways\n(n = %d)", nrow(g)), y = "NES") +
    theme_classic(base_size = 14) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11, colour = "black"),
          axis.title = element_text(size = 15))
}
row_l <- (nes_curve(B[[1]]) | lolli(B[[1]])) + plot_layout(widths = c(1, 1.9))
row_t <- (nes_curve(B[[2]]) | lolli(B[[2]])) + plot_layout(widths = c(1, 1.9))
panelL <- row_l / row_t
ggsave(file.path(ROOT, "figs/Fig4l_panelL_full.pdf"), panelL, width = 10, height = 8.4, device = "pdf")
ggsave(file.path(ROOT, "figs/Fig4l_panelL_full.png"), panelL, width = 10, height = 8.4, dpi = 300, bg = "white")
cat("saved: figs/Fig4l_panelL_full.{pdf,png}\n")
