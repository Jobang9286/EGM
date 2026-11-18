# ---- DESeq2 ----
suppressMessages(library(DESeq2))
# ---- ggplot2 / dplyr ----
suppressMessages({library(ggplot2); library(dplyr)})
# ---- gene lists, condition vectors, count matrices ----
suppressMessages(library(DEGreport)); suppressMessages(library(ggsci))
fig_tsd <- "reproduce/data/tissue"
## 원 스크립트에서 genelist_*2 정의만 뽑아 실행
fig_src <- readLines(file.path(fig_tsd, "final", "time_DEGreport_code.R"))
fig_gl  <- new.env()
for (ln in grep("^genelist_(liver|brain|kidney)2 <- ", fig_src, value = TRUE)) eval(parse(text = ln), envir = fig_gl)

fig_tissue <- list(
  liver  = list(f = "Liver_final/counts_liver_for_time.csv",
                cond = c(rep('bbb',7), rep('aaa',6), rep('ccc',7), rep('aaa',6), rep('bbb',4)),
                gl = "genelist_liver2"),
  brain  = list(f = "Brain_final/bam/counts_brain_time.csv",
                cond = c(rep('bbb',8), rep('aaa',6), rep('ccc',8), rep('aaa',6), rep('bbb',4)),
                gl = "genelist_brain2"),
  kidney = list(f = "public/kidney/bam/bam/counts_kidney_time.csv",
                cond = c(rep('bbb',8), rep('aaa',6), rep('ccc',7), rep('aaa',5), rep('bbb',4)),
                gl = "genelist_kidney2"))
for (k in names(fig_tissue)) {
  m <- read.csv(file.path(fig_tsd, fig_tissue[[k]]$f)); rownames(m) <- m$Geneid; m <- m[,-1]
  cat(sprintf("%-7s counts %d x %d | condition %d | gene list %d | 매칭 %d\n", k, nrow(m), ncol(m),
      length(fig_tissue[[k]]$cond), length(get(fig_tissue[[k]]$gl, fig_gl)),
      sum(get(fig_tissue[[k]]$gl, fig_gl) %in% rownames(m))))
  fig_tissue[[k]]$mat <- m
}
# ---- DESeq2 + rlog, liver ----
fig_build <- function(k, transform = "rlog") {
  ti <- fig_tissue[[k]]
  cd <- data.frame(row.names = colnames(ti$mat), condition = ti$cond)
  d  <- DESeq(DESeqDataSetFromMatrix(ti$mat, cd, ~ condition), quiet = TRUE)
  rl <- if (transform == "rlog") rlogTransformation(d) else vst(d, blind = TRUE)
  list(dds = d, mat = assay(rl)[get(ti$gl, fig_gl), , drop = FALSE])
}
fig_t_liv <- system.time(fig_liv <- fig_build("liver", "rlog"))
# ---- brain/kidney; degPatterns normalised z-scores ----
fig_brn <- fig_build("brain",  "rlog")
fig_kid <- fig_build("kidney", "rlog")
## 원본과 동일하게 degPatterns 로 normalized 값을 얻되, cluster 는 하나로 합쳐 통째로 그린다
fig_norm <- function(obj) {
  cl <- degPatterns(obj$mat, metadata = colData(obj$dds), time = "condition",
                    minc = 0, plot = FALSE, cutoff = 0.7)
  n <- cl$normalized
  cat("  clusters:", paste(names(table(n$cluster)), table(n$cluster), sep="=", collapse=" "), "\n")
  n
}
cat("liver\n");  fig_n_liv <- fig_norm(fig_liv)
cat("brain\n");  fig_n_brn <- fig_norm(fig_brn)
cat("kidney\n"); fig_n_kid <- fig_norm(fig_kid)
# ---- row plot function ----
fig_row <- function(n, tissue, show_x = FALSE, show_y = FALSE) {
  n$condition <- factor(as.character(n$condition), levels = c("aaa","bbb","ccc"),
                        labels = c("Young","Mid","Old"))
  n$tis <- tissue
  n$top <- sprintf("All %d genes", length(unique(n$genes)))
  ggplot(n, aes(condition, value)) +
    geom_line(aes(group = genes, colour = condition), alpha = .06, linewidth = .3) +
    geom_jitter(aes(colour = condition), width = .13, size = .55, alpha = .45) +
    geom_boxplot(aes(colour = condition), width = .5, fill = NA,
                 outlier.shape = NA, linewidth = .55) +
    facet_grid(tis ~ top) +
    scale_color_nejm(name = NULL, labels = c("1-3","6-18","24-27")) +
    labs(x = if (show_x) "Age(month)" else NULL,
         y = if (show_y) "Z-score of gene abundance" else NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "none",
          axis.text.x = if (show_x) element_text() else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(),
          strip.background = element_rect(fill = "grey90", colour = "black", linewidth = .4))
}
# ---- shared legend ----
suppressMessages(library(cowplot))
fig_r1 <- fig_row(fig_n_liv, "Liver")
fig_r2 <- fig_row(fig_n_kid, "Kidney", show_y = TRUE)
fig_r3 <- fig_row(fig_n_brn, "Brain", show_x = TRUE)
fig_leg <- get_legend(fig_r1 + theme(legend.position = "bottom") +
                      guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 2))))
# ---- outer y-axis label ----
fig_ylab  <- ggdraw() + draw_label("Z-score of gene abundance", angle = 90, size = 9)
# ---- stacked figure (strip text transparent) and save ----
## element_blank() 는 회색 박스까지 없앰 -> 글자만 투명하게 해서 박스는 자리로 남김
fig_nostrip <- function(p) p + theme(strip.text = element_text(colour = NA))
fig_r1 <- fig_nostrip(fig_row(fig_n_liv, "Liver"))
fig_r2 <- fig_nostrip(fig_row(fig_n_kid, "Kidney"))
fig_r3 <- fig_nostrip(fig_row(fig_n_brn, "Brain", show_x = TRUE))
fig_stack <- plot_grid(fig_r1, fig_r2, fig_r3, ncol = 1, align = "v", axis = "lr",
                       rel_heights = c(1, 1, 1.18))
fig_body  <- plot_grid(fig_ylab, fig_stack, nrow = 1, rel_widths = c(0.055, 1))
fig_final <- plot_grid(fig_body, fig_leg, ncol = 1, rel_heights = c(1, 0.08))
ggsave("reproduce/out/degreport_3tissue_stacked_allgenes.pdf",
       fig_final, width = 3.0, height = 5.4)
