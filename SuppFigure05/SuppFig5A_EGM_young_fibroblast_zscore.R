## ---- RRHO2 inputs -----------------------------------------------------------
suppressMessages(library(RRHO2))

fig_dir <- "reproduce"
fig_rank <- function(f) {
  x <- read.csv(file.path(fig_dir, "data", f))
  x <- x[!duplicated(x[, 1]), ]
  rownames(x) <- x[, 1]
  x
}

fig_egm  <- fig_rank("egm_vs_no_egm_for_RRHO.csv")
fig_fb   <- fig_rank("fb_young_vs_old_for_RRHO.csv")

## ---- RRHO2 (fibroblast x EGM) -----------------------------------------------
fig_rrho <- function(a, b, la, lb) {
  cg <- intersect(rownames(a), rownames(b))
  RRHO2_initialize(as.data.frame(a[cg, ]), as.data.frame(b[cg, ]),
                   labels = c(la, lb), log10.ind = TRUE)
}

fig_o_fb <- fig_rrho(fig_fb,   fig_egm, "Young fibroblast", "EGM")

## ---- ear fibroblast raw counts ----------------------------------------------
suppressMessages(library(DESeq2))

fig_ef <- "reproduce/data/ear_fib_counts"
fig_ff <- list.files(fig_ef, "_htseq_count\\.txt$", full.names = TRUE)
fig_l  <- lapply(fig_ff, function(p) { x <- read.delim(p, header = FALSE, row.names = 1)
                                       colnames(x) <- sub("^Fib_|_htseq_count\\.txt$", "", basename(p)); x })
stopifnot(all(sapply(fig_l, function(x) identical(rownames(x), rownames(fig_l[[1]])))))
fig_efm <- do.call(cbind, fig_l)
fig_efm <- fig_efm[!grepl("^__", rownames(fig_efm)), ]
# 23,366 genes x 18 samples

## ---- sample table -----------------------------------------------------------
colnames(fig_efm) <- sub("_htseq_count\\.txt$", "", colnames(fig_efm))
fig_cd_ef <- data.frame(age = factor(ifelse(grepl("^29", colnames(fig_efm)), "m29", "m3"),
                                     levels = c("m29", "m3")), row.names = colnames(fig_efm))

## ---- the 65 EGM-Young fibroblast genes --------------------------------------
fig_65 <- fig_o_fb$genelist_uu$gene_list_overlap_uu

## ---- helpers ----------------------------------------------------------------
suppressMessages({library(ggplot2); library(dplyr)})

## 유전자별 z-score -> 그룹별 평균 (점 1개 = 유전자 1개), 선은 같은 유전자를 그룹 간 연결
## (per-gene z-score -> group mean; one point = one gene; lines join the same gene)
fig_zdf <- function(dds, genes, grp, grp_levels, grp_labels) {
  v <- assay(vst(dds, blind = TRUE))
  g <- intersect(genes, rownames(v)); v <- v[g, , drop = FALSE]
  z <- t(scale(t(v)))                                     # 유전자별 z (per-gene z)
  gr <- factor(grp, levels = grp_levels, labels = grp_labels)
  d <- do.call(rbind, lapply(levels(gr), function(k)
        data.frame(gene = g, group = k, z = rowMeans(z[, gr == k, drop = FALSE]))))
  d$group <- factor(d$group, levels = grp_labels)
  attr(d, "n") <- length(g); d
}

fig_zplot <- function(d, title, pal, leg_labels) {
  ggplot(d, aes(group, z)) +
    geom_line(aes(group = gene), colour = "#9ecae1", alpha = .35, linewidth = .3) +
    geom_boxplot(aes(colour = group), width = .45, outlier.shape = NA, fill = NA, linewidth = .6) +
    geom_jitter(aes(colour = group), width = .12, size = 1, alpha = .8) +
    scale_colour_manual(values = pal, labels = leg_labels, name = "Age (month)") +
    labs(title = title, x = "Age", y = "Z-score of gene abundance") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = .5, size = 10),
          panel.grid.minor = element_blank())
}

## ---- DESeq2 object (size factors only) + z-score table ----------------------
fig_dds_ef2 <- estimateSizeFactors(DESeqDataSetFromMatrix(as.matrix(fig_efm), fig_cd_ef, ~ age))
fig_d_fb <- fig_zdf(fig_dds_ef2, fig_65, fig_cd_ef$age, c("m3","m29"), c("3 month","29 month"))

## ---- plot -------------------------------------------------------------------
fig_pal_fb <- c(`3 month` = "#E8735C", `29 month` = "#E8A33C")

fig_p_fb <- fig_zplot(fig_d_fb, "EGM-Young fibroblast correlated factors (n=65)",
                      fig_pal_fb, c("3", "29"))

## ---- final cosmetic settings ------------------------------------------------
fig_p_fb$layers[[3]]$aes_params$alpha <- 0.5   # geom_jitter: 0.8 -> 0.5
fig_p_fb$layers[[1]]$aes_params$colour <- "#E8735C"   # 연결선: 파란색 -> 3 month 색 (lines: blue -> 3-month colour)
fig_p_fb$layers[[1]]$aes_params$alpha     <- 0.2    # 연결선 투명도 (line alpha)
fig_p_fb$layers[[2]]$aes_params$linewidth <- 0.6    # 박스는 원래대로 (box outline back to 0.6)
fig_p_fb <- fig_p_fb + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.1))
ggsave("reproduce/out/zscore_EGM_fibroblast_65.pdf", fig_p_fb, width = 4.0, height = 2.9)
