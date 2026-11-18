# ---------------------------------------------------------------------------
# P-BP specific chromatin: GSEA NES bar
#   source("Fig1I-J_NES_bars.R")
#   p <- nes_bar("ATAC");  p                     # 화면 확인
#   nes_bar("ATAC", save = TRUE)                 # figures/ 에 png+pdf+tsv 저장
# term 선정(want_*)만 사람이 정하고 NES/p 값은 fgsea 결과에서 그대로 가져온다.
# ---------------------------------------------------------------------------
suppressMessages({library(tidyverse); library(forcats)})

ROOT <- "chromatin"
RANK <- file.path(ROOT, "rank"); FIG <- file.path(ROOT, "figures")

## ---- 선정한 term (fgsea 결과의 term 문자열과 정확히 일치해야 함) ----------
WANT <- list(
  ATAC = tribble(
    ~cat,                      ~term,                                     ~label,
    "Cell cycle / translation", "mitotic cell cycle checkpoint signaling", "Mitotic checkpoint signaling",
    "Cell cycle / translation", "cytoplasmic translation",                 "Cytoplasmic translation",
    "Cell cycle / translation", "mitotic nuclear division",                "Mitotic nuclear division",
    "Cell cycle / translation", "spindle organization",                    "Spindle organization",
    "Genome maintenance",       "DNA-templated DNA replication maintenance of fidelity","DNA replication fidelity",
    "Genome maintenance",       "regulation of telomere maintenance",      "Telomere maintenance",
    "Genome maintenance",       "base-excision repair",                    "Base-excision repair",
    "Genome maintenance",       "DNA damage checkpoint signaling",         "DNA damage checkpoint"),
  H3K27me3 = tribble(
    ~cat,                       ~term,                                       ~label,
    "Development",              "embryonic organ morphogenesis",             "Embryonic morphogenesis",
    "Development",              "cell fate commitment",                      "Cell fate commitment",
    "Development",              "mesenchymal cell differentiation",          "Mesenchymal cell differentiation",
    "Development",              "stem cell differentiation",                 "Stem cell differentiation",
    "SASP/Senescence",          "regulation of response to oxidative stress", "Regulation of oxidative stress",
    "SASP/Senescence",          "regulation of cellular senescence",          "Regulation of senescence",
    "SASP/Senescence",          "interleukin-6-mediated signaling pathway",   "IL-6 signaling pathway",
    "SASP/Senescence",          "response to interleukin-1",                  "Response to IL-1"))

PAL <- list(ATAC     = c("Cell cycle / translation" = "#C79A2E", "Genome maintenance" = "#9B4A62"),
            H3K27me3 = c("Development" = "#7B5EA7", "SASP/Senescence" = "#4E958A"))

# promoter 랭킹의 fgsea_<mark>.rds 를 쓴다.
# fgsea_<mark>_filt.rds 는 다른 랭킹(값이 전혀 다름)이라 여기서 쓰지 않는다.
RDS <- c(ATAC = "fgsea_atac.rds", H3K27ac = "fgsea_h3k27ac.rds",
         H3K27me3 = "fgsea_h3k27me3.rds")

nes_data <- function(mark) {
  fg <- readRDS(file.path(RANK, RDS[[mark]])) |> as_tibble()
  w  <- WANT[[mark]]
  d  <- w |> left_join(fg, by = "term")
  if (any(is.na(d$NES))) stop(sprintf("fgsea 결과에 없는 term: %s",
                                      paste(d$term[is.na(d$NES)], collapse = "; ")))
  d |> mutate(cat = fct_rev(factor(cat, levels = unique(w$cat)))) |>
    arrange(cat, NES) |> mutate(label = fct_inorder(label))
}

nes_bar <- function(mark, save = FALSE, width = 4.2, height = 3.6, xmax = NULL) {
  d   <- nes_data(mark)
  pal <- PAL[[mark]]
  if (is.null(xmax)) xmax <- ceiling(max(d$NES) * 10) / 10 + 0.1
  p <- ggplot(d, aes(NES, label, fill = cat)) +
    geom_col(width = .70) +
    scale_fill_manual(values = pal, name = NULL, breaks = names(pal)) +
    scale_x_continuous(limits = c(0, xmax), breaks = 0:floor(xmax), expand = c(0, 0)) +
    labs(x = "NES", y = NULL,
         title = sprintf("GO Terms of\nP-BP specific %s", mark)) +
    theme_classic(base_size = 12) +
    theme(axis.line  = element_line(colour = "black", linewidth = .5),
          axis.ticks = element_line(colour = "black", linewidth = .5),
          axis.ticks.length = unit(3.5, "pt"),
          axis.text  = element_text(size = 11, colour = "black"),
          axis.title.x = element_text(size = 12, margin = margin(t = 5)),
          plot.title = element_text(size = 12.5, hjust = .5, lineheight = 1.15),
          legend.position = "bottom", legend.text = element_text(size = 10.5),
          legend.key.size = unit(4, "mm"), legend.margin = margin(t = -2),
          plot.margin = margin(6, 14, 4, 6))
  if (save) {
    dir.create(FIG, showWarnings = FALSE)
    ggsave(file.path(FIG, sprintf("%s_4FBPspecific_NES.pdf", mark)), p, width = width, height = height)
    ggsave(file.path(FIG, sprintf("%s_4FBPspecific_NES.png", mark)), p, width = width, height = height, dpi = 300)
    write_tsv(d |> select(cat, term, label, NES, pval, padj, size),
              file.path(FIG, sprintf("%s_4FBPspecific_NES.tsv", mark)))
  }
  p
}
