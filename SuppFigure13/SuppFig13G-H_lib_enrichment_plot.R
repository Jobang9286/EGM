# ---------------------------------------------------------------------------
# SuppFig13G-H_lib_enrichment_plot.R  --  figure functions of the EGM enrichment templates
#   plot_module_enrichment.R
#     (ONT_COL, LAB_SHORT, shorten, drop_redundant, cap1, relabel, plot_go, plot_kegg, blank_text)
#   plot_organ_enrichment.R
#     (label_width_in; organ_fig -> enrich_fig below, same geometry: fixed panel widths,
#      fixed label-column widths, gap, auto height)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(stringr); library(patchwork); library(ggh4x)
})
ONT_COL <- c(BP = "#4C9F4C", CC = "#E58F2E", MF = "#4A78B5")
## coloured category strips + legend for plot_kegg(strip_fill = TRUE),
## Set2 colours: Metabolism teal, Cellular Processes orange,
## Organismal Systems pink, Environmental Information lilac-blue, Genetic Information green.
KCAT_COL <- c("Metabolism" = "#66C2A5", "Genetic Information Processing" = "#A6D854",
              "Environmental Information Processing" = "#8DA0CB", "Cellular Processes" = "#FC8D62",
              "Organismal Systems" = "#E78AC3", "Human Diseases" = "#FFD92F", "Drug Development" = "#E5C494",
              "KEGG" = "#B3B3B3", "Reactome" = "#B3B3B3")
KCAT_LAB <- c("Metabolism" = "Metabolism", "Genetic Information Processing" = "Genetic information",
              "Environmental Information Processing" = "Environmental information", "Cellular Processes" = "Cellular process",
              "Organismal Systems" = "Organismal systems", "Human Diseases" = "Human diseases", "Drug Development" = "Drug development",
              "KEGG" = "KEGG", "Reactome" = "Reactome")

## ---- LAB_SHORT / shorten -------------------------------------------------
## 짧은 라벨 (그림에 표시할 이름만 축약)
LAB_SHORT <- c(
  "transforming growth factor beta receptor superfamily signaling pathway" = "TGF-beta receptor signaling",
  "reactive oxygen species metabolic process" = "ROS metabolic process",
  "intrinsic apoptotic signaling pathway"     = "Intrinsic apoptotic signaling",
  "extracellular matrix structural constituent" = "ECM structural constituent",
  "extracellular matrix organization"         = "ECM organization",
  "Cytokine-cytokine receptor interaction"    = "Cytokine-cytokine receptor int.",
  "Viral protein interaction with cytokine and cytokine receptor" = "Viral protein-cytokine int.",
  "transmembrane receptor protein kinase activity" = "Transmembrane receptor kinase",
  "Golgi apparatus subcompartment"            = "Golgi subcompartment",
  "active monoatomic ion transmembrane transporter activity" = "Active ion transporter",
  "oxidoreduction-driven active transmembrane transporter activity" = "Oxidoreduction-driven transporter",
  "nitric-oxide synthase regulator activity" = "NOS regulator activity",
  "NADH dehydrogenase (ubiquinone) activity" = "NADH dehydrogenase activity",
  "Biosynthesis of unsaturated fatty acids"  = "Unsaturated FA biosynthesis",
  "actin polymerization or depolymerization" = "Actin (de)polymerization",
  "positive regulation of supramolecular fiber organization" = "Pos. reg. of fiber organization",
  "regulation of actin filament-based process" = "Reg. of actin filament process",
  "regulation of cellular component size"    = "Reg. of cell component size",
  "positive regulation of protein polymerization" = "Pos. reg. of protein polymerization",
  "structural constituent of cytoskeleton"   = "Cytoskeletal structural constituent",
  "Fc gamma R-mediated phagosome formation"  = "FcgR phagosome formation",
  "Regulation of actin cytoskeleton"         = "Actin cytoskeleton regulation",
  "ribonucleoprotein complex binding"        = "RNP complex binding",
  "telomere maintenance via telomerase"      = "Telomere maintenance (telomerase)",
  "G-rich strand telomeric DNA binding"      = "G-rich telomeric DNA binding",
  "mitotic spindle assembly checkpoint signaling" = "Mitotic spindle checkpoint",
  "mitotic sister chromatid segregation"     = "Mitotic chromatid segregation",
  "mitotic cell cycle phase transition"      = "Mitotic cell cycle transition",
  "Progesterone-mediated oocyte maturation"  = "Oocyte maturation (progesterone)",
  "plus-end-directed microtubule motor activity" = "Plus-end microtubule motor",
  "Amplification of signal from the kinetochores" = "Kinetochore signal amplification",
  "Polymerase switching on the C-strand of the telomere" = "Telomere C-strand pol switching",
  "Activation of the pre-replicative complex" = "Pre-replicative complex activation",
  "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)" = "IGF transport/uptake by IGFBPs",
  "Collagen biosynthesis and modifying enzymes" = "Collagen biosynthesis",
  "Degradation of the extracellular matrix"  = "ECM degradation",
  "Modulation of host responses by IFN-stimulated genes" = "Host response by ISGs",
  "Cholesterol biosynthesis via desmosterol (Bloch pathway)" = "Cholesterol biosynth. (Bloch)",
  "Regulation of actin dynamics for phagocytic cup formation" = "Actin dynamics, phagocytic cup",
  "RHO GTPases Activate WASPs and WAVEs"     = "RHO GTPases activate WASP/WAVE",
  "Processing of Capped Intron-Containing Pre-mRNA" = "Capped pre-mRNA processing",
  "Telomere C-strand (Lagging Strand) Synthesis" = "Telomere C-strand synthesis")

# KEGG/GO 공통 축약 규칙: 뒤에 붙는 " signaling pathway" -> " signaling" 등
shorten <- function(x) {           # base:: 명시 (전역환경에 sub 등 동명 객체가 있어도 안전)
  x <- base::sub("^oxidoreductase activity, acting on paired donors.*", "Oxidoreductase, paired donors", x)
  x <- base::gsub("endoplasmic reticulum", "ER", x)
  x <- base::sub("^establishment of protein localization", "Protein localization", x)
  x <- base::sub("^[Cc]atalytic activity, acting on (a )?", "Catalytic activity on ", x)
  x <- base::gsub("phosphatidylinositol", "PtdIns", x)
  x <- base::sub("Non-alcoholic fatty liver disease", "NAFLD", x)
  x <- base::gsub("interferon", "IFN", x)
  x <- base::sub("NF-kappaB", "NF-kB", x); x <- base::sub("NF-kappa B", "NF-kB", x)
  x <- base::sub("Human papillomavirus", "HPV", x)
  x <- base::sub("Herpes simplex virus 1", "HSV-1", x)
  x <- base::sub("Epstein-Barr virus", "EBV", x)
  x <- base::sub("innate immune response", "innate immunity", x)
  x <- base::sub("adaptive immune response", "adaptive immunity", x)
  x <- base::sub("viral genome replication", "viral replication", x)
  x <- base::sub("^positive regulation of ", "Pos. reg. of ", x)
  x <- base::sub("^negative regulation of ", "Neg. reg. of ", x)
  x <- base::sub("^regulation of ", "Reg. of ", x)
  x <- base::sub(" signal transduction$", " signaling", x)
  x <- base::sub(" signaling pathway$", " signaling", x)
  x <- base::sub("^Toll-like receptor", "TLR", x)
  x <- base::sub("^NOD-like receptor", "NLR", x)
  x <- base::sub("^RIG-I-like receptor", "RLR", x)
  x <- base::sub(" - animal$", "", x)
  unname(x)
}

## ---- drop_redundant ------------------------------------------------------
drop_redundant <- function(df, cut = 0.8) {          # 유전자 집합 Jaccard >= cut 이면 중복으로 제거
  keep <- integer(0); sets <- strsplit(df$geneID, "/")
  for (i in seq_len(nrow(df))) {
    dup <- FALSE
    for (j in keep) {
      inter <- length(intersect(sets[[i]], sets[[j]]))
      uni   <- length(union(sets[[i]], sets[[j]]))
      if (uni > 0 && inter / uni >= cut) { dup <- TRUE; break }
    }
    if (!dup) keep <- c(keep, i)
  }
  df[keep, , drop = FALSE]
}

## ---- cap1 / relabel ------------------------------------------------------
# 첫 글자만 대문자로. 단 mRNA/rRNA/snoRNA/mTOR 처럼 이미 대문자를 품은 단어는 그대로 둠
cap1 <- function(x) {
  first <- base::sub(" .*$", "", x)
  keep  <- base::grepl("[A-Z]", first)
  ifelse(keep, x, paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))))
}
relabel <- function(x) {                       # LAB_SHORT 은 대소문자 무시하고 매칭
  i <- match(tolower(x), tolower(names(LAB_SHORT)))
  shorten(ifelse(is.na(i), x, LAB_SHORT[i]))
}

## ---- plot_go / plot_kegg / blank_text ------------------------------------
## ---- 패널 ---------------------------------------------------------
plot_go <- function(g, wrap = 34, free_x = TRUE, panel_w = NULL, txt = 8.5) {
  g <- g |> arrange(ONT, logq) |>
    mutate(lab = str_wrap(Description, wrap), lab = factor(lab, levels = unique(lab)))
  ggplot(g, aes(logq, lab, fill = ONT)) +
    geom_col(width = 0.72) +
    # free_x = TRUE -> BP/CC/MF 블록마다 x축 따로 (facet_wrap2 라야 x가 자유로워짐)
    (if (free_x)
       facet_wrap2(vars(ONT), ncol = 1, scales = "free", strip.position = "right",
                   strip = strip_themed(
                     background_y = elem_list_rect(fill = unname(ONT_COL[sort(unique(g$ONT))])),
                     text_y = elem_list_text(colour = "white", face = "bold", angle = -90)))
     else
       facet_grid2(ONT ~ ., scales = "free_y", space = "free_y",
                   strip = strip_themed(
                     background_y = elem_list_rect(fill = unname(ONT_COL[sort(unique(g$ONT))])),
                     text_y = elem_list_text(colour = "white", face = "bold", angle = -90)))) +
    (if (!is.null(panel_w)) force_panelsizes(cols = unit(panel_w, "in")) else NULL) +
    scale_fill_manual(values = ONT_COL, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                       breaks = scales::pretty_breaks(n = 4)) +
    labs(x = expression(-log[10](q)), y = NULL, title = "GO enrichment") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text.y = element_text(size = txt, colour = "black"),
          strip.text.y = element_text(size = 9))
}

# group: "full" = KEGG category 그대로, "coarse" = 2개로 통합, "none" = 패널 분할 없음
plot_kegg <- function(k, wrap = 32, panel_w = NULL, txt = 8.5, dot_range = c(1.2, 4.6),
                      x_breaks = 2, n_size = NULL, group = c("full", "coarse", "none"), strip_fill = FALSE) {
  group <- match.arg(group)
  if (group != "none" && nrow(k) <= 3 && !any(k$category %in% "Reactome")) group <- "none"
  COARSE <- c("Environmental Information Processing" = "Signaling",
              "Organismal Systems"                   = "Signaling",
              "Cellular Processes"                   = "Cellular / metabolic",
              "Metabolism"                           = "Cellular / metabolic",
              "Genetic Information Processing"       = "Cellular / metabolic",
              "Human Diseases"                       = "Cellular / metabolic",
              "Drug Development"                     = "Cellular / metabolic")
  cat_lv <- c("Metabolism", "Genetic Information Processing", "Environmental Information Processing",
              "Cellular Processes", "Organismal Systems", "Human Diseases", "Drug Development")
  # 스트립이 좁아 긴 카테고리명이 잘리므로 표시용으로 줄인다
  CAT_SHORT <- c("Environmental Information Processing" = "Environmental",
                 "Genetic Information Processing"       = "Genetic info.",
                 "Organismal Systems"                   = "Organismal",
                 "Cellular Processes"                   = "Cellular",
                 "Human Diseases"                       = "Disease")
  lv <- intersect(cat_lv, unique(k$category)); if (!length(lv)) lv <- unique(k$category)
  if ("Reactome" %in% k$category && !"Reactome" %in% lv) lv <- c(lv, "Reactome")   # keep Reactome as its own strip below the KEGG categories
  k <- k |> mutate(category = factor(category, levels = lv))
  if (group == "coarse" && !any(k$category %in% c("KEGG", "Reactome")))
    k <- k |> mutate(category = factor(unname(COARSE[as.character(category)]),
                                       levels = c("Signaling", "Cellular / metabolic")))
  k <- k |> arrange(category, FoldEnrichment) |>
    mutate(lab = str_wrap(Description, wrap), lab = factor(lab, levels = unique(lab)))
  if (group == "none") k <- k |> arrange(FoldEnrichment) |>
    mutate(lab = factor(str_wrap(Description, wrap), levels = str_wrap(Description, wrap)))
  cats <- levels(droplevels(k$category))
  ggplot(k, aes(FoldEnrichment, lab, size = Count, colour = pvalue)) +
    geom_point() +
    (if (strip_fill)   # invisible points carrying the category fill -> legend swatches
       list(geom_point(aes(fill = category), shape = 22, alpha = 0, size = 0, show.legend = TRUE),
            scale_fill_manual(values = KCAT_COL, labels = KCAT_LAB, name = "KEGG category", drop = TRUE),
            guides(fill = guide_legend(order = 0, override.aes = list(alpha = 1, size = 4, shape = 22, colour = NA))))
     else NULL) +
    (if (group == "none" || nlevels(droplevels(k$category)) < 2) NULL else if (strip_fill)
       facet_grid2(category ~ ., scales = "free_y", space = "free_y",
                   strip = strip_themed(background_y = elem_list_rect(fill = unname(KCAT_COL[cats]), colour = NA),
                                        text_y = elem_list_text(colour = NA)))
     else
       facet_grid2(category ~ ., scales = "free_y", space = "free_y",
                   labeller = labeller(category = function(x)
                     str_wrap(ifelse(x %in% names(CAT_SHORT), CAT_SHORT[x], x), 10)))) +
    (if (!is.null(panel_w)) force_panelsizes(cols = unit(panel_w, "in")) else NULL) +
    scale_colour_gradientn(colours = c("#D7191C", "#FDAE61", "#FFFFBF", "#ABD9E9", "#2C7BB6"),
                           labels = scales::label_scientific(digits = 2), name = "pvalue") +
    scale_size_continuous(range = dot_range, name = "Number",
                          breaks = if (is.null(n_size)) waiver() else scales::breaks_pretty(n = n_size)) +
    guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) +  # pvalue 위, Number 아래
    scale_x_continuous(expand = expansion(mult = c(0.22, 0.22)),
                       breaks = scales::pretty_breaks(n = x_breaks)) +
    coord_cartesian(clip = "off") +                                  # 양끝 점이 패널 경계에 잘리지 않게
    labs(x = "Enrichment Score", y = NULL,
         title = if (any(k$category %in% "Reactome")) "KEGG / Reactome Enrichment" else "KEGG Enrichment") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text.y = element_text(size = txt, colour = "black"),
          strip.text.y = element_text(size = 7.5, angle = -90),
          strip.background = element_rect(fill = "grey92", colour = "grey60"),
          legend.key.height = unit(0.5, "cm")) +
    (if (strip_fill) theme(strip.text.y = element_text(size = 7.5, angle = -90, margin = margin(l = 2.5, r = 2.5)),
                           strip.background = element_rect(colour = NA)) else NULL)
}

# 글자 없는 버전 (일러스트에서 텍스트를 직접 얹을 때). 패널/막대/점 위치는 그대로 둔다.
blank_text <- function(p) {
  # 글자를 투명색으로 두면 라벨 버전과 패널 위치·크기가 정확히 같아 겹쳐 쓰기 좋다
  inv <- element_text(colour = NA)
  p & theme(axis.text.x = inv, axis.text.y = inv, axis.title.x = inv, axis.title.y = inv,
            plot.title = inv, strip.text.x = inv, strip.text.y = inv,
            legend.text = inv, legend.title = inv,
            axis.ticks = element_line(colour = "black"))
}

## ---- label_width_in ------------------------------------------------------
label_width_in <- function(p) {
  g  <- ggplotGrob(p)
  ax <- which(grepl("^axis-l", g$layout$name))
  if (!length(ax)) return(0)
  sum(grid::convertWidth(g$widths[unique(g$layout$l[ax])], "in", valueOnly = TRUE))
}

## ---- enrich_fig: organ_fig geometry, but takes the two data frames directly ----
# g: GO table with ONT, Description, Count, qvalue, logq ; k: KEGG(+Reactome) table with
# category, Description, Count, FoldEnrichment, pvalue. Same defaults as organ_fig.
enrich_fig <- function(g, k, width = 8.2, height = NULL, txt = 8.5, panel_titles = FALSE,
                       dot_range = c(1.2, 4.6), x_breaks = 2, n_size = 3, kegg_group = "full",
                       panel_go = 1.0, panel_kegg = 0.85, gap = 0.25,
                       lab_go = 1.9, lab_kegg = 2.6, wrap_go = 34, wrap_kegg = 48,
                       title = NULL, stem = NULL, outdir = NULL, strip_fill = FALSE, txt_base = NULL) {
  if (is.null(height)) height <- 0.33 * max(nrow(g), nrow(k)) + 0.9
  p1 <- plot_go(g, panel_w = panel_go, txt = txt, wrap = wrap_go)
  p2 <- plot_kegg(k, dot_range = dot_range, x_breaks = x_breaks, n_size = n_size,
                  group = kegg_group, panel_w = panel_kegg, txt = txt, wrap = wrap_kegg, strip_fill = strip_fill)
  pad1 <- max(0, lab_go   - label_width_in(p1)) * 72
  pad2 <- max(0, lab_kegg - label_width_in(p2)) * 72
  p <- (p1 + theme(plot.margin = margin(6, gap * 36, 4, 6 + pad1))) |
    (p2 + theme(plot.margin = margin(6, 6, 4, gap * 36 + pad2)))
  if (!panel_titles) p <- p & theme(plot.title = element_blank())
  if (!is.null(txt_base))   # scale every text element that is not set explicitly (axis titles, ticks, legend)
    p <- p & theme(text = element_text(size = txt_base), axis.text.x = element_text(size = txt_base * 0.8),
                   legend.text = element_text(size = txt_base * 0.8), legend.title = element_text(size = txt_base * 0.9),
                   strip.text.y = element_text(size = txt_base * 0.8))
  if (!is.null(title))
    p <- p + plot_annotation(title = title, theme = theme(plot.title = element_text(face = "bold", size = 13)))
  if (!is.null(stem) && !is.null(outdir)) {
    pb <- blank_text(p)
    ggsave(file.path(outdir, paste0(stem, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
    suppressWarnings(ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white"))
    ggsave(file.path(outdir, paste0(stem, "_blank.png")), pb, width = width, height = height, dpi = 300, bg = "white")
    suppressWarnings(ggsave(file.path(outdir, paste0(stem, "_blank.pdf")), pb, width = width, height = height, bg = "white"))
    write.csv(g, file.path(outdir, paste0(stem, "_GO_used.csv")), row.names = FALSE)
    write.csv(k |> dplyr::select(any_of(c("category", "subcategory", "ID", "Description", "Count", "FoldEnrichment",
                                   "pvalue", "p.adjust", "qvalue"))),
              file.path(outdir, paste0(stem, "_KEGG_used.csv")), row.names = FALSE)
  }
  p
}
