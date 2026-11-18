REPO_DIR <- ".."                     # EGM_github repository root, relative to Figure5/

# ---------------------------------------------------------------------------
# Monkey 장기별 enrichment 그림
#   그림 함수(plot_go / plot_kegg / 라벨 축약)는 EGM 모듈 스크립트를 그대로 재사용한다.
#   source("Fig5M-N_03_plot_organ_enrichment.R")
#   p <- organ_fig("Kidney"); p
#   organ_fig("Kidney", save = TRUE)
# term/pathway 는 사용자가 고른 PICK_GO / PICK_KEGG 목록을 그대로 쓴다 (자동 선택 없음).
# ---------------------------------------------------------------------------
source(file.path(REPO_DIR, "shared/enrichment_plotting/plot_module_enrichment.R"))
suppressPackageStartupMessages(library(readxl))

MDIR <- "EGM_script/3_monkey_GO_tables"
ODIR <- file.path(MDIR, "plots")          # PDIR(모듈용)과 겹치지 않게 별도 이름
dir.create(ODIR, showWarnings = FALSE)

DIRECTION <- c(Kidney = "up", Liver = "up", Brain = "up",
               Spleen_up = "up", Spleen_down = "down")
organ_of <- function(key) base::sub("_(up|down)$", "", key)   # 그림 키 -> 장기명

## ---- 장기별로 고른 term / pathway -----------------------------------------
PICK_GO <- list(
  Kidney = c("generation of precursor metabolites and energy",
             "proton transmembrane transport",
             "proton motive force-driven mitochondrial ATP synthesis",
             "ATP metabolic process",
             "ATP biosynthetic process",
             "lipid catabolic process"),
  Liver = c("generation of precursor metabolites and energy",
            "fatty acid catabolic process",
            "mitochondrial ATP synthesis coupled electron transport",
            "fatty acid metabolic process",
            "lipid catabolic process",
            "fatty acid beta-oxidation"),
  Brain = c("modulation of chemical synaptic transmission",
            "synapse organization",
            "learning or memory",
            "regulation of neuron projection development",
            "cognition",
            "dendrite development",
            "neuron projection guidance"),
  Spleen_up = c("blood vessel morphogenesis",
                "external encapsulating structure organization",
                "ameboidal-type cell migration",
                "cell-matrix adhesion",
                "response to wounding",
                "endothelial cell migration",
                "epithelial cell differentiation"),
  Spleen_down = c("activation of immune response",
                  "T cell activation",
                  "regulation of response to biotic stimulus",
                  "immune response-regulating signaling pathway",
                  "leukocyte mediated immunity",
                  "regulation of innate immune response",
                  "regulation of immune effector process",
                  "immune response-activating signaling pathway",
                  "cell killing",
                  "leukocyte migration",
                  "cytokine-mediated signaling pathway")
)

PICK_KEGG <- list(
  Kidney = c("Oxidative phosphorylation",
             "Citrate cycle (TCA cycle)",
             "Glycolysis / Gluconeogenesis",
             "Fatty acid degradation",
             "Proximal tubule bicarbonate reclamation",
             "Aldosterone-regulated sodium reabsorption"),
  Liver = c("Oxidative phosphorylation",
            "Fatty acid degradation",
            "Peroxisome",
            "PPAR signaling pathway",
            "Citrate cycle (TCA cycle)",
            "Mitophagy - animal",
            "Cholesterol metabolism"),
  Brain = c("Neuroactive ligand-receptor interaction",
            "Axon guidance",
            "Glutamatergic synapse",
            "GABAergic synapse",
            "Cholinergic synapse",
            "Growth hormone synthesis, secretion and action",
            "Long-term potentiation",
            "Dopaminergic synapse",
            "Calcium signaling pathway",
            "cAMP signaling pathway"),
  Spleen_up = c("Integrin signaling",
                "ECM-receptor interaction",
                "Adherens junction",
                "Signaling pathways regulating pluripotency of stem cells",
                "Longevity regulating pathway - multiple species",
                "Fatty acid biosynthesis",
                "FoxO signaling pathway",
                "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
                "Leukocyte transendothelial migration"),
  Spleen_down = c("NOD-like receptor signaling pathway",
                  "Th17 cell differentiation",
                  "Chemokine signaling pathway",
                  "Toll-like receptor signaling pathway",
                  "JAK-STAT signaling pathway",
                  "NF-kappa B signaling pathway")
)

## 장기 그림은 표에 있는 term 문구를 그대로 쓴다. 다만 아래 목록에 있는 것만
## 한 줄에 들어가도록 줄여 쓴다 (Spleen down 처럼 긴 면역 term 이 두 줄로 접히는 경우).
LAB_ORGAN <- c(
  "regulation of response to biotic stimulus"    = "Response to biotic stimulus",
  "immune response-regulating signaling pathway" = "Immune response-regulating",
  "regulation of innate immune response"         = "Reg. of innate immunity",
  "regulation of immune effector process"        = "Reg. of immune effector",
  "immune response-activating signaling pathway" = "Immune response-activating",
  "cytokine-mediated signaling pathway"          = "Cytokine-mediated signaling",
  "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate" =
    "Glycosaminoglycan biosynthesis")

raw_label <- function(x) {
  i <- match(tolower(x), tolower(names(LAB_ORGAN)))
  cap1(ifelse(is.na(i), x, LAB_ORGAN[i]))
}

## ---- 데이터 (FULL_*.xlsx) --------------------------------------------------
read_tab <- function(organ, db) {          # db: "GO" | "KEGG"
  f <- file.path(MDIR, sprintf("FULL_%s%s_%s.xlsx", db, DIRECTION[[organ]], organ_of(organ)))
  if (!file.exists(f)) stop("파일 없음: ", f)
  suppressMessages(read_excel(f, sheet = 1))
}

organ_go <- function(organ) {              # 이 표는 GO BP 전용
  d <- read_tab(organ, "GO"); want <- PICK_GO[[organ]]
  miss <- setdiff(want, d$Description)
  if (length(miss)) stop(sprintf("%s GO: 표에 없는 term -> %s", organ, paste(miss, collapse = "; ")))
  d |> filter(Description %in% want) |>
    transmute(ONT = "BP", ID, Description = raw_label(Description), Count, qvalue,
              logq = -log10(qvalue)) |> arrange(qvalue)
}

organ_kegg <- function(organ) {
  d <- read_tab(organ, "KEGG"); want <- PICK_KEGG[[organ]]
  miss <- setdiff(want, d$Description)
  if (length(miss)) stop(sprintf("%s KEGG: 표에 없는 pathway -> %s", organ, paste(miss, collapse = "; ")))
  d |> filter(Description %in% want) |> mutate(Description = raw_label(Description)) |> arrange(qvalue)
}

## ---- 그림 (모듈 템플릿 그대로) ---------------------------------------------
## y축 라벨 칸 폭을 인치로 고정 -> 장기가 달라도 패널 왼쪽 모서리 위치가 같아진다
label_width_in <- function(p) {
  g  <- ggplotGrob(p)
  ax <- which(grepl("^axis-l", g$layout$name))
  if (!length(ax)) return(0)
  sum(grid::convertWidth(g$widths[unique(g$layout$l[ax])], "in", valueOnly = TRUE))
}

# panel_go / panel_kegg: 패널(플롯 영역) 폭을 인치로 고정 -> 장기가 달라도 패널 크기와
# 두 패널 사이 간격(gap)이 같게 유지된다. gap 은 인치.
organ_fig <- function(organ, title = FALSE, width = 8.2, height = NULL, save = FALSE, txt = 8.5, panel_titles = FALSE,
                      dot_range = c(1.2, 4.6), x_breaks = 2, n_size = 3, kegg_group = "full",
                      panel_go = 1.0, panel_kegg = 0.85, gap = 0.25,
                      lab_go = 1.9, lab_kegg = 2.6, wrap_go = 34, wrap_kegg = 48) {
  g <- organ_go(organ); k <- organ_kegg(organ)
  # 항목 수가 달라도 한 줄 높이가 같도록 저장 높이를 자동 계산
  if (is.null(height)) height <- 0.33 * max(nrow(g), nrow(k)) + 0.9
  p1 <- plot_go(g, panel_w = panel_go, txt = txt, wrap = wrap_go)
  p2 <- plot_kegg(k, dot_range = dot_range, x_breaks = x_breaks, n_size = n_size,
                  group = kegg_group, panel_w = panel_kegg, txt = txt, wrap = wrap_kegg)
  # 라벨 칸을 lab_go / lab_kegg 인치로 맞춰 두 장기의 패널 위치를 동일하게 만든다
  pad1 <- max(0, lab_go   - label_width_in(p1)) * 72
  pad2 <- max(0, lab_kegg - label_width_in(p2)) * 72
  p <- (p1 + theme(plot.margin = margin(6, gap * 36, 4, 6 + pad1))) |
    (p2 + theme(plot.margin = margin(6, 6, 4, gap * 36 + pad2)))               # 왼쪽 여백 = gap
  if (!panel_titles) p <- p & theme(plot.title = element_blank())   # "GO enrichment"/"KEGG Enrichment" 제거
  if (title)
    p <- p + plot_annotation(title = sprintf("%s — %s in EGM", organ, toupper(DIRECTION[[organ]])),
                             theme = theme(plot.title = element_text(face = "bold", size = 13)))
  if (save) {
    pb <- blank_text(p)
    ggsave(file.path(ODIR, sprintf("%s_GO_KEGG.png", organ)), p, width = width, height = height,
           dpi = 300, bg = "white")
    suppressWarnings(ggsave(file.path(ODIR, sprintf("%s_GO_KEGG.pdf", organ)), p,
                            width = width, height = height, bg = "white"))
    ggsave(file.path(ODIR, sprintf("%s_GO_KEGG_blank.png", organ)), pb, width = width,
           height = height, dpi = 300, bg = "white")
    suppressWarnings(ggsave(file.path(ODIR, sprintf("%s_GO_KEGG_blank.pdf", organ)), pb,
                            width = width, height = height, bg = "white"))
    write.csv(g, file.path(ODIR, sprintf("%s_GO_used.csv", organ)), row.names = FALSE)
    write.csv(k |> select(category, subcategory, ID, Description, Count, FoldEnrichment,
                          pvalue, p.adjust, qvalue),
              file.path(ODIR, sprintf("%s_KEGG_used.csv", organ)), row.names = FALSE)
  }
  p
}
