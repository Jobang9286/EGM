# ------------------------------------------------------------------
# Module enrichment figure: GO (BP/CC/MF) + KEGG dotplot
#   source("plot_module_enrichment.R")
#   p <- make_fig("M1"); p          # 화면 출력만 (저장은 save = TRUE)
#
# term 선정 = GO_SEL / KEGG_SEL 목록 그대로 사용.
# 새 모듈 후보를 뽑을 때는 candidates("M3") 사용.
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(stringr); library(patchwork); library(ggh4x)
})

## Callers may redefine TDIR / PDIR after sourcing.
TDIR <- "EGM_script/1_module_heatmap_trajectory_script/tables"
PDIR <- "EGM_script/1_module_heatmap_trajectory_script/GO_results/plots"
dir.create(PDIR, showWarnings = FALSE, recursive = TRUE)

MODULE_LAB <- c(M1 = "Aged / ECM senescence",
                M3 = "Aged inflammatory / IFN",
                M5 = "Cholesterol & lipid biosynthesis",
                M6 = "Proteostasis & vesicle trafficking",
                M9 = "Cytoskeleton / actin",
                M2 = "Young identity (mid)",
                M11 = "Cell cycle / regeneration")
MODULE_COL <- c(M1 = "#B23A32", M3 = "#E08214", M5 = "#1B9E8A", M6 = "#4CB98A",
                M9 = "#3B7DBF", M2 = "#6A51A3", M11 = "#B8A0D8")
ONT_COL <- c(BP = "#4C9F4C", CC = "#E58F2E", MF = "#4A78B5")

## ---- 확정된 term 목록 ---------------------------------------------
GO_SEL <- list(
  M1 = list(
    BP = c("extracellular matrix organization",
           "tissue remodeling",
           "collagen fibril organization",
           "response to hypoxia",
           "collagen metabolic process",
           "reactive oxygen species metabolic process",
           "response to oxidative stress",
           "leukocyte migration",
           "intrinsic apoptotic signaling pathway",
           "transforming growth factor beta receptor superfamily signaling pathway",
           "cellular senescence"),
    CC = c("extracellular matrix",
           "fibrillar collagen trimer",
           "collagen trimer",
           "lysosomal membrane",
           "secondary lysosome",
           "lysosomal lumen",
           "autophagosome",
           "oxidoreductase complex",
           "focal adhesion",
           "secretory granule membrane"),
    MF = c("extracellular matrix structural constituent",
           "growth factor binding",
           "collagen binding",
           "integrin binding",
           "cytokine receptor binding",
           "cytokine binding",
           "protease binding",
           "peptidase regulator activity",
           "chemokine activity",
           "peptidase inhibitor activity")),
  M3 = list(
    BP = c("response to virus",
           "response to interferon-beta",
           "interferon-mediated signaling pathway",
           "positive regulation of innate immune response",
           "negative regulation of viral genome replication",
           "canonical NF-kappaB signal transduction",
           "response to interferon-alpha",
           "negative regulation of cytokine production"),
    CC = c("extracellular matrix",
           "basement membrane",
           "early endosome",
           "Golgi stack",
           "postsynaptic density",
           "Golgi apparatus subcompartment",
           "beta-catenin-TCF complex",
           "sarcolemma"),
    MF = c("transmembrane receptor protein kinase activity",
           "extracellular matrix structural constituent",
           "growth factor binding",
           "cell adhesion molecule binding",
           "fibronectin binding",
           "integrin binding",
           "double-stranded RNA binding",
           "cytokine receptor activity")),
  M5 = list(
    BP = c("cholesterol biosynthetic process",
           "cholesterol metabolic process",
           "steroid biosynthetic process",
           "isoprenoid biosynthetic process",
           "Golgi vesicle transport",
           "aerobic electron transport chain",
           "phospholipid biosynthetic process",
           "fatty acid metabolic process"),
    CC = c("transport vesicle",
           "respiratory chain complex",
           "peroxisome",
           "cytoplasmic side of endoplasmic reticulum membrane",
           "smooth endoplasmic reticulum"),
    MF = c("electron transfer activity",
           "active monoatomic ion transmembrane transporter activity",
           "oxidoreduction-driven active transmembrane transporter activity",
           "nitric-oxide synthase regulator activity",
           "GO:0016717")),                      # oxidoreductase, paired donors (이름이 너무 길어 ID로 지정)
  M6 = list(
    BP = c("lipid transport",
           "vesicle organization",
           "establishment of protein localization to membrane",
           "endosomal transport",
           "vesicle fusion",
           "macroautophagy",
           "protein folding"),
    CC = c("clathrin-coated vesicle",
           "early endosome",
           "vacuolar membrane",
           "endocytic vesicle",
           "lysosomal membrane"),
    MF = c("catalytic activity, acting on a glycoprotein",
           "G protein activity",
           "phosphatidylinositol-4-phosphate binding",
           "GTP binding",
           "phosphatidylinositol binding")),
  M9 = list(
    BP = c("actin filament organization",
           "actin polymerization or depolymerization",
           "regulation of cellular component size",
           "regulation of actin filament-based process",
           "positive regulation of supramolecular fiber organization",
           "positive regulation of protein polymerization",
           "filopodium assembly"),
    CC = c("ruffle",
           "cell-substrate junction",
           "actin filament",
           "cell leading edge",
           "actin filament bundle",
           "lamellipodium"),
    MF = c("actin binding",
           "actin filament binding",
           "structural constituent of cytoskeleton",
           "GTPase binding",
           "protein folding chaperone")),
  M2 = list(
    BP = c("ribosome biogenesis",
           "rRNA metabolic process",
           "DNA replication",
           "double-strand break repair",
           "telomere maintenance",
           "telomere maintenance via telomerase",
           "regulation of telomere maintenance"),
    CC = c("spliceosomal complex",
           "preribosome",
           "nuclear chromosome",
           "small-subunit processome",
           "telomerase holoenzyme complex",
           "chromosome, telomeric region"),
    MF = c("mRNA binding",
           "snoRNA binding",
           "single-stranded DNA binding",
           "histone binding",
           "telomerase RNA binding",
           "telomeric DNA binding")),
  M11 = list(
    BP = c("chromosome segregation",
           "nuclear division",
           "mitotic cell cycle phase transition",
           "mitotic spindle assembly checkpoint signaling",
           "spindle organization",
           "positive regulation of cell cycle process",
           "DNA replication"),
    CC = c("condensed chromosome",
           "kinetochore",
           "spindle",
           "mitotic spindle",
           "spindle pole",
           "midbody"),
    MF = c("microtubule binding",
           "microtubule motor activity",
           "catalytic activity, acting on DNA",
           "ATP hydrolysis activity",
           "DNA endonuclease activity",
           "plus-end-directed microtubule motor activity"))
)

KEGG_SEL <- list(
  M1 = c("Lysosome biogenesis", "HIF-1 signaling pathway", "TNF signaling pathway",
         "FoxO signaling pathway", "PI3K-Akt signaling pathway", "Efferocytosis",
         "IL-17 signaling pathway", "Integrin signaling", "Toll-like receptor signaling pathway",
         "Cytokine-cytokine receptor interaction", "NOD-like receptor signaling pathway",
         "Chemokine signaling pathway", "JAK-STAT signaling pathway", "Autophagy - animal",
         "Peroxisome", "p53 signaling pathway", "ECM-receptor interaction",
         "Cellular senescence", "TGF-beta signaling pathway", "NF-kappa B signaling pathway"),
  M3 = c("Integrin signaling", "PI3K-Akt signaling pathway", "Focal adhesion",
         "Rap1 signaling pathway", "ECM-receptor interaction", "Cadherin signaling",
         "JAK-STAT signaling pathway", "Wnt signaling pathway",
         "RIG-I-like receptor signaling pathway", "Toll-like receptor signaling pathway",
         "Human papillomavirus infection", "Herpes simplex virus 1 infection",
         "Epstein-Barr virus infection", "Influenza A"),
  M5 = c("Steroid biosynthesis", "Terpenoid backbone biosynthesis", "Oxidative phosphorylation",
         "Thermogenesis", "Endocytosis", "Fatty acid metabolism", "Sphingolipid signaling pathway",
         "Biosynthesis of unsaturated fatty acids", "Insulin signaling pathway",
         "mTOR signaling pathway"),
  M6 = c("Endocytosis", "Oxidative phosphorylation", "Lysosome biogenesis",
         "Protein processing in endoplasmic reticulum", "Autophagy - animal",
         "Thermogenesis", "Synaptic vesicle cycle", "Non-alcoholic fatty liver disease"),
  M9 = c("Regulation of actin cytoskeleton", "Tight junction",
         "Fc gamma R-mediated phagosome formation", "Cadherin signaling", "IgSF CAM signaling",
         "Motor proteins", "Cytoskeleton in muscle cells", "Focal adhesion",
         "Salmonella infection"),  # 액틴 조작 경로라 Human Diseases 예외
  M2 = c("Nucleocytoplasmic transport", "Ribosome biogenesis in eukaryotes", "DNA replication",
         "Cell cycle", "mRNA surveillance pathway", "Nucleotide metabolism", "Proteasome",
         "Nucleotide excision repair", "One carbon pool by folate", "Mismatch repair"),
  M11 = c("Cell cycle", "Oocyte meiosis", "Motor proteins",
          "Progesterone-mediated oocyte maturation", "Fanconi anemia pathway")
)

## 모듈별로 KEGG 에 덧붙일 Reactome pathway (KEGG 가 얇은 모듈용)
REACT_SEL <- list(
  M1 = c("Extracellular matrix organization",
         "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)",
         "Collagen formation",
         "Collagen biosynthesis and modifying enzymes",
         "Degradation of the extracellular matrix"),
  M3 = c("Modulation of host responses by IFN-stimulated genes",
         "ISG15 antiviral mechanism",
         "Collagen degradation",
         "Extracellular matrix organization",
         "Signaling by PDGF"),
  M5 = c("Cholesterol biosynthesis",
         "Lanosterol biosynthesis",
         "Cholesterol biosynthesis via desmosterol (Bloch pathway)",
         "eNOS activation"),
  M9 = c("Regulation of actin dynamics for phagocytic cup formation",
         "RHO GTPase Effectors",
         "RHO GTPases Activate WASPs and WAVEs",
         "Smooth Muscle Contraction"),
  M2 = c("Processing of Capped Intron-Containing Pre-mRNA",
         "mRNA Splicing - Major Pathway",
         "snRNP Assembly",
         "Extension of Telomeres",
         "Telomere C-strand (Lagging Strand) Synthesis",
         "Telomere Maintenance"),
  M11 = c("M Phase", "Mitotic Prometaphase", "Cell Cycle Checkpoints",
          "Amplification of signal from the kinetochores", "RHO GTPase Effectors",
          "Activation of the pre-replicative complex",
          "Polymerase switching on the C-strand of the telomere"))

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

## ---- 후보 탐색용 (다음 모듈에서 사용) -----------------------------
SEN_KEY <- paste(
  "senescen","aging","ageing","DNA damage","telomer","cell cycle arrest","p53","apopto","autophag",
  "lysosom","oxidative stress","reactive oxygen","superoxide","oxidoreduct",
  "inflamm","cytokine","interleukin","chemokine","chemotax","NF-kappa","tumor necrosis","TGF",
  "leukocyte","macrophage","myeloid","immune","complement","acute-phase",
  "extracellular matrix","collagen","matrix","remodel","wound","fibro","angiogen","hypoxia",
  "integrin","adhesion","peptidase","protease","proteolysis","growth factor","secretory granule",
  sep = "|")

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

candidates <- function(mod, key = SEN_KEY, n = 18, cut = 0.8, drop_cat = "Human Diseases") {
  for (o in c("BP", "CC", "MF")) {
    f <- file.path(TDIR, sprintf("%s_GO_%s.csv", mod, o))
    if (!file.exists(f)) next
    cat(sprintf("\n########## %s %s 후보 ##########\n", mod, o))
    print(read.csv(f) |> arrange(qvalue) |>
            filter(grepl(key, Description, ignore.case = TRUE)) |> drop_redundant(cut) |> head(n) |>
            transmute(Description, Count, q = signif(qvalue, 2)))
  }
  f <- file.path(TDIR, sprintf("%s_KEGG.csv", mod))
  if (file.exists(f)) {
    cat(sprintf("\n########## %s KEGG 후보 ##########\n", mod))
    print(read.csv(f) |> filter(!category %in% drop_cat) |> arrange(qvalue) |> head(28) |>
            transmute(cat = substr(category, 1, 18), Description, Count,
                      FE = round(FoldEnrichment, 2), q = signif(qvalue, 2)))
  }
}

## ---- 데이터 뽑기 ---------------------------------------------------
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

go_sel <- function(mod) {
  bind_rows(lapply(c("BP", "CC", "MF"), function(o) {
    want <- GO_SEL[[mod]][[o]]; if (is.null(want)) return(NULL)
    d <- read.csv(file.path(TDIR, sprintf("%s_GO_%s.csv", mod, o)))
    miss <- setdiff(tolower(want[!grepl("^GO:", want)]), tolower(d$Description))
    miss <- c(miss, setdiff(want[grepl("^GO:", want)], d$ID))
    if (length(miss)) warning(sprintf("%s %s: 표에 없는 term -> %s", mod, o, paste(miss, collapse = "; ")))
    d |> filter(tolower(Description) %in% tolower(want) | ID %in% want) |> arrange(qvalue) |>
      transmute(ONT = o, ID, Description = cap1(relabel(Description)), Count, qvalue,
                logq = -log10(qvalue))
  }))
}

kegg_sel <- function(mod) {
  d <- read.csv(file.path(TDIR, sprintf("%s_KEGG.csv", mod)))
  want <- KEGG_SEL[[mod]]
  miss <- setdiff(want, d$Description)
  if (length(miss)) warning(sprintf("%s KEGG: 표에 없는 pathway -> %s", mod, paste(miss, collapse = "; ")))
  out <- d |> filter(Description %in% want) |> mutate(Description = relabel(Description)) |> arrange(qvalue)
  rw <- REACT_SEL[[mod]]
  if (!is.null(rw)) {                                   # Reactome 수혈
    f <- file.path(TDIR, sprintf("%s_Reactome.csv", mod))
    r <- read.csv(f)
    rmiss <- setdiff(rw, r$Description)
    if (length(rmiss)) warning(sprintf("%s Reactome: 표에 없는 pathway -> %s", mod, paste(rmiss, collapse = "; ")))
    r <- r |> filter(Description %in% rw) |> mutate(Description = relabel(Description),
                                                    category = "Reactome", subcategory = NA_character_) |>
      arrange(qvalue)
    out <- bind_rows(out |> mutate(category = "KEGG"), r |> select(any_of(names(out))))
  }
  out
}


## ---- 축약(half) 버전 ----------------------------------------------
GO_KEEP  <- list(M1 = c("Cellular senescence"),        # half 에서도 반드시 유지 (축약된 라벨 기준)
                 M3 = c("Double-stranded RNA binding", "Cytokine receptor activity", "Early endosome"),
                 M5 = c("Aerobic electron transport chain", "Smooth ER"),
                 M6 = c("Macroautophagy", "Protein folding"),
                 M9 = c("Filopodium assembly", "Lamellipodium"),
                 M2 = c("Telomere maintenance", "Telomere maintenance (telomerase)",
                        "Reg. of telomere maintenance",
                        "Chromosome, telomeric region", "Telomerase RNA binding",
                        "Telomeric DNA binding"),
                 M11 = c("DNA replication", "Pos. reg. of cell cycle process", "Midbody",
                         "Plus-end microtubule motor", "DNA endonuclease activity"))
KEGG_KEEP <- list(M1 = c("Cellular senescence"),
                  M3 = c("RLR signaling", "TLR signaling", "JAK-STAT signaling"),
                  M5 = c("Unsaturated FA biosynthesis", "Insulin signaling", "mTOR signaling"),
                  M6 = c("NAFLD", "Synaptic vesicle cycle"),
                  M9 = c("Focal adhesion", "Salmonella infection"),
                  M2 = c("Mismatch repair", "One carbon pool by folate"),
                  M11 = c("Fanconi anemia pathway"))

half_go <- function(g, mod, n = 5) {
  keep <- GO_KEEP[[mod]]
  g |> group_by(ONT) |> arrange(qvalue, .by_group = TRUE) |>
    filter(row_number() <= n | Description %in% keep) |> ungroup()
}
half_kegg <- function(k, mod, n = 10) {
  keep <- KEGG_KEEP[[mod]]
  k <- k |> arrange(qvalue)
  if (any(k$category %in% "Reactome")) {       # KEGG+Reactome 혼합 모듈
    kg <- k |> filter(category == "KEGG")
    nk <- sum(kg$Description %in% keep)
    kg <- bind_rows(kg |> filter(Description %in% keep),
                    kg |> filter(!Description %in% keep) |> head(max(0, n - nk))) |> arrange(qvalue)
    return(bind_rows(kg, k |> filter(category == "Reactome")))
  }
  bind_rows(k |> filter(Description %in% keep),
            k |> filter(!Description %in% keep) |> head(n - sum(k$Description %in% keep))) |>
    arrange(qvalue)
}

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
                      x_breaks = 2, n_size = NULL, group = c("full", "coarse", "none")) {
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
  k <- k |> mutate(category = factor(category, levels = lv))
  if (group == "coarse" && !any(k$category %in% c("KEGG", "Reactome")))
    k <- k |> mutate(category = factor(unname(COARSE[as.character(category)]),
                                       levels = c("Signaling", "Cellular / metabolic")))
  k <- k |> arrange(category, FoldEnrichment) |>
    mutate(lab = str_wrap(Description, wrap), lab = factor(lab, levels = unique(lab)))
  if (group == "none") k <- k |> arrange(FoldEnrichment) |>
    mutate(lab = factor(str_wrap(Description, wrap), levels = str_wrap(Description, wrap)))
  ggplot(k, aes(FoldEnrichment, lab, size = Count, colour = pvalue)) +
    geom_point() +
    (if (group == "none" || nlevels(droplevels(k$category)) < 2) NULL else
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
          legend.key.height = unit(0.5, "cm"))
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

make_fig <- function(mod, half = FALSE, free_x = TRUE, width = 15, height = 8.5,
                     n_go = 5, n_kegg = 10, panel_w = NULL, txt = 8.5, title = FALSE,
                     dot_range = c(1.2, 4.6), kegg_group = "coarse", x_breaks = 2,
                     save = FALSE, suffix = if (half) "_half" else "") {
  g <- go_sel(mod); k <- kegg_sel(mod)
  if (half) { g <- half_go(g, mod, n_go); k <- half_kegg(k, mod, n_kegg) }
  # panel_w 를 주면 패널 폭 고정(라벨 공간 확보), NULL 이면 상대 폭 배분
  p <- plot_go(g, free_x = free_x, panel_w = panel_w, txt = txt) |
    plot_kegg(k, panel_w = panel_w, txt = txt, dot_range = dot_range,
                    x_breaks = x_breaks, group = kegg_group)
  if (is.null(panel_w)) p <- p + plot_layout(widths = c(1, 1.05))
  if (title) p <- p + plot_annotation(title = sprintf("%s — %s", mod, MODULE_LAB[[mod]]),
                                      theme = theme(plot.title = element_text(
                                        face = "bold", size = 14, colour = MODULE_COL[[mod]])))
  if (save) {
    ggsave(file.path(PDIR, sprintf("%s_GO_KEGG%s.png", mod, suffix)), p, width = width, height = height,
           dpi = 300, bg = "white")
    suppressWarnings(ggsave(file.path(PDIR, sprintf("%s_GO_KEGG%s.pdf", mod, suffix)), p,
                            width = width, height = height, bg = "white"))
    write.csv(g, file.path(PDIR, sprintf("%s_GO_used%s.csv", mod, suffix)), row.names = FALSE)
    write.csv(k[, c("category","subcategory","ID","Description","Count","FoldEnrichment",
                    "pvalue","p.adjust","qvalue")],
              file.path(PDIR, sprintf("%s_KEGG_used%s.csv", mod, suffix)), row.names = FALSE)
  }
  p
}
