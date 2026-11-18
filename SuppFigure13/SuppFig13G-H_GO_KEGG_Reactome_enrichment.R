# ---------------------------------------------------------------------------
# Human EGM vs Control: volcano + up/down GO (BP/CC/MF), KEGG, Reactome enrichment,
# drawn with the EGM template figure functions (SuppFig13G-H_lib_enrichment_plot.R).
#
# Enrichment recipe = run_module_enrichment.R (clusterProfiler enrichGO SYMBOL BP/CC/MF,
# enrichKEGG hsa via ENTREZ, ReactomePA enrichPathway; BH, p<0.05, q<0.05). UNIVERSE below:
#   "orgdb"  = template default (whole annotation as background)  <- used for the figures
#   "tested" = the 14,548 genes tested by DESeq2 (monkey functional-study convention; stricter:
#              the DOWN side then keeps only 20 BP terms and no apoptosis/p53 terms).
#              Those tables are kept in tables/strict_universe/.
# DEG definition = this project's module definition: padj < 0.05 & |apeglm lfc| > 0.5.
#
# Term selection: PICK_* lists below. NULL = automatic top-N by q-value after
# drop_redundant() (Jaccard >= 0.8) and, for KEGG, without Human Diseases / Drug Development.
# Fill PICK_* with exact Description strings to curate (the script stops if one is missing).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA)
  library(dplyr); library(ggplot2); library(ggrepel); library(patchwork); library(stringr); library(ggh4x)
})
find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  cands <- c(if (length(f)) dirname(dirname(normalizePath(f[1]))), getwd(), dirname(getwd()))
  hit <- cands[dir.exists(file.path(cands, "data"))]
  if (!length(hit)) stop("cannot locate the data/ directory; run from the project root")
  hit[1]
}
ROOT <- find_root(); D <- file.path(ROOT, "data"); FIG <- file.path(ROOT, "figures"); TAB <- file.path(ROOT, "tables")
dir.create(TAB, showWarnings = FALSE)
source("SuppFig13G-H_lib_enrichment_plot.R")
select <- dplyr::select; filter <- dplyr::filter   # AnnotationDbi masks these when clusterProfiler loads first

PADJ <- 0.05; LFC <- 0.5
UNIVERSE <- "orgdb"   # "orgdb" (template default) | "tested"
# Curated term lists: UP = rejuvenation-positive programmes, DOWN = inflammaging / stress
# programmes that EGM switches off. Set any entry to NULL to fall back to automatic top-N.
PICK_GO <- list(
  up   = c("mitotic nuclear division", "chromosome segregation", "DNA replication", "mitotic cell cycle checkpoint signaling", "wound healing",
           "cholesterol biosynthetic process", "isoprenoid biosynthetic process"),
           # lipid BP = the two clean mevalonate/cholesterol terms only. "positive regulation of lipid localization"
           # (C3, EDN1, IL1B, PLA2G4A, PTGES, NFKB1 ...), "unsaturated fatty acid biosynthetic process" (EDN1, IL1B, PTGS2, PLA2G4A ...,
           # i.e. prostaglandin/eicosanoid genes) and "LDL particle clearance" (PCSK9, MYLIP, CD36, IL19 mixed in) are inflammation-driven
           # or mixed, so they are NOT used. Lipid-storage / fat-cell terms are likewise inflammation-driven and left out.
  down = c("regulation of apoptotic signaling pathway", "intrinsic apoptotic signaling pathway",
           "intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress",
           "extrinsic apoptotic signaling pathway", "signal transduction by p53 class mediator",
           "response to endoplasmic reticulum stress", "endoplasmic reticulum unfolded protein response"))
           # no virus / antiviral terms; apoptosis, p53 and ER-stress terms instead (template default universe)
PICK_KEGG <- list(
  up   = c("Cell cycle", "DNA replication", "Homologous recombination", "Steroid biosynthesis"),
  down = c("p53 signaling pathway", "Ferroptosis", "Integrated stress response (ISR) signaling pathway"))   # ISR category missing in clusterProfiler's KEGG table -> set below
PICK_REACT <- list(
  up   = c("Cell Cycle Checkpoints", "Telomere Maintenance", "Cholesterol biosynthesis", "Regulation of cholesterol biosynthesis by SREBP (SREBF)"),
           # "Cell Cycle, Mitotic" (>500 genes) is excluded by maxGSSize under the default universe -> Cell Cycle Checkpoints
  down = c("Unfolded Protein Response (UPR)", "PERK regulates gene expression",
           "ATF4 activates genes in response to endoplasmic reticulum stress", "KEAP1-NFE2L2 pathway", "Interferon alpha/beta signaling"))
N_GO <- 7; N_KEGG <- 8; N_REACT <- 4
# one-line GO labels for the panel (template LAB_ORGAN convention: shorten only what would wrap at 34 chars)
LAB_PANEL <- c("Mitotic cell cycle checkpoint signaling" = "Mitotic checkpoint signaling",
               "Intrinsic apoptotic signaling pathway in response to ER stress" = "ER stress-induced apoptosis",
               "Signal transduction by p53 class mediator" = "p53-mediated signaling",
               "Unsaturated fatty acid biosynthetic process" = "Unsaturated FA biosynthesis",
               "Low-density lipoprotein particle clearance" = "LDL particle clearance",
               "Regulation of cholesterol biosynthesis by SREBP (SREBF)" = "Reg. of cholesterol biosynthesis by SREBP")
one_line <- function(x) { i <- match(x, names(LAB_PANEL)); ifelse(is.na(i), x, LAB_PANEL[i]) }
KEGG_STRIP <- list(up = "KEGG", down = "KEGG")   # plain "KEGG" / "Reactome" strips on both panels
H_SCALE <- 0.85   # both panels 15 % shorter than the template auto height (0.33 x rows + 0.9 in)
TXT <- 10; TXT_BASE <- 11.5; LAB_GO <- 2.2; LAB_KEGG <- 3.0; FIG_W <- 8.8   # larger text; label columns widened to match
STRIP_FILL <- FALSE  # TRUE = coloured category strips + legend (organ-panel look), kept available in lib   # "category" = KEGG class strips (organ template); "KEGG" = one plain strip (module template)

## ---- DE table -------------------------------------------------------------
sig <- read.csv(file.path(D, "EGM_vs_Control_signature.csv")) |> filter(named, !is.na(padj))
uni <- unique(sig$id)
up  <- sig$id[sig$padj < PADJ & sig$lfc >  LFC]
dn  <- sig$id[sig$padj < PADJ & sig$lfc < -LFC]
cat(sprintf("tested %d named genes | up %d | down %d  (padj<%.2f, |lfc|>%.1f)\n", length(uni), length(up), length(dn), PADJ, LFC))

## ---- volcano --------------------------------------------------------------
CAP <- 100
v <- sig |> mutate(y = -log10(pmax(padj, 1e-300)), capped = y > CAP, y = pmin(y, CAP),
                   cls = case_when(padj < PADJ & lfc >  LFC ~ "Up in EGM",
                                   padj < PADJ & lfc < -LFC ~ "Down in EGM", TRUE ~ "n.s."),
                   cls = factor(cls, c("Up in EGM", "Down in EGM", "n.s.")))
lab_n <- 12
lab <- bind_rows(v |> filter(cls == "Up in EGM")   |> arrange(desc(y), desc(abs(lfc))) |> head(lab_n),
                 v |> filter(cls == "Down in EGM") |> arrange(desc(y), desc(abs(lfc))) |> head(lab_n))
cnt <- table(v$cls)
volc <- ggplot(v, aes(lfc, y)) +
  geom_hline(yintercept = -log10(PADJ), linetype = 2, linewidth = .3, colour = "grey40") +
  geom_vline(xintercept = c(-LFC, LFC), linetype = 2, linewidth = .3, colour = "grey40") +
  geom_point(data = v |> filter(cls == "n.s."), colour = "grey82", size = .7, alpha = .7) +
  geom_point(data = v |> filter(cls != "n.s."), aes(colour = cls, shape = capped), size = 1.1, alpha = .8) +
  geom_text_repel(data = lab, aes(label = id, colour = cls), size = 2.6, fontface = "italic",
                  max.overlaps = 40, box.padding = .25, segment.size = .25, show.legend = FALSE) +
  scale_colour_manual(values = c("Up in EGM" = "#B2182B", "Down in EGM" = "#2166AC", "n.s." = "grey82"),
                      labels = c(sprintf("Up in EGM (n = %d)", cnt[["Up in EGM"]]),
                                 sprintf("Down in EGM (n = %d)", cnt[["Down in EGM"]]), "n.s."), name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  scale_y_continuous(breaks = seq(0, CAP, 25), labels = function(b) ifelse(b == CAP, paste0("≥", CAP), b)) +
  labs(x = expression(log[2]~"fold change (EGM / Control, apeglm)"), y = expression(-log[10]~"adjusted"~italic(P)),
       title = "EGM vs Control, primary human fibroblasts",
       subtitle = sprintf("DESeq2, n = 3 vs 3; %d genes tested; triangles = padj below 1e-%d", length(uni), CAP)) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(.02, .02), legend.justification = c(0, 0), legend.background = element_blank(),
        legend.key.size = unit(.4, "cm"), plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"))
ggsave(file.path(FIG, "volcano_EGM_vs_Control.png"), volc, width = 5.6, height = 5.0, dpi = 300, bg = "white")
ggsave(file.path(FIG, "volcano_EGM_vs_Control.pdf"), volc, width = 5.6, height = 5.0, bg = "white")
write.csv(v |> filter(cls != "n.s.") |> dplyr::select(id, lfc, lfc_raw, baseMean, pvalue, padj, cls) |> arrange(cls, desc(abs(lfc))),
          file.path(TAB, "DEG_EGM_vs_Control_padj0.05_lfc0.5.csv"), row.names = FALSE)
cat("volcano written\n")

## ---- enrichment -----------------------------------------------------------
P_CUT <- 0.05; Q_CUT <- 0.05
n_terms <- function(res) if (is.null(res)) 0L else nrow(as.data.frame(res))
uni_eg <- suppressMessages(suppressWarnings(bitr(uni, "SYMBOL", "ENTREZID", org.Hs.eg.db)))$ENTREZID
U_sym <- if (UNIVERSE == "tested") uni else NULL; U_eg <- if (UNIVERSE == "tested") uni_eg else NULL
run_enrich <- function(genes, tag) {
  out <- list(); summ <- list()
  for (ont in c("BP", "CC", "MF")) {
    ego <- tryCatch(enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = ont, universe = U_sym,
                             pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, minGSSize = 10, maxGSSize = 500),
                    error = function(e) { message(tag, " GO_", ont, ": ", conditionMessage(e)); NULL })
    d <- as.data.frame(ego); out[[paste0("GO_", ont)]] <- d
    if (nrow(d)) write.csv(d, file.path(TAB, sprintf("EGM_%s_GO_%s.csv", tag, ont)), row.names = FALSE)
    summ[[paste0("GO_", ont)]] <- nrow(d); cat(sprintf("  %s GO_%s: %d terms\n", tag, ont, nrow(d)))
  }
  eg <- suppressMessages(suppressWarnings(bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)))$ENTREZID
  cat(sprintf("  %s ENTREZ mapped %d/%d\n", tag, length(eg), length(genes)))
  ek <- tryCatch(enrichKEGG(eg, organism = "hsa", keyType = "kegg", universe = U_eg,
                            pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, minGSSize = 5, maxGSSize = 600),
                 error = function(e) { message(tag, " KEGG: ", conditionMessage(e)); NULL })
  if (n_terms(ek) > 0) ek <- setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  d <- as.data.frame(ek); out$KEGG <- d
  if (nrow(d)) write.csv(d, file.path(TAB, sprintf("EGM_%s_KEGG.csv", tag)), row.names = FALSE)
  summ$KEGG <- nrow(d); cat(sprintf("  %s KEGG: %d pathways\n", tag, nrow(d)))
  er <- tryCatch(enrichPathway(eg, organism = "human", universe = U_eg, pAdjustMethod = "BH",
                               pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, readable = TRUE),
                 error = function(e) { message(tag, " Reactome: ", conditionMessage(e)); NULL })
  d <- as.data.frame(er); out$Reactome <- d
  if (nrow(d)) write.csv(d, file.path(TAB, sprintf("EGM_%s_Reactome.csv", tag)), row.names = FALSE)
  summ$Reactome <- nrow(d); cat(sprintf("  %s Reactome: %d pathways\n", tag, nrow(d)))
  attr(out, "summary") <- unlist(summ); out
}
cat("enrichment: up\n");   E_up <- run_enrich(up, "up")
cat("enrichment: down\n"); E_dn <- run_enrich(dn, "down")
write.csv(data.frame(direction = rep(c("up", "down"), each = 5), db = names(attr(E_up, "summary")),
                     n_terms = c(attr(E_up, "summary"), attr(E_dn, "summary")), n_genes = rep(c(length(up), length(dn)), each = 5)),
          file.path(TAB, "enrichment_summary_counts.csv"), row.names = FALSE)

## ---- term selection (template rules) ---------------------------------------
sel_go <- function(d, want, n) {
  if (!nrow(d)) return(d); d$Description <- str_squish(d$Description)
  if (is.null(want)) d <- d |> arrange(qvalue) |> drop_redundant(0.8) |> head(n)
  else { miss <- setdiff(tolower(want), tolower(d$Description)); if (length(miss)) stop("GO: not in table -> ", paste(miss, collapse = "; "))
         d <- d |> filter(tolower(Description) %in% tolower(want)) }
  d |> arrange(qvalue) |> transmute(ONT = "BP", ID, Description = one_line(cap1(relabel(Description))), Count, qvalue, logq = -log10(qvalue))
}
sel_kegg <- function(d, want, n) {
  if (!nrow(d)) return(d); d$Description <- str_squish(d$Description)
  d <- d |> mutate(category = ifelse(is.na(category) | category == "", "Environmental Information Processing", category),
                   subcategory = ifelse(is.na(subcategory) | subcategory == "", "Signal transduction", subcategory))
  if (is.null(want)) d <- d |> filter(!category %in% c("Human Diseases", "Drug Development")) |> arrange(qvalue) |> drop_redundant(0.8) |> head(n)
  else { miss <- setdiff(tolower(want), tolower(d$Description)); if (length(miss)) stop("KEGG: not in table -> ", paste(miss, collapse = "; "))
         d <- d |> filter(tolower(Description) %in% tolower(want)) }
  d |> mutate(Description = relabel(Description)) |> arrange(qvalue)
}
sel_react <- function(d, want, n) {
  if (!nrow(d)) return(d); d$Description <- str_squish(d$Description)
  if (is.null(want)) d <- d |> arrange(qvalue) |> drop_redundant(0.8) |> head(n)
  else { miss <- setdiff(tolower(want), tolower(d$Description)); if (length(miss)) stop("Reactome: not in table -> ", paste(miss, collapse = "; "))
         d <- d |> filter(tolower(Description) %in% tolower(want)) }
  d |> mutate(Description = one_line(relabel(Description)), category = "Reactome", subcategory = NA_character_) |> arrange(qvalue)
}
make_dir_fig <- function(E, tag, pick_go, pick_kegg, pick_react) {
  g <- sel_go(E$GO_BP, pick_go, N_GO)
  k <- sel_kegg(E$KEGG, pick_kegg, N_KEGG)
  r <- sel_react(E$Reactome, pick_react, N_REACT)
  if (nrow(k) && KEGG_STRIP[[tag]] == "KEGG") k <- k |> mutate(category = "KEGG")
  if (nrow(r)) k <- bind_rows(if (nrow(k)) k |> mutate(category = as.character(category)) else NULL, r |> dplyr::select(any_of(names(k)), everything()))
  if (!nrow(g) || !nrow(k)) { message(tag, ": nothing to plot (GO ", nrow(g), ", KEGG/Reactome ", nrow(k), ")"); return(invisible(NULL)) }
  cat(sprintf("\n[%s] GO BP shown:\n", tag)); print(g |> transmute(Description, Count, q = signif(qvalue, 2)), row.names = FALSE)
  cat(sprintf("[%s] KEGG/Reactome shown:\n", tag)); print(k |> transmute(category, Description, Count, FE = round(FoldEnrichment, 2), q = signif(qvalue, 2)), row.names = FALSE)
  enrich_fig(g, k, height = (0.33 * max(nrow(g), nrow(k)) + 0.9) * H_SCALE, width = FIG_W,
             txt = TXT, txt_base = TXT_BASE, lab_go = LAB_GO, lab_kegg = LAB_KEGG,
             stem = sprintf("enrichment_EGM_%s_GO_KEGG", tag), outdir = FIG, strip_fill = STRIP_FILL)
}
p_up <- make_dir_fig(E_up, "up",   PICK_GO$up,   PICK_KEGG$up,   PICK_REACT$up)
p_dn <- make_dir_fig(E_dn, "down", PICK_GO$down, PICK_KEGG$down, PICK_REACT$down)
## Panel titles are not drawn.
saveRDS(list(E_up = E_up, E_dn = E_dn, up = up, dn = dn, universe = uni), file.path(D, "enrichment_EGM_vs_Control.rds"))
cat("\ndone\n")
