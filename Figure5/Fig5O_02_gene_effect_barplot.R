## ============================================================================
## Fig5O_02_gene_effect_barplot.R
## Panel D (Fig4g-analog): 7 up / 7 down curated genes per organ (Brain/Liver/
## Kidney/Spleen), red/cyan staircase with adjP stars.
##
## The UP-regulated gene set for each organ is drawn directly from the GO-BP /
## KEGG terms that are displayed in that organ's enrichment figure
## (Fig5M-N_03_plot_organ_enrichment.R, PICK_GO /
## PICK_KEGG lists) -- i.e. the exact same terms shown in <organ>.pdf /
## plots/<organ>_GO_KEGG.pdf. The DOWN-regulated set keeps classic senescence
## markers (TP53, CDKN1A) wherever they are themselves significantly
## down-regulated in that tissue, and fills remaining slots with the largest
## |logFC| significant down genes overall (CDKN2A does not pass expression
## filtering / is not present in the macaque annotation in this dataset, so it
## never qualifies).
##
## Inputs:
##   repro_native_clock_egm/outputs/egm_effect_all.csv   limma EGM-vs-Saline DE
##     (logFC / adj.P.Val / AveExpr / tissue / symbol), from egm_effect.R
##   3_monkey_GO_tables/FULL_GOup_<Organ>.xlsx            enrichGO BP results
##   3_monkey_GO_tables/FULL_KEGGup_<Organ>.xlsx           enrichKEGG results
##   (both provide a `geneID` column: "/"-separated HGNC gene symbols that are
##    members of that term, i.e. the padj<0.05 DESeq2 DEGs used as enrichGO/
##    enrichKEGG input for that organ x "up" direction)
##
## Gene selection:
##   1. candidate pool = union of geneID members across PICK_GO[[organ]] and
##      PICK_KEGG[[organ]] (the terms actually displayed in the organ figure)
##   2. keep genes that are also DE in the independent limma table
##      (adj.P.Val < 0.05, AveExpr >= 2 -- excludes near-zero/"on-off" genes
##      whose huge fold change is a dropout artefact, e.g. LCT in Kidney,
##      CHGB in Kidney), same direction as candidate pool ("up")
##   3. rank by |logFC|, take top 7; LOC-prefixed and common blood-
##      contamination genes (ALB, HBB, HBA1/2) excluded
##   4. DOWN side: same filters but logFC < 0; TP53 (Liver) / CDKN1A (Spleen)
##      force-included when they themselves qualify (adj.P<0.05, down); the
##      remaining slots filled by largest |logFC|
##
## x-axis ticks are set by hand per organ (not auto):
## Brain -3/0/3/6, Liver -3/0/3,
## Kidney -3/-1.5/0/1.5, Spleen -2/0/2.
##
## Outputs (OUT folder):
##   panelD_v5_geneEffect.csv   final gene table (tissue, dir, symbol, logFC, adjP, basis)
##   panelD_v5_geneEffect.png / .pdf   the figure
## ============================================================================
suppressMessages({ library(dplyr); library(ggplot2); library(ggh4x); library(readxl) })

GODIR <- "EGM_script/3_monkey_GO_tables"
RP    <- "monkey/reproduce data/repro_native_clock_egm"
OUT   <- "3_monkey_GO_tables/panelD_gene_effect_v5"

## ---- 1. terms displayed in each organ's enrichment figure (same as
## Fig5M-N_03_plot_organ_enrichment.R PICK_GO / PICK_KEGG; Spleen here
## = the "up" direction figure, matching this panel's "up-regulated" focus) ----
PICK_GO <- list(
  Brain  = c("modulation of chemical synaptic transmission", "synapse organization",
             "learning or memory", "regulation of neuron projection development",
             "cognition", "dendrite development", "neuron projection guidance"),
  Liver  = c("generation of precursor metabolites and energy", "fatty acid catabolic process",
             "mitochondrial ATP synthesis coupled electron transport", "fatty acid metabolic process",
             "lipid catabolic process", "fatty acid beta-oxidation"),
  Kidney = c("generation of precursor metabolites and energy", "proton transmembrane transport",
             "proton motive force-driven mitochondrial ATP synthesis", "ATP metabolic process",
             "ATP biosynthetic process", "lipid catabolic process"),
  Spleen = c("blood vessel morphogenesis", "external encapsulating structure organization",
             "ameboidal-type cell migration", "cell-matrix adhesion", "response to wounding",
             "endothelial cell migration", "epithelial cell differentiation"))

PICK_KEGG <- list(
  Brain  = c("Neuroactive ligand-receptor interaction", "Axon guidance", "Glutamatergic synapse",
             "GABAergic synapse", "Cholinergic synapse", "Growth hormone synthesis, secretion and action",
             "Long-term potentiation", "Dopaminergic synapse", "Calcium signaling pathway",
             "cAMP signaling pathway"),
  Liver  = c("Oxidative phosphorylation", "Fatty acid degradation", "Peroxisome",
             "PPAR signaling pathway", "Citrate cycle (TCA cycle)", "Mitophagy - animal",
             "Cholesterol metabolism"),
  Kidney = c("Oxidative phosphorylation", "Citrate cycle (TCA cycle)", "Glycolysis / Gluconeogenesis",
             "Fatty acid degradation", "Proximal tubule bicarbonate reclamation",
             "Aldosterone-regulated sodium reabsorption"),
  Spleen = c("Integrin signaling", "ECM-receptor interaction", "Adherens junction",
             "Signaling pathways regulating pluripotency of stem cells",
             "Longevity regulating pathway - multiple species", "Fatty acid biosynthesis",
             "FoxO signaling pathway",
             "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
             "Leukocyte transendothelial migration"))

## GO/KEGG file on disk is "Spleen" for both up and down directions (organ
## name only, no _up/_down suffix in the filename)
FILE_ORGAN <- c(Brain = "Brain", Liver = "Liver", Kidney = "Kidney", Spleen = "Spleen")

get_candidates <- function(organ) {
  go   <- read_excel(file.path(GODIR, sprintf("FULL_GOup_%s.xlsx", FILE_ORGAN[[organ]])), sheet = 1)
  kegg <- read_excel(file.path(GODIR, sprintf("FULL_KEGGup_%s.xlsx", FILE_ORGAN[[organ]])), sheet = 1)
  go_genes   <- unlist(strsplit(go$geneID[go$Description %in% PICK_GO[[organ]]], "/"))
  kegg_genes <- unlist(strsplit(kegg$geneID[kegg$Description %in% PICK_KEGG[[organ]]], "/"))
  sort(unique(c(go_genes, kegg_genes)))
}

## ---- 2. limma EGM-vs-Saline DE table (egm_effect.R output) ----
de <- read.csv(file.path(RP, "outputs/egm_effect_all.csv"))
de <- de[!is.na(de$symbol) & de$symbol != "", ]
de_best <- de[order(de$adj.P.Val), ]
de_best <- de_best[!duplicated(de_best[, c("tissue", "symbol")]), ]

EXCLUDE_ARTIFACT <- c("ALB", "HBB", "HBA1", "HBA2")   # common blood-contamination markers

pick_up <- function(tissue, candidates, n = 7, minExpr = 2) {
  sub <- de_best[de_best$tissue == tissue & de_best$symbol %in% candidates &
                   de_best$logFC > 0 & de_best$adj.P.Val < 0.05 & de_best$AveExpr >= minExpr, ]
  sub <- sub[!grepl("^LOC", sub$symbol) & !sub$symbol %in% EXCLUDE_ARTIFACT, ]
  sub <- sub[order(-abs(sub$logFC)), ]
  head(sub, n)
}
pick_down <- function(tissue, force_genes = character(0), n = 7, minExpr = 2, maxAbsFC = 8) {
  sub <- de_best[de_best$tissue == tissue & de_best$logFC < 0 &
                   de_best$adj.P.Val < 0.05 & de_best$AveExpr >= minExpr &
                   abs(de_best$logFC) <= maxAbsFC, ]
  sub <- sub[!grepl("^LOC", sub$symbol) & !sub$symbol %in% EXCLUDE_ARTIFACT, ]
  sub <- sub[order(-abs(sub$logFC)), ]
  forced <- sub[sub$symbol %in% force_genes, ]
  rest   <- sub[!sub$symbol %in% force_genes, ]
  out <- rbind(forced, head(rest, n - nrow(forced)))
  out[order(-abs(out$logFC)), ]
}

organs <- c("Brain", "Liver", "Kidney", "Spleen")
cand <- setNames(lapply(organs, get_candidates), organs)

sel <- bind_rows(
  data.frame(tissue = "Brain",  dir = "up",   symbol = pick_up("Brain",  cand$Brain)$symbol),
  data.frame(tissue = "Brain",  dir = "down", symbol = pick_down("Brain")$symbol),
  data.frame(tissue = "Liver",  dir = "up",   symbol = pick_up("Liver",  cand$Liver)$symbol),
  data.frame(tissue = "Liver",  dir = "down", symbol = pick_down("Liver", force_genes = "TP53")$symbol),
  data.frame(tissue = "Kidney", dir = "up",   symbol = pick_up("Kidney", cand$Kidney)$symbol),
  data.frame(tissue = "Kidney", dir = "down", symbol = pick_down("Kidney")$symbol),
  data.frame(tissue = "Spleen", dir = "up",   symbol = pick_up("Spleen", cand$Spleen)$symbol),
  data.frame(tissue = "Spleen", dir = "down", symbol = pick_down("Spleen", force_genes = "CDKN1A")$symbol)
)
sel$basis <- ifelse(sel$dir == "up", "Displayed GO/KEGG term member (this organ's enrichment figure)",
                     "Senescence marker (if significant) + largest |logFC| down genes")

d <- sel %>% left_join(de_best %>% select(tissue, symbol, logFC, adjP = adj.P.Val),
                        by = c("tissue", "symbol"))
stopifnot(!anyNA(d$logFC), all(sign(d$logFC) == ifelse(d$dir == "up", 1, -1)))
write.csv(d, file.path(OUT, "panelD_v5_geneEffect.csv"), row.names = FALSE)

## ---- 3. plot (text/legend/strip removed,
## per-organ hand-picked x breaks) ----
d$tissue <- factor(d$tissue, levels = c("Brain", "Liver", "Kidney", "Spleen"))
d <- d %>% arrange(tissue, logFC)
d$ykey <- factor(paste(d$tissue, d$symbol, sep = "___"), levels = paste(d$tissue, d$symbol, sep = "___"))
d$dirlab <- factor(ifelse(d$logFC > 0, "Upregulated", "Downregulated"), levels = c("Upregulated", "Downregulated"))
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
d$star <- star(d$adjP); d$hj <- ifelse(d$logFC > 0, -0.25, 1.25)

manual_breaks <- list(Brain = c(-3, 0, 3, 6), Liver = c(-3, 0, 3),
                       Kidney = c(-3, -1.5, 0, 1.5), Spleen = c(-2, 0, 2))
sx <- lapply(split(d, d$tissue), function(z) {
  lo <- min(z$logFC); hi <- max(z$logFC); pad <- 0.28 * (hi - lo)
  scale_x_continuous(limits = c(lo - pad, hi + pad),
                      breaks = manual_breaks[[as.character(z$tissue[1])]], expand = c(0, 0))
})

p <- ggplot(d, aes(logFC, ykey, fill = dirlab)) +
  geom_col(width = 0.74, color = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey30") +
  geom_text(aes(label = star, hjust = hj), size = 3.5, vjust = 0.78, color = "grey20") +
  facet_wrap(~tissue, scales = "free", nrow = 1) +
  facetted_pos_scales(x = sx) +
  scale_y_discrete(labels = function(x) sub(".*___", "", x)) +
  scale_fill_manual(values = c(Upregulated = "#E64B35", Downregulated = "#4DBBD5"), name = NULL) +
  labs(x = "Effect on EGM treatment (log2FC)", y = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "grey30"), axis.ticks = element_line(color = "grey30"),
        panel.spacing = unit(1.1, "lines"),
        axis.text.y = element_text(size = 9.8), axis.title.x = element_text(size = 12))

ggsave(file.path(OUT, "panelD_v5_geneEffect.png"), p, width = 16, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT, "panelD_v5_geneEffect.pdf"), p, width = 16, height = 7, bg = "white")
print(p)
cat("saved to", OUT, "\n")
