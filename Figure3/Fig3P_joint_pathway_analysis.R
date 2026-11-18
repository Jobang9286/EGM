# ============================================================
# Joint pathway analysis (RNA-seq + lipidomics), liver EGM vs control
#   1) DESeq2 DEG        : counts_liver.csv (3 vs 3)
#   2) Lipid list        : n_liver_final2.csv, EGM > control (56 lipids)
#   3) MetaboAnalystR    : integrated KEGG pathways (mmu), hypergeometric + degree centrality
#   4) ggplot2           : Pathway impact vs -log10(p) bubble plot
# ============================================================
suppressPackageStartupMessages({
  library(DESeq2); library(MetaboAnalystR); library(ggplot2); library(ggrepel)
})

base_dir <- file.path(getwd(), "LCMS")
jpa_dir  <- file.path(base_dir, "joint_analysis")
out_dir  <- file.path(jpa_dir, "results")      # MetaboAnalystR writes its files here
res_dir  <- file.path(jpa_dir, "resources")    # mmu_genes.sqlite (xialab.ca/resources/sqlite)
dir.create(out_dir, showWarnings = FALSE)

## ---- parameters (TUNABLE) ----
deg_padj <- 0.05
deg_lfc  <- 2          # |log2FC| cutoff; 2 gave the closest match to the MetaboAnalyst web figure
rerun_jpa <- FALSE     # TRUE -> call the xialab.ca API again; FALSE -> reuse results/ if present

## ---- 1. RNA-seq DEG (DESeq2) ----
cnt <- read.csv(file.path(jpa_dir, "counts_liver.csv"), row.names = 1, check.names = FALSE)
colnames(cnt) <- sub("\\.bam$", "", colnames(cnt))
coldata <- data.frame(group = factor(ifelse(grepl("EGM", colnames(cnt)), "EGM", "control"),
                                     levels = c("control", "EGM")),
                      row.names = colnames(cnt))
dds <- DESeqDataSetFromMatrix(cnt, coldata, ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds, quiet = TRUE)
deg_tab <- as.data.frame(results(dds, contrast = c("group", "EGM", "control")))
deg_tab <- cbind(gene = rownames(deg_tab), deg_tab)
write.csv(deg_tab, file.path(out_dir, "DESeq2_EGM_vs_control.csv"), row.names = FALSE)

deg <- subset(deg_tab, !is.na(padj) & padj < deg_padj & abs(log2FoldChange) > deg_lfc)
write.table(deg[, c("gene", "log2FoldChange")], file.path(out_dir, "input_genes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## ---- 2. Lipid list (EGM-up) ----
lip <- read.csv(file.path(base_dir, "n_liver_final2.csv"), check.names = FALSE)
lip$diff <- rowMeans(lip[, grep("EGM", names(lip))]) - rowMeans(lip[, grep("Con", names(lip))])
cmpd <- data.frame(name = lip$name, fc = lip$diff)[lip$diff > 0, ]
# names the joint-pathway mapper misses without lipid = TRUE -> use KEGG standard names
cmpd$name[cmpd$name == "Eicosatrienoic acid"] <- "Dihomo-gamma-linolenic acid"
cmpd$name[cmpd$name == "Eicosenoic acid"]     <- "11Z-Eicosenoic acid"
write.table(cmpd, file.path(out_dir, "input_lipids.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## ---- 3. MetaboAnalystR joint pathway analysis ----
run_jpa <- function(genes_df, cmpd_df, outdir, topo = "dc", enrich = "hyper",
                    libOpt = "integ", integOpt = "query") {
  old <- setwd(outdir); on.exit(setwd(old))
  gene_txt <- paste(paste(genes_df[[1]], signif(genes_df[[2]], 5), sep = "\t"), collapse = "\n")
  cmpd_txt <- paste(paste(cmpd_df[[1]], signif(cmpd_df[[2]], 5), sep = "\t"), collapse = "\n")
  mSet <- InitDataObjects("conc", "pathinteg", FALSE)
  assign("url.pre", paste0(res_dir, "/"), envir = .GlobalEnv)   # local gene-ID sqlite
  mSet <- SetOrganism(mSet, "mmu")
  mSet <- PerformGeneMapping(mSet, gene_txt, "mmu", "symbol")
  mSet <- PerformCmpdMapping(mSet, cmpd_txt, "mmu", "name")
  mSet <- CreateMappingResultTable(mSet)
  mSet <- PrepareIntegData(mSet)
  PerformIntegPathwayAnalysis(mSet, topo, enrich, libOpt, integOpt)   # sends lists to xialab.ca/api
}
path_csv <- file.path(out_dir, "MetaboAnalyst_result_pathway.csv")
if (rerun_jpa || !file.exists(path_csv)) invisible(run_jpa(deg[, c("gene", "log2FoldChange")], cmpd, out_dir))

pm <- read.csv(path_csv, row.names = 1, check.names = FALSE)
pa <- data.frame(pathway = rownames(pm), total = pm$Total, hits = pm$Hits,
                 p = pm[["Raw p"]], fdr = pm$FDR, impact = pm$Impact)
pa$nlp <- -log10(pa$p)
pa$sig <- pa$p < 0.05
write.csv(pa, file.path(jpa_dir, "joint_pathway_table.csv"), row.names = FALSE)

## ---- 4. Plot (TUNABLE) ----
# no grid; significant terms get a number next to the dot (key in joint_pathway_significant.csv)
font <- "Arial"
sig_by  <- "p"         # "fdr" (FDR-adjusted) or "p" (raw p)
sig_cut <- 0.05
lab_size <- 5          # number text size (mm)
ink <- c(text = "black", line = "black")
lw <- 0.6              # line width for axes / ticks; point outline = lw * 0.875
ramp <- c("#FFE66D", "#FDB42F", "#F57B26", "#E03B24", "#B5121B")   # yellow -> red (as MetaboAnalyst)

nlp_lab <- expression(paste("−", log[10], "(", italic(p), ")"))   # text minus -> stays in Arial
size_rng <- c(1.1, 4) * 1.2
fs <- 2                # text size multiplier

pa$pt_size <- scales::area_pal(size_rng)(scales::rescale(pa$hits, from = range(pa$hits)))  # = scale_size mapping
is_sig <- pa[[sig_by]] < sig_cut
sig_tab <- pa[is_sig, ][order(pa$p[is_sig]), ]
sig_tab$no <- seq_len(nrow(sig_tab))
write.csv(sig_tab[, c("no", "pathway", "hits", "total", "p", "fdr", "impact")],
          file.path(jpa_dir, "joint_pathway_significant.csv"), row.names = FALSE)
pa$num <- ""                                   # empty label = dot that numbers are repelled from
pa$num[match(sig_tab$pathway, pa$pathway)] <- sig_tab$no
num_left <- c("Biosynthesis of unsaturated fatty acids",    # numbers placed on the left of the dot
              "Nitrogen metabolism", "alpha-Linolenic acid metabolism")
pa$nudge <- ifelse(pa$num == "", 0, ifelse(pa$pathway %in% num_left, -0.07, 0.07))

pa <- pa[order(pa$nlp), ]                      # most significant drawn on top
g <- ggplot(pa, aes(impact, nlp)) +
  geom_point(aes(size = hits, fill = nlp), shape = 21, colour = ink["line"], stroke = lw * 0.875) +
  geom_text_repel(aes(label = num, point.size = pt_size), family = font, size = lab_size,
                  colour = ink["text"], nudge_x = pa$nudge, min.segment.length = Inf,
                  box.padding = 0.05, point.padding = 0.05, force = 0.5, force_pull = 3,
                  max.overlaps = Inf, seed = 1) +
  scale_fill_gradientn(colours = ramp, name = nlp_lab,
                       guide = guide_colourbar(barheight = unit(28, "mm"), barwidth = unit(4.5, "mm"),
                                               order = 1)) +
  scale_size(range = size_rng, limits = range(pa$hits), breaks = c(2, 4, 6, 8), name = "Hits",
             guide = guide_legend(order = 2, override.aes = list(fill = "#FDB42F"))) +
  scale_x_continuous(limits = c(0, 2), breaks = c(0, 1, 2), labels = c("0", "1.0", "2.0"), expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2), expand = expansion(mult = c(0.03, 0.03))) +
  labs(x = "Pathway impact", y = nlp_lab) +
  theme_classic(base_size = 11 * fs, base_family = font) +
  theme(axis.line  = element_line(colour = ink["line"], linewidth = lw),
        axis.ticks = element_line(colour = ink["line"], linewidth = lw),
        text       = element_text(colour = ink["text"]),
        axis.text  = element_text(colour = ink["text"]),
        axis.title = element_text(colour = ink["text"]),
        legend.title = element_text(colour = ink["text"], size = 9 * fs),
        legend.text  = element_text(colour = ink["text"], size = 8.5 * fs),
        legend.spacing.y = unit(2, "mm"),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(jpa_dir, "joint_pathway_impact.png"), g, width = 7, height = 5.2, dpi = 300, bg = "white")
print(g)   # RStudio Plots pane
