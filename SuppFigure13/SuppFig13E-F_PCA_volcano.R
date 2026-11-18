# PCA (panel E style) + label-free volcano (panel F style) for the human EGM data.
# Style: open diamonds / plain PC1-PC2 for the PCA;
# -log10(P) volcano with Control grey, EGM blue, n.s. black, dot-dash cut-offs, gene counts, square legend on top.
# No titles.
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(dplyr)})
ROOT <- "human_egm_261007/egm_score"; D <- file.path(ROOT, "data"); FIG <- file.path(ROOT, "figures")
COL <- c(Control = "#8C8C8C", EGM = "#5B9BD5")
box_theme <- theme_classic(base_size = 13) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = .8), axis.line = element_blank(),
        axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 15),
        legend.background = element_blank(), legend.key = element_blank())

## ---- PCA: DESeq2 vst (blind), top 500 variable genes ---------------------
dds <- readRDS(file.path(D, "deseq2_objects.rds"))$dds
vsd <- vst(dds, blind = TRUE)
pc  <- plotPCA(vsd, intgroup = "group", ntop = 500, returnData = TRUE); pv <- round(100 * attr(pc, "percentVar"))
pc$group <- factor(pc$group, c("Control", "EGM"))
write.csv(pc |> transmute(sample = name, group, PC1, PC2), file.path(D, "PCA_vst_top500.csv"), row.names = FALSE)
cat(sprintf("PC1 %d%%, PC2 %d%% of variance\n", pv[1], pv[2]))
LIM <- ceiling(max(abs(c(pc$PC1, pc$PC2))) / 10) * 10; BRK <- seq(-LIM, LIM, 10)   # axes end exactly on a labelled tick (+/-30)
pca <- ggplot(pc, aes(PC1, PC2, colour = group)) +
  geom_point(shape = 5, size = 5.5, stroke = 1.3) +
  scale_colour_manual(values = COL) +
  scale_x_continuous(limits = c(-LIM, LIM), breaks = BRK, expand = c(0, 0)) + scale_y_continuous(limits = c(-LIM, LIM), breaks = BRK, expand = c(0, 0)) +
  labs(x = "PC1", y = "PC2") + box_theme +
  theme(legend.position = "none", plot.margin = margin(8, 12, 8, 8))   # no legend inside the plot
ggsave(file.path(FIG, "PCA_vst_top500.png"), pca, width = 3.8, height = 3.8, dpi = 300, bg = "white")
ggsave(file.path(FIG, "PCA_vst_top500.pdf"), pca, width = 3.8, height = 3.8, bg = "white")

## ---- volcano without gene labels, |lfc| > 1.5 cut-offs, no gene-count annotations ----
## counts (padj < 0.05, apeglm lfc, named genes): |lfc|>0.5 2,118 down / 2,348 up; |lfc|>1 741 / 1,006; |lfc|>1.5 see script output
## NB the enrichment / heatmap work uses |lfc| > 0.5; the volcano lines are display thresholds only.
sig <- read.csv(file.path(D, "EGM_vs_Control_signature.csv")) |> filter(named, !is.na(padj))
PADJ <- 0.05; LFC <- 1.5; CAP <- 100; XR <- 8
v <- sig |> mutate(y = pmin(-log10(pmax(pvalue, 1e-300)), CAP),
                   cls = factor(case_when(padj < PADJ & lfc > LFC ~ "EGM", padj < PADJ & lfc < -LFC ~ "Control", TRUE ~ "ns"), c("Control", "EGM", "ns")))
p_line <- -log10(max(v$pvalue[v$padj < PADJ], na.rm = TRUE))   # raw-P position of the padj = 0.05 boundary
cat(sprintf("|lfc| > %.1f: down %d / up %d\n", LFC, sum(v$cls == "Control"), sum(v$cls == "EGM")))
volc <- ggplot(v, aes(lfc, y)) +
  geom_hline(yintercept = p_line, linetype = "dotdash", linewidth = .45, colour = "grey25") +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dotdash", linewidth = .45, colour = "grey25") +
  geom_point(data = v |> filter(cls == "ns"), colour = "black", size = .6, alpha = .8) +
  geom_point(data = v |> filter(cls != "ns"), aes(colour = cls), size = 1.3, alpha = .85) +
  scale_colour_manual(values = COL, breaks = c("Control", "EGM"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1))) +
  scale_x_continuous(limits = c(-XR, XR), breaks = seq(-8, 8, 4), oob = scales::squish) +   # genes beyond +/-8 sit on the edge
  scale_y_continuous(breaks = seq(0, CAP, 25), labels = function(b) ifelse(b == CAP, paste0("\u2265", CAP), b), expand = expansion(mult = c(.02, .04))) +
  labs(x = expression(Log[2]*FoldChange), y = "-log(Pvalue)") + box_theme +
  theme(legend.position = "top", legend.text = element_text(size = 14), legend.key.size = unit(.5, "cm"), legend.margin = margin(0, 0, -4, 0))
ggsave(file.path(FIG, "volcano_EGM_vs_Control_nolabel.png"), volc, width = 4.6, height = 4.6, dpi = 300, bg = "white")
ggsave(file.path(FIG, "volcano_EGM_vs_Control_nolabel.pdf"), volc, width = 4.6, height = 4.6, bg = "white")
write.csv(v |> filter(cls != "ns") |> transmute(gene = id, direction = ifelse(cls == "EGM", "up in EGM", "down in EGM"), lfc, lfc_raw, baseMean, pvalue, padj) |> arrange(direction, desc(abs(lfc))),
          file.path(ROOT, "tables", sprintf("DEG_EGM_vs_Control_padj0.05_lfc%s.csv", LFC)), row.names = FALSE)
