## ============================ Step 1 =========================================
## lipids increased in EGM liver: mean(EGM) - mean(Con) > 0 (values are log2-scale)
df <- read.csv("LCMS/n_liver_final2.csv", check.names = FALSE)
rownames(df) <- df$name
mat <- as.matrix(df[, -1])
mat <- mat[, c("P_Con_Liver_1","P_Con_Liver_2","P_EGM_Liver_1","P_EGM_Liver_2")]

con_mean <- rowMeans(mat[, c("P_Con_Liver_1","P_Con_Liver_2")])
egm_mean <- rowMeans(mat[, c("P_EGM_Liver_1","P_EGM_Liver_2")])
log2fc <- egm_mean - con_mean   # values already look log2-scale

res <- data.frame(name = rownames(mat), con_mean, egm_mean, log2fc = log2fc)
res <- res[order(-res$log2fc), ]
print(res, digits = 3)

cat("\n=== UP in EGM (log2fc > 0) count:", sum(res$log2fc > 0), "===\n")
cat("=== DOWN in EGM (log2fc < 0) count:", sum(res$log2fc < 0), "===\n")
## The 56 lipids with log2fc > 0 (egm_up_lipids.txt) are the cmpd.vec input of Step 2.

## ============================ Step 2 =========================================
## ORA with MetaboAnalystR (lipid_pathway library, hypergeometric test).
library(MetaboAnalystR)

cmpd.vec <- c("PC(16:0/18:3)","PS(20:0/20:6)","PE(16:0/18:1)","PS(18:0/20:4)","PE(18:0/20:7)",
"PI(18:0/18:4)","PI(18:0/20:5)","PE(18:0/18:1)","PG(16:0/18:2)","PC(16:0/16:1)",
"PC(18:0/20:5)","Eicosenoic acid","PE(20:0/20:5)","PI(18:0/20:4)","PE(16:0/18:2)",
"Eicosadienoic acid","Palmitoleic acid","PC(16:0/16:0)","PE(18:0/18:5e)","LNAPE(20:4/N-18:0)",
"LNAPE(18:1/N-18:0)","alpha-Linolenic acid","PE(20:0/20:6e)","PC(16:0/18:1)","LPE(16:0)",
"PC(18:0/18:4)","Eicosatrienoic acid","LPE(18:0)","PE(19:0/18:3e)","LPE(18:1)",
"LPC(16:0)","PE(18:0/18:2)","LPC(18:1)","SM(10:1;2O/30:0)","Cer-NDS(d40:1)",
"PC(18:0/18:3)","PE(18:0/18:3)","PE(18:0/18:5)","PC(16:0/18:2)","PE(18:0/20:4)",
"PE(18:0/20:6)","LPE(22:6)","Oleic acid","PE(18:0/20:5)","PC(18:0/20:6)",
"PC(16:0/18:0)","PE(18:0/18:4)","Linoleic acid","Docosatetraenoic acid","PC(18:0/20:4)",
"LPC(20:4)","PC(18:0/18:1)","PC(18:0/18:2)","PE(18:0/20:2)","LPC(18:0)","PE(18:0/20:5e)")

setwd("LCMS/enrichment")

mSet <- InitDataObjects("conc", "msetora", FALSE)
mSet <- Setup.MapData(mSet, cmpd.vec)
mSet <- CrossReferencing(mSet, "name", lipid = TRUE)
mSet <- CreateMappingResultTable(mSet)
mSet <- SetMetabolomeFilter(mSet, FALSE)
mSet <- SetCurrentMsetLib(mSet, "lipid_pathway", 2)
mSet <- CalculateHyperScore(mSet, "hyperg")
mSet <- PlotEnrichDotPlot(mSet, "ora", "ora_dot_local_", "png", 150, width=NA)

cat("\nDone.\n")

## ============================ Step 3 =========================================
## read ORA table, top 25 by raw P
library(ggplot2)
library(dplyr)

res <- read.csv("msea_ora_result.csv", check.names = FALSE)
colnames(res)[1] <- "MetaboliteSet"

res$EnrichmentRatio <- res$hits / res$expected
res$neg_log10_p <- -log10(res$`Raw p`)

top25 <- res %>%
  arrange(`Raw p`) %>%
  head(25) %>%
  arrange(neg_log10_p)   # ascending so coord_flip reads top-to-bottom by significance

top25$MetaboliteSet <- factor(top25$MetaboliteSet, levels = top25$MetaboliteSet)

## display labels
library(ggh4x)

top25_disp <- top25
top25_disp$Label <- as.character(top25_disp$MetaboliteSet)
top25_disp$Label[top25_disp$Label == "Synthesis, secretion, and inactivation of Glucagon-like Peptide-1 (GLP-1)"] <- "GLP-1 synthesis, secretion & inactivation"
top25_disp$Label <- factor(top25_disp$Label, levels = top25_disp$Label)

## shortened labels
shorten_map <- c(
  "Transcriptional regulation of white adipocyte differentiation" = "Regulation of white adipocyte differentiation",
  "Defective SLC27A4 causes ichthyosis prematurity syndrome (IPS)" = "Defective SLC27A4 (IPS)",
  "Regulation of cholesterol biosynthesis by SREBP (SREBF)" = "Cholesterol biosynthesis by SREBP",
  "Activation of gene expression by SREBF (SREBP)" = "Gene expression activation by SREBP",
  "Transcriptional activation of mitochondrial biogenesis" = "Mitochondrial biogenesis (transcriptional)",
  "Regulation of lipid metabolism by PPARalpha" = "Lipid metabolism by PPARalpha",
  "Linoleic acid metabolism affected by SARS-CoV-2" = "Linoleic acid metabolism (SARS-CoV-2)"
)

top25_disp$Label <- as.character(top25_disp$Label)
for (nm in names(shorten_map)) {
  top25_disp$Label[top25_disp$Label == nm] <- shorten_map[[nm]]
}
# restore correct significance order (was ordered ascending neg_log10_p earlier)
top25_disp <- top25_disp[order(top25_disp$neg_log10_p), ]
top25_disp$Label <- factor(top25_disp$Label, levels = top25_disp$Label)

## final plot
p <- ggplot(top25_disp, aes(x = neg_log10_p, y = Label)) +
  geom_point(aes(size = EnrichmentRatio, color = `Raw p`)) +
  scale_color_gradient(low = "red", high = "orange", name = "P-value") +
  scale_size_continuous(name = "Enrichment Ratio", range = c(2, 8)) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  labs(title = "Overview of Enriched Metabolite Sets (Top 25)",
       x = "-log10 (p-value)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(family = "Arial", colour = "black"),
    plot.title = element_text(size = 13, family = "Arial", colour = "black"),
    axis.text.y = element_text(size = 13, family = "Arial", colour = "black"),
    axis.text.x = element_text(size = 11, family = "Arial", colour = "black"),
    axis.title.x = element_text(colour = "black"),
    legend.text = element_text(colour = "black"),
    legend.title = element_text(colour = "black")
  ) +
  force_panelsizes(cols = unit(3.3, "in"), rows = unit(6.2, "in"))

ggsave("enrichment_dotplot_ggplot2.pdf", p,
       width = 10, height = 7.2, device = cairo_pdf)
ggsave("enrichment_dotplot_ggplot2.png", p,
       width = 10, height = 7.2, dpi = 300, type = "cairo")
cat("saved.\n")
