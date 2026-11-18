suppressMessages({library(Seurat); library(ggplot2); library(Matrix)})

## ---- load object, define time points ----------------------------------------
root   <- "reproduce"
outdir <- file.path(root, "out")

## ---- 1. load -----------------------------------------------------------------
e <- new.env(); load("SOT/sce_seurat.RData", envir = e)   # ~10 s, ~3.5 GB
so   <- e$sce.seurat
expr <- GetAssayData(so, assay = "originalexp", layer = "data")
pt   <- e$pseudotimeED[colnames(so), 1]

## ---- 2. timepoints ------------------------------------------------------------
# osk_d4 has 0 cells and is dropped; iPSC (31.2) and ESC (31.6) are nearly identical
lab <- c(mef="MEF", osk_d0="Day0", osk_d1="Day1", osk_d2="Day2", osk_d3="Day3",
         osk_d5="Day5", osk_d7="Day7", osk_d8="Day8", ips="iPSC", esc="ESC")
tp <- factor(unname(lab[as.character(so$Day)]), levels = unname(lab))
names(tp) <- colnames(so); so$Timepoint <- droplevels(tp)

## ---- palette ----------------------------------------------------------------
cols10  <- c("#ACA4C9","#F7A58A","#EBEEA2","#4E7BAB","#E74A48","#916367","#935DA7","#F6C883","#FF9F1C","#7FBF70")

## ---- A. panel M ---------------------------------------------------------------
# Module scores are already in meta.data (AddModuleScore(ctrl = 3) in code_SOT1.R).
Idents(so) <- "Timepoint"
ms <- so$young_ear_signature1
pM <- VlnPlot(so, features = "young_ear_signature1", pt.size = 0, cols = cols10) +
  geom_boxplot(outlier.size = 0, fill = "white", width = 0.1) +
  annotate("text", x = 5, y = max(ms) * 1.10, label = "****", size = 5) +
  scale_y_continuous(limits = c(min(ms) * 1.05, max(ms) * 1.22)) +
  labs(title = "EGM-Young Module scoring\n(n=65 genes)", x = "Timepoints", y = "MS value") +
  theme(legend.position = "none", plot.title = element_text(hjust = .5, size = 11, face = "plain"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
ggsave(file.path(outdir, "panelM_EGM_Young_module_scoring.pdf"), pM, width = 3.4, height = 3.4)
# **** = Day3 vs MEF, wilcoxon p = 1.06e-22
wilcox.test(ms[so$Timepoint == "Day3"], ms[so$Timepoint == "MEF"])$p.value
