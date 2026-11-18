# ==== packages, counts, DESeq2 (HGPS vs HGPS_EGM, three tissues pooled) ====
library(tximport)
library(biomaRt)
library(Rsubread)
library(DESeq2)
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("edgeR")
library(gplots)
library(pheatmap)
library(RRHO2)
library(UpSetR)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
setwd("39_fibroblasts/2_RNA_seq_tissue/only_egm_effect_in_tissue/")

egm_mat <- as.data.frame(read.csv("counts_egm.csv"))
row.names(egm_mat) <- egm_mat[,"Geneid"]
egm_mat <- as.data.frame(egm_mat[,-1])

col1 <- colnames(egm_mat)
col1

condition = c(rep('HGPS',9),  rep('HGPS_EGM',9))
condition


exp.data <- data.frame(
  row.names=colnames(egm_mat),
  condition = condition)

exp.data

dds <- DESeqDataSetFromMatrix(countData = egm_mat,
                              colData = exp.data,
                              design = ~ condition)

egm_dds <- DESeq(dds)
# ==== fgsea on lipid GO terms (gene selection for 6O) ====
res <- results(egm_dds, contrast = c("condition", "HGPS_EGM", "HGPS"))

write.csv(res, "res.csv")

res <- res[!is.na(res$pvalue) & !is.na(res$log2FoldChange), ]

res

pval_threshold <- 0.05

logfc_threshold <- 0.1

filtered_res <- res[res$pvalue < pval_threshold & abs(res$log2FoldChange) > logfc_threshold, ]

ranks <- filtered_res$log2FoldChange
names(ranks) <- rownames(filtered_res)

ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)


mouse_gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")

filtered_gobp <- mouse_gobp %>%
  filter(grepl("BIOSYNTHESIS|TRANSPORT|METABOLISM", gs_name, ignore.case = TRUE))

filtered_gobp

mouse_gobp_pathways <- split(x = filtered_gobp$gene_symbol, f = filtered_gobp$gs_name)

length(mouse_gobp_pathways)

fgsea_res <- fgsea(
  pathways = mouse_gobp_pathways,
  stats = ranks,
  nperm = 1000
)

significant_pathways <- fgsea_res[fgsea_res$pval < 0.05, ]

significant_pathways

significant_pathways$leadingEdge

write.csv(significant_pathways_mat, "./pathway.csv")

selected_pathways <- c("GOBP_RESPONSE_TO_LIPID", "GOBP_FATTY_ACID_METABOLIC_PROCESS")

filtered_fgsea_res <- fgsea_res %>%
  filter(pathway %in% selected_pathways) %>%
  mutate(pathway = factor(pathway, pathway[order(NES, decreasing = TRUE)]))
# ==== heatmaps (6O: first genes list; 6P: genes_br) ====
rld <- rlog(egm_dds)

egm_dds

genes <- c("Meiosin", "Cxcl2", "Tnip3", "Acta1", "Foxp3", "Ppbp",      "Fos",       "Ffar2",     "Ccr7",      "Nr1d1",     "Atp1a2",    "Heyl",      "Nr4a3",    
"Ppargc1b",  "Map4k1",    "Nr4a1",     "E2f1",      "Adcyap1r1", "Mmp2",      "Tbxas1",    "Ezh2",      "Paqr7",  "Ucp2",  "Irf8", "Acsl1", "Acsl3", "Acsl5", "Elovl3", "Elovl5", "Acsm1", "Acsm3", "Amacr", "Dhrs11", "Hsd3b2", "Fitm2", "Nudt8","Lep","Slc27a6")
mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("condition")])

pheatmap(mat, annotation_col = anno, border_color = "white", cluster_rows = T, cluster_cols = F, scale = "row")

rld <- rlog(egm_dds)

genes <- c("Acsl1", "Acsl3", "Acsl5", "Acsm1", "Acsm3", "Acta1", "Adcyap1r1", "Amacr", 
           "Arsa", "Arsg", "Atp1a2", "Ccr7", "Crybg2", "Cxcl2", "Dhrs11", "E2f1", 
           "Elovl3", "Elovl5", "Ezh2", "Ffar2", "Fitm2", "Fndc11", "Fos", "Foxp3", 
           "Gja10", "Heyl", "Hsd3b2", "Hsd3b5", "Irf8", "Lep", "Map4k1", "Meiosin", 
           "Mmp2", "Nr1d1", "Nr4a1", "Nr4a3", "Nudt8", "Padi4", "Paqr7", "Ppargc1b", 
           "Ppbp", "Slc16a8", "Slc27a6", "Tbxas1", "Tnip3", "Ucp2")

genes_br <- c("Lman2l", "Tmeff2", "Nr4a2", "Agps", "Cyp7b1", "Plppr4", "Plppr5", "Ugt8a", "Elovl6", "Marcksl1", "Miat", "St7", "Akt2", "Rab30", "Agpat5", "Hhip", "Fabp7", "Dhrs13", "Elovl7", "Rab3c", "Atg12", "Tmx3")

mat <- assay(rld)[genes_br, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("condition")])

pheatmap(mat, annotation_col = anno, border_color = "white", cluster_rows = T, cluster_cols = F, scale = "row")
