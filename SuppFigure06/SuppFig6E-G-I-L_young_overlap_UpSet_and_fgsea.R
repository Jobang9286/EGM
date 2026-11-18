# ==== packages, liver counts, liver DESeq2 ====
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
setwd("39_fibroblasts/2_RNA_seq_tissue/")

#####################liver final###########################

liv_mat <- as.data.frame(read.csv("Liver_final/counts_liver_revised.csv"))
row.names(liv_mat) <- liv_mat[,"Geneid"]
liv_mat <- as.data.frame(liv_mat[,-1])

liv_col1 <- colnames(liv_mat)
liv_col1

condition = c(rep('Liver_HGPS_control',3),  rep('Liver_HGPS_EGM',3), rep('Liver_WT_YOUNG',6), rep('Liver_WT_AGED',7), rep('Liver_WT_YOUNG',6))
condition = c(rep('Liver_HGPS_control',3),  rep('Liver_HGPS_EGM',3), rep('Liver_WT_1M',6), rep('Liver_WT_24M',3), rep('Liver_WT_27M',4), rep('Liver_WT_3M',6))


condition


liv_mat <- as.data.frame(read.csv("Liver_final/counts_liver_only_egm.csv"))
row.names(liv_mat) <- liv_mat[,"Geneid"]
liv_mat <- as.data.frame(liv_mat[,-1])

liv_col1 <- colnames(liv_mat)
liv_col1

condition = c(rep('Liver_HGPS_control',3),  rep('Liver_HGPS_EGM',3))
condition

exp.data <- data.frame(
  row.names=colnames(liv_mat),
  condition = condition)

exp.data

liver_dds <- DESeqDataSetFromMatrix(countData = liv_mat,
                              colData = exp.data,
                              design = ~ condition)

liver_dds <- DESeq(liver_dds)
# ==== liver fgsea NES bars (Supp 6L liver) ====
#################################################################################################################################

res <- results(liver_dds, contrast = c("condition", "Liver_HGPS_EGM", "Liver_HGPS_control"))
ranks <- res$log2FoldChange
names(ranks) <- rownames(res)
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

mouse_gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")

filtered_gobp <- mouse_gobp %>%
  filter(grepl("DNA_REPAIR|CYCLE|DEVELOPMENT|DIVISION|METABOLISM|LIVER|WOUND|STRESS|SIGNALING", gs_name, ignore.case = TRUE))

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

selected_pathways <- c("GOBP_CELL_CYCLE", "GOBP_DNA_REPAIR", "GOBP_POSITIVE_REGULATION_OF_CELL_DEVELOPMENT", "GOBP_CELL_DIVISION", "GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY")

filtered_fgsea_res <- fgsea_res %>%
  filter(pathway %in% selected_pathways) %>%
  mutate(pathway = factor(pathway, pathway[order(NES, decreasing = TRUE)]))

ggplot(filtered_fgsea_res, aes(x = pathway, y = NES, fill = -log10(pval))) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  coord_flip() +
  labs(
    title = "Top Enriched Pathways",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_gradient(
    low = "pink",    # 연한 빨간색
    high = "red"     # 진한 빨간색
  ) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

ggplot(filtered_fgsea_res, aes(x = pathway, y = NES, fill = -log10(pval))) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +  # 바 테두리
  labs(
    title = "Top Enriched Pathways",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_gradient(
    low = "pink",    # 연한 빨간색
    high = "red"     # 진한 빨간색
  ) +
  scale_y_continuous(
    breaks = seq(0, 2, by = 0.5),  # y축 눈금 표시 간격 설정
    limits = c(0, 2)               # y축 범위 설정
  ) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),        # y축 레이블 색상
    axis.ticks.y = element_line(color = "black", size = 0.5),     # y축 눈금선 색상
    axis.line.y = element_line(color = "black", size = 0.8),      # y축 축선 색상
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), # x축 레이블 스타일
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # 그래프 바깥쪽 테두리
  )

# ==== liver UpSet + hypergeometric test (Supp 6E) ====
###################################################################################################

res_condition1 <- results(liver_dds, contrast=c("condition","Liver_WT_YOUNG","Liver_WT_AGED"))
res_condition2 <- results(liver_dds, contrast=c("condition","Liver_HGPS_EGM","Liver_HGPS_control"))
res_condition3 <- results(liver_dds, contrast=c("condition","Liver_WT_MID","Liver_WT_AGED"))
res_condition4 <- results(liver_dds, contrast=c("condition","Liver_HGPS_EGM","Liver_WT_YOUNG"))
res_condition5 <- results(liver_dds, contrast=c("condition","Liver_HGPS_EGM","Liver_WT_MID"))
res_condition6 <- results(liver_dds, contrast=c("condition","Liver_HGPS_EGM","Liver_WT_AGED"))



Liver_Y_vs_A <- rownames(subset(res_condition1, pvalue < 0.05 & log2FoldChange > 0.5))
Liver_EGM_vs_CON <- rownames(subset(res_condition2, pvalue < 0.05 & log2FoldChange > 0.5))

deg_list <- list(
  Condition1 = Liver_Y_vs_A,
  Condition2 = Liver_EGM_vs_CON
)

upset(fromList(deg_list), nsets = length(deg_list))

install.packages("ComplexUpset")
library(ComplexUpset)

# 데이터 준비
upset_data <- fromList(deg_list)

# ComplexUpset 플롯 생성
upset_plot <- ComplexUpset::upset(
  upset_data,
  intersect = colnames(upset_data),
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 0))))

upset_plot + 
  theme(
    # Intersection Size 축 (Y축) 눈금선
    panel.grid.major.y = element_line(color = "black", linetype = "dotted", linewidth = 0.5),
    # Set Size 축 (X축) 눈금선
    panel.grid.major.x = element_line(color = "black", linetype = "dotted", linewidth = 0.5)
  )

all_genes <- unique(c(Liver_Y_vs_A, Liver_EGM_vs_CON))

# 3. 바이너리 매트릭스 생성
binary_matrix <- data.frame(
  Condition1 = as.integer(all_genes %in% Liver_Y_vs_A),
  Condition2 = as.integer(all_genes %in% Liver_EGM_vs_CON),
  row.names = all_genes
)

write.csv(binary_matrix, "final/binary_matrix_liver_final_1204.csv")

N <- length(rownames(res_condition1))  # 분석된 전체 유전자 수 (일반적으로 res_condition1, res_condition2에서 동일)

# 각 조건의 DEG 유전자 수
K <- length(Liver_Y_vs_A)        # Condition1 (Liver_Y_vs_A)
n <- length(Liver_EGM_vs_CON)    # Condition2 (Liver_EGM_vs_CON)

# 교집합 유전자 수 계산
x <- length(intersect(Liver_Y_vs_A, Liver_EGM_vs_CON))

# 하이퍼기하분포 검정 p-value 계산
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)

# 결과 출력
cat("Number of overlapping genes:", x, "\n")
cat("P-value:", p_value, "\n")

# ==== brain counts, DESeq2, UpSet + hypergeometric test (Supp 6I) ====
#################################################################################################################################


brain_mat <- as.data.frame(read.csv("Brain_final/bam/counts_brain_revised.csv"))
row.names(brain_mat) <- brain_mat[,"Geneid"]
brain_mat <- as.data.frame(brain_mat[,-1])

col1 <- colnames(brain_mat)
col1

condition = c(rep('Brain_HGPS_control',3),  rep('Brain_HGPS_EGM',3), rep('Brain_WT_YOUNG',4), rep('Brain_WT_AGED',8), rep('Brain_WT_YOUNG',4))
condition

brain_mat <- as.data.frame(read.csv("Brain_final/bam/counts_brain_only_egm.csv"))
row.names(brain_mat) <- brain_mat[,"Geneid"]
brain_mat <- as.data.frame(brain_mat[,-1])

col1 <- colnames(brain_mat)
col1

condition = c(rep('Brain_HGPS_control',3),  rep('Brain_HGPS_EGM',3))
condition


exp.data <- data.frame(
  row.names=colnames(brain_mat),
  condition = condition)

exp.data

dds <- DESeqDataSetFromMatrix(countData = brain_mat,
                              colData = exp.data,
                              design = ~ condition)

brain_dds <- DESeq(dds)

vsd <- vst(brain_dds)
head(colData(vsd),7)

vsd <- vst(brain_dds, blind = T)

vsd_mat <- assay(vsd)

vsd_cor <- cor(vsd_mat)

vsd_cor

pheatmap(vsd_cor)

plotPCA(vsd, "condition")

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

library(ggplot2)
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(c(-20, 20)) + 
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


res_condition1 <- results(brain_dds, contrast=c("condition","Brain_WT_YOUNG","Brain_WT_AGED"))
res_condition2 <- results(brain_dds, contrast=c("condition","Brain_HGPS_EGM","Brain_HGPS_control"))
res_condition3 <- results(brain_dds, contrast=c("condition","Brain_WT_MID","Brain_WT_AGED"))
res_condition4 <- results(brain_dds, contrast=c("condition","Brain_HGPS_EGM","Brain_WT_YOUNG"))
res_condition5 <- results(brain_dds, contrast=c("condition","Brain_HGPS_EGM","Brain_WT_MID"))
res_condition6 <- results(brain_dds, contrast=c("condition","Brain_HGPS_EGM","Brain_WT_AGED"))

Brain_Y_vs_A <- rownames(subset(res_condition1, pvalue < 0.05 & log2FoldChange > 0.1))
Brain_EGM_vs_CON <- rownames(subset(res_condition2, pvalue < 0.05 & log2FoldChange > 0.1))
Brain_M_vs_A <- rownames(subset(res_condition3, pvalue < 0.05 & abs(log2FoldChange) > 0.5))
Brain_EGM_vs_A <- rownames(subset(res_condition6, pvalue < 0.05 & log2FoldChange > 1))

deg_list <- list(
  Condition1 = Brain_Y_vs_A,
  Condition2 = Brain_EGM_vs_CON
)

upset(fromList(deg_list), nsets = length(deg_list))

install.packages("ComplexUpset")
library(ComplexUpset)

# 데이터 준비
upset_data <- fromList(deg_list)

# ComplexUpset 플롯 생성
upset_plot <- ComplexUpset::upset(
  upset_data,
  intersect = colnames(upset_data),
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 5))))

upset_plot + 
  theme(
    # Intersection Size 축 (Y축) 눈금선
    panel.grid.major.y = element_line(color = "black", linetype = "dotted", linewidth = 0.5),
    # Set Size 축 (X축) 눈금선
    panel.grid.major.x = element_line(color = "black", linetype = "dotted", linewidth = 0.5)
  )

all_genes <- unique(c(Brain_Y_vs_A, Brain_EGM_vs_CON))

# 3. 바이너리 매트릭스 생성
binary_matrix <- data.frame(
  Condition1 = as.integer(all_genes %in% Brain_Y_vs_A),
  Condition2 = as.integer(all_genes %in% Brain_EGM_vs_CON),
  row.names = all_genes
)

write.csv(binary_matrix, "final/binary_matrix_brain_final_1204.csv")

N <- length(rownames(res_condition1))  # 분석된 전체 유전자 수 (일반적으로 res_condition1, res_condition2에서 동일)

# 각 조건의 DEG 유전자 수
K <- length(Brain_Y_vs_A)        # Condition1 (Liver_Y_vs_A)
n <- length(Brain_EGM_vs_CON)    # Condition2 (Liver_EGM_vs_CON)

# 교집합 유전자 수 계산
x <- length(intersect(Brain_Y_vs_A, Brain_EGM_vs_CON))

# 하이퍼기하분포 검정 p-value 계산
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)

# 결과 출력
cat("Number of overlapping genes:", x, "\n")
cat("P-value:", p_value, "\n")

# ==== brain fgsea NES bars (Supp 6L brain) ====
#################################################################################################################

res <- results(brain_dds, contrast = c("condition", "Brain_HGPS_EGM", "Brain_HGPS_control"))
res
ranks <- res$log2FoldChange
names(ranks) <- rownames(res)
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

mouse_gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")

filtered_gobp <- mouse_gobp %>%
  filter(grepl("PATHWAY|NEURAL|RESPONSE", gs_name, ignore.case = TRUE))


filtered_gobp

mouse_gobp_pathways <- split(x = filtered_gobp$gene_symbol, f = filtered_gobp$gs_name)

length(mouse_gobp_pathways)

fgsea_res <- fgsea(
  pathways = mouse_gobp_pathways,
  stats = ranks,
  nperm = 2000
)

significant_pathways <- fgsea_res[fgsea_res$pval < 0.05, ]

significant_pathways

selected_pathways <- c("GOBP_NEURAL_CREST_CELL_DIFFERENTIATION", "GOBP_NEURAL_CREST_CELL_MIGRATION", "GOBP_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "GOBP_RESPONSE_TO_LECTIN")

filtered_fgsea_res <- fgsea_res %>%
  filter(pathway %in% selected_pathways) %>%
  mutate(pathway = factor(pathway, pathway[order(NES, decreasing = TRUE)]))

ggplot(filtered_fgsea_res, aes(x = pathway, y = NES, fill = -log10(pval))) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  coord_flip() +
  labs(
    title = "Top Enriched Pathways",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_gradient(
    low = "pink",    # 연한 빨간색
    high = "red"     # 진한 빨간색
  ) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

ggplot(filtered_fgsea_res, aes(x = pathway, y = NES, fill = -log10(pval))) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +  # 바 테두리
  labs(
    title = "Top Enriched Pathways",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_gradient(
    low = "pink",    # 연한 빨간색
    high = "red"     # 진한 빨간색
  ) +
  scale_y_continuous(
    breaks = seq(-2, 2, by = 0.5),  # y축 눈금 표시 간격 설정
    limits = c(-2, 2)               # y축 범위 설정
  ) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks.y = element_line(color = "black", size = 0.5),     
    axis.line.y = element_line(color = "black", size = 0.8),      
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# ==== kidney counts, DESeq2, UpSet + hypergeometric test (Supp 6G); kidney-only DESeq2 for GSEA ====
##############all#########################################################################################

kid_mat <- as.data.frame(read.csv("Kidney_final/total/counts_kidney_time.csv"))
kid_mat <- as.data.frame(read.csv("Kidney_final/total/counts_kidney_time.csv"))

row.names(kid_mat) <- kid_mat[,"Geneid"]
kid_mat <- as.data.frame(kid_mat[,-1])

col1 <- colnames(kid_mat)
col1

condition = c(rep('Kidney_HGPS_control_1',1),  rep('Kidney_HGPS_control',2), rep('Kidney_HGPS_EGM',3),rep('Kidney_WT_YOUNG',4), rep('Kidney_WT_AGED',7), rep('Kidney_WT_YOUNG',3), rep('Kidney_WT_6M',4))
condition

condition = c(rep('Kidney_HGPS_control',3), rep('Kidney_HGPS_EGM',3), rep('Kidney_WT_1M',4), rep('Kidney_WT_24M',3), rep('Kidney_WT_27M',4), rep('Kidney_WT_3M',3), rep('Kidney_WT_6M',4))
condition

condition = c(rep('Kidney_WT_OLD',3), rep('Kidney_WT_YOUNG',5), rep('Kidney_HGPS_control',3), rep('Kidney_HGPS_EGM',3))
condition

exp.data <- data.frame(
  row.names=colnames(kid_mat),
  condition = condition)

exp.data

kid_dds <- DESeqDataSetFromMatrix(countData = kid_mat,
                              colData = exp.data,
                              design = ~ condition)

kid_dds

kid_time_dds <- DESeq(kid_dds)

vsd <- vst(kid_dds)

head(colData(vsd),7)

vsd <- vst(kid_dds, blind = T)

vsd_mat <- assay(vsd)

vsd_cor <- cor(vsd_mat)

vsd_cor

pheatmap(vsd_cor)

plotPCA(vsd, "condition")

res_condition1 <- results(kid_dds, contrast=c("condition","Kidney_WT_YOUNG","Kidney_WT_AGED"))
res_condition2 <- results(kid_dds, contrast=c("condition","Kidney_HGPS_EGM","Kidney_HGPS_control"))
res_condition3 <- results(kid_dds, contrast=c("condition","Brain_WT_MID","Brain_WT_AGED"))
res_condition4 <- results(kid_dds, contrast=c("condition","Brain_HGPS_EGM","Brain_WT_YOUNG"))
res_condition5 <- results(kid_dds, contrast=c("condition","Brain_HGPS_EGM","Brain_WT_MID"))
res_condition6 <- results(kid_dds, contrast=c("condition","Kidney_HGPS_EGM","Kidney_WT_AGED"))

Kidney_Y_vs_A <- rownames(subset(res_condition1, pvalue < 0.05 & abs(log2FoldChange) > 0.5))
Kidney_EGM_vs_CON <- rownames(subset(res_condition2, pvalue < 0.05 & abs(log2FoldChange) > 0.5))
Brain_M_vs_A <- rownames(subset(res_condition3, pvalue < 0.05 & abs(log2FoldChange) > 0.5))
Kidney_EGM_vs_A <- rownames(subset(res_condition6, pvalue < 0.05 & abs(log2FoldChange) > 2))

Kidney_Y_vs_A <- rownames(subset(res_condition1, pvalue < 0.05 & log2FoldChange > 0.2))
Kidney_EGM_vs_CON <- rownames(subset(res_condition2, pvalue < 0.05 & log2FoldChange > 1))
Kidney_EGM_vs_A <- rownames(subset(res_condition6, pvalue < 0.05 & log2FoldChange > 2))


deg_list <- list(
  Condition1 = Kidney_Y_vs_A,
  Condition2 = Kidney_EGM_vs_CON
)

upset(fromList(deg_list), nsets = length(deg_list))

install.packages("ComplexUpset")
library(ComplexUpset)

# 데이터 준비
upset_data <- fromList(deg_list)

# ComplexUpset 플롯 생성
upset_plot <- ComplexUpset::upset(
  upset_data,
  intersect = colnames(upset_data),
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 3))))

upset_plot + 
  theme(
    # Intersection Size 축 (Y축) 눈금선
    panel.grid.major.y = element_line(color = "black", linetype = "dotted", linewidth = 0.5),
    # Set Size 축 (X축) 눈금선
    panel.grid.major.x = element_line(color = "black", linetype = "dotted", linewidth = 0.5)
  )

all_genes <- unique(c(Kidney_Y_vs_A, Kidney_EGM_vs_CON))

# 3. 바이너리 매트릭스 생성
binary_matrix <- data.frame(
  Condition1 = as.integer(all_genes %in% Kidney_Y_vs_A),
  Condition2 = as.integer(all_genes %in% Kidney_EGM_vs_CON),
  row.names = all_genes
)

write.csv(binary_matrix, "final/binary_matrix_kidney_final_1206.csv")

N <- length(rownames(res_condition1))  # 분석된 전체 유전자 수 (일반적으로 res_condition1, res_condition2에서 동일)

# 각 조건의 DEG 유전자 수
K <- length(Kidney_Y_vs_A)        # Condition1 (Liver_Y_vs_A)
n <- length(Kidney_EGM_vs_CON)    # Condition2 (Liver_EGM_vs_CON)

# 교집합 유전자 수 계산
x <- length(intersect(Kidney_Y_vs_A, Kidney_EGM_vs_CON))

# 하이퍼기하분포 검정 p-value 계산
p_value <- phyper(x - 1, K, N - K, n, lower.tail = FALSE)

# 결과 출력
cat("Number of overlapping genes:", x, "\n")
cat("P-value:", p_value, "\n")

kid_mat <- as.data.frame(read.csv("Kidney_final/total/counts_kidney_only_egm.csv"))
row.names(kid_mat) <- kid_mat[,"Geneid"]
kid_mat <- as.data.frame(kid_mat[,-1])

kid_col1 <- colnames(kid_mat)
kid_col1

condition = c(rep('Liver_HGPS_control',3),  rep('Liver_HGPS_EGM',3))
condition

exp.data <- data.frame(
  row.names=colnames(kid_mat),
  condition = condition)

exp.data

kid_dds <- DESeqDataSetFromMatrix(countData = kid_mat,
                                    colData = exp.data,
                                    design = ~ condition)

kid_dds <- DESeq(kid_dds)
