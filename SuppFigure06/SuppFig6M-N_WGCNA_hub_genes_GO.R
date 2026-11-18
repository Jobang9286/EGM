library(WGCNA)
library(readr)
library(dplyr)
library(ggplot2)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

egm_mat <- as.data.frame(read.csv("39_fibroblasts/2_RNA_seq_tissue/only_egm_effect_in_tissue/counts_egm_kid_liv.csv"))
row.names(egm_mat) <- egm_mat[,"Geneid"]
egm_mat <- as.data.frame(egm_mat[,-1])

col1 <- colnames(egm_mat)
col1

condition = c(rep('HGPS',6),  rep('HGPS_EGM',6))
condition


exp.data <- data.frame(
  row.names=colnames(egm_mat),
  condition = condition)

exp.data

dds <- DESeqDataSetFromMatrix(countData = egm_mat,
                              colData = exp.data,
                              design = ~ condition)

dds

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)
rld <- rlog(dds)  # 로그 변환
normalizedCounts <- assay(rld)  # 정규화된 발현 데이터

# WGCNA 입력 데이터 준비
datExpr <- t(normalizedCounts)  # 샘플이 행, 유전자가 열로 변환

# 결측치 확인 및 제거
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 샘플 클러스터링: 이상치 확인
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")

powers <- c(1:20)  # 테스트할 파워 값
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 결과 시각화
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit",
     type = "b", main = "Soft Thresholding Power Selection")
abline(h = 0.9, col = "red")  # R^2 기준선

softPower <- 6  # 최적의 파워 값

# Adjacency 및 TOM 계산
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 유전자 클러스터링
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 모듈 탐색
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
dynamicColors <- labels2colors(dynamicMods)

dynamicColors

# 결과 시각화
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    main = "Gene dendrogram and module colors")

plotDendroAndColors(
  geneTree,                  # 계층적 클러스터링 결과
  colors = dynamicColors,    # 모듈 색상
  groupLabels = "Module Colors", 
  main = "Gene Dendrogram with Module Colors",
  dendroLabels = FALSE,      # 유전자 레이블 제거
  hang = 0.03,               # 선의 길이 조정
  addGuide = TRUE,           # 가이드 선 추가
  guideHang = 0.05           # 가이드 선의 길이 설정
)

selectedColors <- c("blue", "turquoise", "brown")  # 관심 모듈 색상
plotDendroAndColors(geneTree, 
                    colors = dynamicColors, 
                    groupLabels = "Dynamic Tree Cut",
                    main = "Gene Dendrogram and Selected Module Colors",
                    dendroLabels = FALSE,  # 유전자 레이블 제거
                    hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)

# 모듈 특징 벡터(MEs) 생성
MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes

# MEs와 traitData 삭제
rm(MEs, traitData)

# Module Eigengenes 재계산
MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes
print(head(MEs))

# traitData 확인 및 condition 열 존재 확인
traitData <- colData(dds)  # DESeq2의 메타데이터 사용
print(head(traitData))

moduleTraitCor <- cor(MEs, as.numeric(traitData$condition), use = "p")

# 결과 확인
print(moduleTraitCor)

traitData$condition <- factor(c(rep("Control", 6), rep("Treatment", 6)))  # 예시
print(traitData$condition)

textMatrix <- signif(moduleTraitCor, 2)  # 소수점 2자리까지 표시
dim(textMatrix) <- dim(moduleTraitCor)  # 행렬의 크기 맞추기

# textMatrix 확인
print(textMatrix)

# 마진 조정
par(mar = c(4, 4, 2, 2))  # 아래, 왼쪽, 위, 오른쪽 여백 설정

# Heatmap 다시 실행
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "Condition",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               main = "Module-Trait Relationships")

##############################################################

# Module Eigengene 계산
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# 표현형 데이터 준비
traitData <- colData(dds)[, "condition", drop = FALSE]

traitData

traitData$condition <- ifelse(traitData$condition == "HGPS_EGM", 1, 0)

# 모듈-표현형 상관관계 분석
moduleTraitCor <- cor(MEs, traitData$condition, use = "p")

moduleTraitCor

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "HGPS_EGM",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               main = "Module-Trait Relationships")

##############################################################
# condition 값을 숫자로 변환 (예: Control=0, Treatment=1)
traitData$condition_numeric <- as.numeric(factor(traitData$condition))

# 모듈-표현형 상관관계 분석
moduleTraitCor <- cor(MEs, traitData$condition_numeric, use = "p")
textMatrix <- signif(moduleTraitCor, 2)

# 시각화

table(dynamicColors)

head(dynamicColors)
length(dynamicColors)

geneModuleTable <- data.frame(
  Gene = colnames(datExpr),  # 유전자 이름
  ModuleColor = dynamicColors  # 모듈 색상
)

# 결과 확인
head(geneModuleTable)

write.csv(geneModuleTable, "./table.csv")

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

# 유전자와 모듈 색상 연결
geneModuleTable <- data.frame(
  Gene = colnames(datExpr),          # 유전자 이름
  ModuleColor = dynamicColors        # 모듈 색상
)

# 모듈 색상별 유전자군 추출
moduleColors <- unique(dynamicColors)  # 모듈 색상 목록
geneLists <- lapply(moduleColors, function(color) {
  geneModuleTable$Gene[geneModuleTable$ModuleColor == color]
})
names(geneLists) <- moduleColors  # 리스트에 색상 이름 부여

# Symbol → Entrez ID 변환 함수
convertToEntrez <- function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
}

# 각 모듈별로 Entrez ID 변환
geneLists_entrez <- lapply(geneLists, function(genes) {
  res <- convertToEntrez(genes)
  return(res$ENTREZID)  # Entrez ID만 추출
})

performGO <- function(entrezIDs, moduleColor) {
  enrichGO(gene = entrezIDs,
           OrgDb = org.Mm.eg.db,
           keyType = "ENTREZID",
           ont = "BP",  # Biological Process
           pAdjustMethod = "BH",
           qvalueCutoff = 0.05,
           readable = TRUE)
}

goResults <- lapply(names(geneLists_entrez), function(color) {
  cat("Running GO analysis for module:", color, "\n")
  performGO(geneLists_entrez[[color]], color)
})
names(goResults) <- names(geneLists_entrez)

for (color in names(goResults)) {
  cat("Checking GO terms for module:", color, "\n")
  if (!is.null(goResults[[color]])) {
    print(head(goResults[[color]]@result))  # GO term 결과 확인
  } else {
    cat("No GO terms found for module:", color, "\n")
  }
}    

for (color in names(goResults)) {
  if (!is.null(goResults[[color]]) && nrow(goResults[[color]]@result) > 0) {
    cat("Plotting barplot for module:", color, "\n")
    barplot(goResults[[color]], showCategory = 10, title = paste("GO Enrichment for", color, "Module"))
  } else {
    cat("No significant GO terms for module:", color, "\n")
  }
}


target_module <- "turquoise"

target_genes <- geneModuleTable$Gene[geneModuleTable$ModuleColor == target_module]
target_genes_entrez <- bitr(target_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
print(head(target_genes_entrez))
target_GO <- enrichGO(
  gene = target_genes_entrez$ENTREZID,  # Entrez ID 입력
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",            # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,   # p-value 유의성 기준
  qvalueCutoff = 0.2,
  readable = TRUE        # 결과를 Symbol로 출력
)

head(target_GO)



barplot(
  target_GO, 
  showCategory = 15,                 # 상위 10개 GO terms 표시
  title = paste("GO Enrichment for", target_module, "Module")
)

# GO 결과를 데이터프레임으로 변환
go_df <- as.data.frame(target_GO)

# Gene Count를 기준으로 내림차순 정렬
go_df <- go_df[order(go_df$Count, decreasing = TRUE), ]

# 상위 10개의 GO Term 선택
top10_GO <- go_df[1:15, ]

# p.adjust 값 구간 설정
top10_GO$p.adjust_group <- cut(top10_GO$p.adjust, 
                               breaks = c(0, 0.01, 0.05, 0.1, 1), 
                               labels = c("<0.01", "0.01-0.05", "0.05-0.1", ">0.1"))

# 수동 색상 설정
ggplot(top10_GO, aes(x = reorder(Description, -Count), y = Count, fill = p.adjust_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("<0.01" = "#40E0D0", "0.01-0.05" = "orange", 
                               "0.05-0.1" = "yellow", ">0.1" = "gray")) +
  labs(title = "GO Enrichment (Gene Count sorted, p.adjust groups)",
       x = "GO Terms",
       y = "Gene Count",
       fill = "p.adjust") +
  theme_minimal()


moduleGenes <- (dynamicColors == "turquoise")  # 예: "blue" 모듈

traitData$condition_numeric <- as.numeric(factor(traitData$condition))

# 변환 결과 확인
print(traitData$condition_numeric)
print(unique(traitData$condition_numeric))

geneSignificance <- abs(cor(datExpr, traitData$condition_numeric, use = "p"))
moduleMembership <- abs(cor(datExpr, MEs$MEblue, use = "p"))

# 기준 설정
threshold_MM <- 0.7  # Module Membership 기준
threshold_GS <- 0.2  # Gene Significance 기준

# 기준을 충족하는 점을 다른 색상으로 설정
point_colors <- ifelse(moduleMembership[moduleGenes] > threshold_MM & 
                         geneSignificance[moduleGenes] > threshold_GS, "red", "black")

# Scatter plot 그리기
plot(moduleMembership[moduleGenes], geneSignificance[moduleGenes],
     xlab = "Module Membership in turquiose Module",
     ylab = "Gene Significance for Condition",
     main = "Module Membership vs. Gene Significance",
     col = point_colors, pch = 19, cex = 0.6)

# 기준선을 추가 (빨간색 점선)
abline(v = threshold_MM, col = "red", lty = 2, lwd = 2)  # X축 기준선
abline(h = threshold_GS, col = "red", lty = 2, lwd = 2)  # Y축 기준선

hubGenes <- names(moduleMembership[moduleGenes])[moduleMembership[moduleGenes] > threshold_MM & 
                                                   geneSignificance[moduleGenes] > threshold_GS]

hubGenes
print(hubGenes)

# moduleMembership과 geneSignificance 구조 확인
print(head(moduleMembership[moduleGenes]))
print(head(geneSignificance[moduleGenes]))

# 이름이 있는지 확인
print(names(moduleMembership[moduleGenes]))
print(names(geneSignificance[moduleGenes]))

names(moduleMembership) <- colnames(datExpr)
names(geneSignificance) <- colnames(datExpr)

# 확인
print(head(names(moduleMembership)))
print(head(names(geneSignificance)))

threshold_MM <- 0.7
threshold_GS <- 0.2

hubGenes <- names(moduleMembership[moduleGenes])[moduleMembership[moduleGenes] > threshold_MM & 
                                                   geneSignificance[moduleGenes] > threshold_GS]

# 결과 확인
print(hubGenes)

# SYMBOL → Entrez ID 변환
hubGenes_entrez <- bitr(hubGenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# 결과 확인
print(hubGenes_entrez)

# GO Enrichment 분석
ego_hub <- enrichGO(
  gene = hubGenes_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",             # Biological Process
  pAdjustMethod = "BH",   # p-value 조정 방식
  pvalueCutoff = 0.05,    # p-value 기준
  qvalueCutoff = 0.2,
  readable = TRUE         # Entrez ID를 SYMBOL로 변환
)

# 결과 확인
head(ego_hub)

# Barplot
barplot(ego_hub, showCategory = 6, title = "GO Enrichment for Hub Genes")

barplot(ego_hub, 
        showCategory = 6, 
        title = "GO Enrichment for Hub Genes", 
        horizontal = TRUE)




go_df <- as.data.frame(ego_hub)

# Gene Count를 기준으로 내림차순 정렬
go_df <- go_df[order(go_df$Count, decreasing = TRUE), ]

# 상위 10개의 GO Term 선택
top10_GO <- go_df[1:8, ]

# p.adjust 값 구간 설정
top10_GO$p.adjust_group <- cut(top10_GO$p.adjust, 
                               breaks = c(0, 0.01, 0.05, 0.1, 1), 
                               labels = c("<0.01", "0.01-0.05", "0.05-0.1", ">0.1"))

# 수동 색상 설정
ggplot(top10_GO, aes(x = reorder(Description, -Count), y = Count, fill = p.adjust_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("<0.01" = "#40E0D0", "0.01-0.05" = "orange", 
                               "0.05-0.1" = "yellow", ">0.1" = "gray")) +
  labs(title = "GO Enrichment (Gene Count sorted, p.adjust groups)",
       x = "GO Terms",
       y = "Gene Count",
       fill = "p.adjust") +
  theme_minimal()

# 특정 GO term의 Description 또는 ID
target_GO_term <- "GO_TERM_NAME"   # 예시: GO_TERM_NAME 대신 실제 GO term 이름
# 또는
target_GO_ID <- "GO:0006631"       # 예시: 실제 GO ID

# 특정 GO term에 해당하는 유전자 리스트 추출
genes_in_target <- ego_hub@result$geneID[ego_hub@result$ID == target_GO_ID]
# 또는 ID를 사용하는 경우
# genes_in_target <- ego_hub@result$geneID[ego_hub@result$ID == target_GO_ID]

# 쉼표로 구분된 유전자 리스트를 분리
genes_in_target_list <- unlist(strsplit(genes_in_target, "/"))

# 결과 확인
print(genes_in_target_list)
