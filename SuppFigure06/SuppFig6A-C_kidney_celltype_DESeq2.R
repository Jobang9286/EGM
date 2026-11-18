# ============================================================
# Kidney 세포 타입별 DEG 분석 및 재생 마커 스크리닝
# 비교: HGPS_EGM (n=3) vs HGPS_control (n=3)
# ============================================================

library(InstaPrism)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# 작업 디렉토리 설정
setwd("EGM_project_for_deconvolution/integrated_analysis_for_lizard_factors/kidney")

# ------------------------------------------------------------
# 1. 데이터 로드
# ------------------------------------------------------------

# InstaPrism 결과 로드
deconv_result <- readRDS("kidney_deconv_result.rds")
Z_array <- readRDS("kidney_Z_array.rds")

# 재생 마커 정의
regen_markers <- data.frame(
  Category = c(
    # EMT_MET: 15개
    rep("EMT_MET", 15),
    # Blastema: 9개 (기존) + 7개 (신규, Krt19 제외) = 16개
    rep("Blastema", 15),
    # Injury_PT: 6개 (신규)
    rep("Injury_PT", 10),
    # ECM: 15개
    rep("ECM", 15),
    # Immune: 4개
    rep("Immune", 4),
    # Wnt/Shh: 15개
    rep("Wnt/Shh", 15),
    # Cytoskeletal: 2개
    rep("Cytoskeletal", 2),
    # ncRNA: 7개
    rep("ncRNA", 7)
  ),
  Gene = c(
    # EMT_MET (15)
    "Twist1", "Snai1", "Snai2", "Zeb1", "Zeb2", "Vim", "Fn1", "Cdh11", "Postn", 
    "Cdh1", "Krt5", "Krt8", "Krt14", "Krt17", "Krt18",
    
    # Blastema (16 = 9 기존 + 7 신규, Krt19 제외)
    "Prrx1", "Msx1", "Msx2", "Sall4", "Pax7", "Sox2", "Sox9", "Fgf8",  # 기존 9개
    "Six2", "Pax2", "Mki67", "Top2a", "Pcna", "Ccnb1", "Cenpf",  # 신규 7개
    
    # Injury_PT (6 신규)
    "Havcr1", "Lcn2", "Spp1", "Klf6", "Atf3","Cxcl1", 'Cxcl2', "Cxcl10", "Il1b","Ccn2",
    
    # ECM (15)
    "Col12a1", "Col20a1", "Col26a1", "Col11a1", "Col1a1", "Col3a1", "Col4a1", 
    "Col4a2", "Col5a1", "Col5a2", "Col6a1", "Col6a3", "Lama1", "Tnc", "Mmp2",
    
    # Immune (4)
    "Mrc1", "Ccl2", "Ccl7", "Cxcl12",
    
    # Wnt/Shh (15)
    "Wnt4", "Wnt9b", "Wnt7a", "Wnt3a", "Wnt2b", "Wnt5a", "Wnt5b", "Wnt6", 
    "Axin2", "Lef1", "Shh", "Ptch1", "Gli1", "Hoxa11", "Hoxd13",
    
    # Cytoskeletal (2)
    "Sptbn5", "Map1b",
    
    # ncRNA (7)
    "Snord26", "Snora63", "Scarna11", "Snora24", "Snord74", "Snord14", "Snu11"
  ),
  stringsAsFactors = FALSE
)
# ------------------------------------------------------------
# 2. 샘플 그룹 정의 (HGPS만 선택)
# ------------------------------------------------------------

sample_names <- dimnames(Z_array)[[1]]
all_sample_groups <- c(
  rep("WT_21M", 3),
  rep("WT_3M", 5),
  rep("HGPS_control", 3),
  rep("HGPS_EGM", 3)
)

# HGPS 샘플만 선택
hgps_idx <- which(all_sample_groups %in% c("HGPS_control", "HGPS_EGM"))
hgps_sample_names <- sample_names[hgps_idx]
hgps_groups <- all_sample_groups[hgps_idx]

# 결과 저장 디렉토리 생성
dir.create("deg_results", showWarnings = FALSE)
dir.create("deg_results/cell_type_specific", showWarnings = FALSE)
dir.create("deg_results/regeneration_markers", showWarnings = FALSE)

# ------------------------------------------------------------
# 3. 세포 타입별 DEG 분석
# ------------------------------------------------------------

cell_types <- dimnames(Z_array)[[3]]
deg_results_list <- list()
regen_results_list <- list()

for (ct in cell_types) {
  cat(paste0("분석 중: ", ct, "...\n"))
  
  tryCatch({
    # 발현 데이터 추출 (genes x samples)
    ct_expr <- t(Z_array[, , ct])
    ct_expr_hgps <- ct_expr[, hgps_idx]
    colnames(ct_expr_hgps) <- hgps_sample_names
    
    # 저발현 유전자 필터링
    keep_genes <- rowSums(ct_expr_hgps > 10) >= 3
    ct_expr_filtered <- ct_expr_hgps[keep_genes, ]
    
    if (nrow(ct_expr_filtered) < 100) {
      cat(paste0("  경고: 유전자 수 부족. 스킵.\n\n"))
      next
    }
    
    # DESeq2 분석
    coldata <- data.frame(
      condition = factor(hgps_groups, levels = c("HGPS_control", "HGPS_EGM")),
      row.names = colnames(ct_expr_filtered)
    )
    
    count_matrix <- round(ct_expr_filtered)
    dds <- DESeqDataSetFromMatrix(count_matrix, coldata, design = ~ condition)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("condition", "HGPS_EGM", "HGPS_control"))
    
    # 결과 저장
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      arrange(padj, desc(abs(log2FoldChange)))
    
    sig_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)
    sig_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)
    cat(paste0("  유의한 DEG: Up=", sig_up, ", Down=", sig_down, "\n"))
    
    deg_results_list[[ct]] <- res_df
    write.csv(res_df, 
              paste0("deg_results/cell_type_specific/", ct, "_DEG_HGPS_EGM_vs_control.csv"),
              row.names = FALSE)
    
    # 재생 마커 스크리닝
    regen_genes <- unique(regen_markers$Gene)
    regen_res <- res_df %>%
      filter(gene %in% regen_genes) %>%
      left_join(regen_markers, by = c("gene" = "Gene")) %>%
      arrange(padj)
    
    if (nrow(regen_res) > 0) {
      regen_results_list[[ct]] <- regen_res
      write.csv(regen_res,
                paste0("deg_results/regeneration_markers/", ct, "_regeneration_markers.csv"),
                row.names = FALSE)
      
      sig_regen <- sum(regen_res$padj < 0.05, na.rm = TRUE)
      cat(paste0("  재생 마커: ", nrow(regen_res), "개 검출, ", sig_regen, "개 유의\n"))
    }
    cat("\n")
    
  }, error = function(e) {
    cat(paste0("  오류: ", e$message, "\n\n"))
  })
}

# ------------------------------------------------------------
# 4. 통합 결과 요약
# ------------------------------------------------------------

# DEG 요약
all_deg_summary <- lapply(names(deg_results_list), function(ct) {
  res <- deg_results_list[[ct]]
  data.frame(
    CellType = ct,
    Total_genes = nrow(res),
    Sig_up = sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE),
    Sig_down = sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE),
    Sig_total = sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm = TRUE)
  )
}) %>% bind_rows() %>% arrange(desc(Sig_total))

write.csv(all_deg_summary, "deg_results/summary_DEG_by_celltype.csv", row.names = FALSE)

# 재생 마커 요약
all_regen_summary <- lapply(names(regen_results_list), function(ct) {
  res <- regen_results_list[[ct]]
  data.frame(
    CellType = ct,
    Total_markers = nrow(res),
    Sig_markers = sum(res$padj < 0.05, na.rm = TRUE),
    Sig_up = sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE),
    Sig_down = sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  )
}) %>% bind_rows() %>% arrange(desc(Sig_markers))

write.csv(all_regen_summary, "deg_results/summary_regeneration_markers_by_celltype.csv", row.names = FALSE)

# 카테고리별 재생 마커
all_regen_combined <- lapply(names(regen_results_list), function(ct) {
  res <- regen_results_list[[ct]]
  res$CellType <- ct
  res
}) %>% bind_rows()

category_summary <- all_regen_combined %>%
  filter(padj < 0.05) %>%
  group_by(Category) %>%
  summarise(
    Sig_markers = n(),
    Mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    Upregulated = sum(log2FoldChange > 0, na.rm = TRUE),
    Downregulated = sum(log2FoldChange < 0, na.rm = TRUE)
  ) %>%
  arrange(desc(Sig_markers))

write.csv(category_summary, "deg_results/summary_regeneration_markers_by_category.csv", row.names = FALSE)

# Top 30 재생 마커
top_regen_markers <- all_regen_combined %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(50) %>%
  select(gene, Category, CellType, log2FoldChange, padj)

write.csv(top_regen_markers, "deg_results/top30_regeneration_markers.csv", row.names = FALSE)

# ------------------------------------------------------------
# 5. 히트맵 생성
# ------------------------------------------------------------

# 유의한 재생 마커 행렬 생성
sig_regen <- all_regen_combined %>%
  filter(padj < 0.05) %>%
  select(gene, Category, CellType, log2FoldChange, padj)

log2fc_matrix <- sig_regen %>%
  select(gene, CellType, log2FoldChange) %>%
  pivot_wider(names_from = CellType, values_from = log2FoldChange) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

padj_matrix <- sig_regen %>%
  select(gene, CellType, padj) %>%
  pivot_wider(names_from = CellType, values_from = padj) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

log2fc_matrix[is.na(log2fc_matrix)] <- 0

# 카테고리별 정렬
gene_category <- sig_regen %>%
  select(gene, Category) %>%
  distinct()

gene_order <- gene_category %>%
  arrange(Category, gene) %>%
  pull(gene)

log2fc_matrix <- log2fc_matrix[gene_order, ]
padj_matrix <- padj_matrix[gene_order, ]

# 통계 유의성 표시 함수
add_significance <- function(j, i, x, y, width, height, fill) {
  padj_val <- padj_matrix[i, j]
  if (!is.na(padj_val)) {
    if (padj_val < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 8, col = "black", fontface = "bold"))
    } else if (padj_val < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 8, col = "black", fontface = "bold"))
    } else if (padj_val < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 8, col = "black", fontface = "bold"))
    }
  }
}

# 행 어노테이션
row_ha <- rowAnnotation(
  Category = gene_category$Category[match(rownames(log2fc_matrix), gene_category$gene)],
  col = list(Category = c(
    "EMT_MET" = "#E64B35FF",
    "ECM" = "#4DBBD5FF",
    "Wnt/Shh" = "#00A087FF",
    "Blastema" = "#3C5488FF",
    "Immune" = "#F39B7FFF",
    "Cytoskeletal" = "#8491B4FF",
    "Injury_PT" = "#91A54667",
    "ncRNA" = "#91D1C2FF"
  )),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  width = unit(0.5, "cm")
)

# 컬러 팔렛트
col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# 히트맵 생성
ht <- Heatmap(
  log2fc_matrix,
  name = "log2FC",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9),
  left_annotation = row_ha,
  rect_gp = gpar(col = "grey90", lwd = 0.3),
  cell_fun = add_significance,
  column_title = "Kidney: Regeneration Markers (HGPS_EGM vs HGPS_control)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "log2FC",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    at = c(-3, -2, -1, 0, 1, 2, 3),
    legend_height = unit(4, "cm")
  ),
  width = unit(6, "cm"),
  height = unit(12, "cm"),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 1)
)

# 히트맵 저장
pdf("deg_results/kidney_regeneration_markers_heatmap.pdf", width = 10, height = 12)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
     
dev.off()
