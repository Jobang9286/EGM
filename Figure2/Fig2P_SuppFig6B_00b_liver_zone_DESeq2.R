# ============================================================================
# Liver Cell Type-Specific DEG Analysis Using InstaPrism
# Date: 2024-11-16
# Purpose: Perform DEG analysis for HGPS EGM vs Control across liver zones
#          and screen for regeneration markers
# ============================================================================

library(InstaPrism)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridis)

# Set working directory
setwd("EGM_project_for_deconvolution/integrated_analysis_for_lizard_factors/liver")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading liver deconvolution results...\n")
# liver/results/InstaPrism/deconv_hgps_zonation.rds 와 동일 파일
liver_deconv <- readRDS("liver_deg_results/deconv_hgps_zonation.rds")

# Load bulk RNA-seq data
cat("Loading bulk RNA-seq data...\n")
# NOTE: bulk count는 liver/counts_liver.csv (comma-separated, gene symbol이 행 이름)
bulk_data <- read.csv("../../liver/counts_liver.csv",
                      header = TRUE, row.names = 1, check.names = FALSE)

# Check sample names
cat("\nSamples in deconvolution:\n")
print(colnames(liver_deconv@Post.ini.ct@theta))

cat("\nSamples in bulk data:\n")
print(colnames(bulk_data))

# Align sample names if needed
sample_names <- colnames(liver_deconv@Post.ini.ct@theta)
bulk_data <- bulk_data[, sample_names]

cat("\nData loaded successfully!\n")
cat("Genes:", nrow(bulk_data), "\n")
cat("Samples:", ncol(bulk_data), "\n")
cat("Cell types:", length(liver_deconv@map), "\n")

# ============================================================================
# 2. LOAD REGENERATION MARKERS
# ============================================================================

cat("\nLoading regeneration markers...\n")

# Regeneration markers data
regen_data <- "Category,Gene
EMT_MET,Twist1
EMT_MET,Snai1
EMT_MET,Snai2
EMT_MET,Zeb1
EMT_MET,Zeb2
EMT_MET,Vim
EMT_MET,Fn1
EMT_MET,Cdh11
EMT_MET,Postn
EMT_MET,Cdh1
EMT_MET,Krt5
EMT_MET,Krt8
EMT_MET,Krt14
EMT_MET,Krt17
EMT_MET,Krt18
EMT_MET,Krt19
Blastema,Prrx1
Blastema,Msx1
Blastema,Msx2
Blastema,Sall4
Blastema,Lix1
Blastema,Pax7
Blastema,Sox2
Blastema,Sox9
Blastema,Fgf8
ECM,Col12a1
ECM,Col20a1
ECM,Col26a
ECM,Col11a1
ECM,Col1a1
ECM,Col3a1
ECM,Col4a1
ECM,Col4a2
ECM,Col5a1
ECM,Col5a2
ECM,Col6a1
ECM,Col6a3
ECM,Lama1
ECM,Tnc
ECM,Mmp2
ECM,Mmp9
ECM,Timp1
ECM,Thbs1
ECM,Serpine1
ECM,Thbs2
Immune,Mrc1
Immune,Ccl2
Immune,Ccl7
Immune,Cxcl12
Wnt/Shh,Wnt4
Wnt/Shh,Wnt9b
Wnt/Shh,Wnt7a
Wnt/Shh,Wnt3a
Wnt/Shh,Wnt2b
Wnt/Shh,Wnt5a
Wnt/Shh,Wnt5b
Wnt/Shh,Wnt6
Wnt/Shh,Axin2
Wnt/Shh,Lef1
Wnt/Shh,Shh
Wnt/Shh,Ptch1
Wnt/Shh,Gli1
Wnt/Shh,HoxA11
Wnt/Shh,HoxD13
Cytoskeletal,Sptbn5
Cytoskeletal,Map1b
ncRNA,Snord26
ncRNA,Snora63
ncRNA,Scarna11
ncRNA,Snora24
ncRNA,Snord74
ncRNA,Snord14
ncRNA,Snu11"

regen_markers <- read.csv(text = regen_data)

cat("Marker categories:\n")
print(table(regen_markers$Category))

# ============================================================================
# 3. RECONSTRUCT CELL TYPE-SPECIFIC EXPRESSION
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("RECONSTRUCTING CELL TYPE-SPECIFIC EXPRESSION\n")
cat(rep("=", 70), "\n", sep = "")

cell_types <- names(liver_deconv@map)
Z_list <- list()

for (ct in cell_types) {
  cat("\nProcessing:", ct, "\n")
  Z_ct <- reconstruct_Z_ct_initial(
    InstaPrism_obj = liver_deconv,
    cell.type.of.interest = ct,
    n.core = 1,
    pseudo.min = 1E-8
  )
  Z_list[[ct]] <- Z_ct
  cat("  Dimensions:", dim(Z_ct)[1], "genes x", dim(Z_ct)[2], "samples\n")
}

cat("\nCell type-specific expression reconstructed!\n")

# ============================================================================
# 4. CREATE SAMPLE METADATA
# ============================================================================

cat("\nCreating sample metadata...\n")

# Sample groups based on file names
sample_info <- data.frame(
  sample = colnames(bulk_data),
  condition = ifelse(grepl("EGM", colnames(bulk_data)), "HGPS_EGM", "HGPS_control"),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample

cat("Sample grouping:\n")
print(table(sample_info$condition))

# ============================================================================
# 5. PERFORM DEG ANALYSIS FOR EACH CELL TYPE
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("PERFORMING DEG ANALYSIS\n")
cat(rep("=", 70), "\n", sep = "")

deg_results <- list()
dds_objects <- list()

for (ct in cell_types) {
  cat("\n", rep("-", 60), "\n", sep = "")
  cat("Analyzing:", ct, "\n")
  cat(rep("-", 60), "\n", sep = "")
  
  # Get cell type-specific expression
  counts_ct <- round(Z_list[[ct]])
  
  # Remove genes with all zeros
  counts_ct <- counts_ct[rowSums(counts_ct) > 0, ]
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts_ct,
    colData = sample_info,
    design = ~ condition
  )
  
  # Set reference level
  dds$condition <- relevel(dds$condition, ref = "HGPS_control")
  
  # Run DESeq2
  dds <- DESeq(dds, quiet = TRUE)
  dds_objects[[ct]] <- dds
  
  # Get results
  res <- results(dds, contrast = c("condition", "HGPS_EGM", "HGPS_control"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$pvalue), ]
  
  # Add significance column
  res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  
  deg_results[[ct]] <- res_df
  
  # Summary
  cat("\nDEG Summary:\n")
  cat("  Total genes tested:", nrow(res_df), "\n")
  cat("  Significant DEGs (padj<0.05, |log2FC|>1):", sum(res_df$significant, na.rm = TRUE), "\n")
  cat("  Upregulated (log2FC>1):", sum(res_df$significant & res_df$log2FoldChange > 1, na.rm = TRUE), "\n")
  cat("  Downregulated (log2FC<-1):", sum(res_df$significant & res_df$log2FoldChange < -1, na.rm = TRUE), "\n")
}

# ============================================================================
# 6. SAVE DEG RESULTS
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("SAVING RESULTS\n")
cat(rep("=", 70), "\n", sep = "")

# Create output directory
dir.create("liver_deg_results", showWarnings = FALSE, recursive = TRUE)

# Save individual DEG results
for (ct in cell_types) {
  ct_clean <- gsub("[^[:alnum:]]", "_", ct)
  outfile <- paste0("liver_deg_results/", ct_clean, "_deg.csv")
  write.csv(deg_results[[ct]], outfile, row.names = FALSE)
  cat("Saved:", outfile, "\n")
}

# Save combined results
combined_deg <- bind_rows(
  lapply(names(deg_results), function(ct) {
    deg_results[[ct]] %>%
      mutate(cell_type = ct) %>%
      select(cell_type, gene, everything())
  })
)
write.csv(combined_deg, "liver_deg_results/all_cell_types_deg.csv", row.names = FALSE)
cat("Saved: liver_deg_results/all_cell_types_deg.csv\n")

# Save RDS objects
saveRDS(deg_results, "liver_deg_results/deg_results_list.rds")
saveRDS(dds_objects, "liver_deg_results/dds_objects_list.rds")
saveRDS(Z_list, "liver_deg_results/Z_celltype_expression.rds")
cat("Saved: RDS objects\n")

# ============================================================================
# 7. REGENERATION MARKER SCREENING
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("REGENERATION MARKER SCREENING\n")
cat(rep("=", 70), "\n", sep = "")

# Screen for regeneration markers in each cell type
regen_screening <- list()

for (ct in cell_types) {
  cat("\n", rep("-", 60), "\n", sep = "")
  cat("Screening:", ct, "\n")
  cat(rep("-", 60), "\n", sep = "")
  
  # Get DEG results for this cell type
  deg_ct <- deg_results[[ct]]
  
  # Match with regeneration markers
  regen_overlap <- merge(
    regen_markers,
    deg_ct,
    by.x = "Gene",
    by.y = "gene",
    all.x = TRUE
  )
  
  # Add cell type column
  regen_overlap$cell_type <- ct
  
  # Count detected markers
  detected <- sum(!is.na(regen_overlap$log2FoldChange))
  significant <- sum(regen_overlap$significant, na.rm = TRUE)
  
  cat("  Total markers:", nrow(regen_overlap), "\n")
  cat("  Detected:", detected, "\n")
  cat("  Significant:", significant, "\n")
  
  # Summary by category
  cat("\n  By category:\n")
  category_summary <- regen_overlap %>%
    group_by(Category) %>%
    summarise(
      total = n(),
      detected = sum(!is.na(log2FoldChange)),
      significant = sum(significant, na.rm = TRUE),
      .groups = "drop"
    )
  print(category_summary)
  
  regen_screening[[ct]] <- regen_overlap
}

# Combine all cell types
regen_combined <- bind_rows(regen_screening)

# Save regeneration marker results
write.csv(regen_combined, "liver_deg_results/regeneration_markers_screening.csv", row.names = FALSE)
cat("\nSaved: liver_deg_results/regeneration_markers_screening.csv\n")

# ============================================================================
# 8. CREATE SUMMARY STATISTICS
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("CREATING SUMMARY STATISTICS\n")
cat(rep("=", 70), "\n", sep = "")

# DEG summary across cell types
deg_summary <- data.frame(
  cell_type = cell_types,
  total_genes = sapply(deg_results, nrow),
  significant_degs = sapply(deg_results, function(x) sum(x$significant, na.rm = TRUE)),
  upregulated = sapply(deg_results, function(x) sum(x$significant & x$log2FoldChange > 1, na.rm = TRUE)),
  downregulated = sapply(deg_results, function(x) sum(x$significant & x$log2FoldChange < -1, na.rm = TRUE))
)

write.csv(deg_summary, "liver_deg_results/deg_summary.csv", row.names = FALSE)
cat("Saved: liver_deg_results/deg_summary.csv\n")

print(deg_summary)

# Regeneration marker summary
regen_summary <- regen_combined %>%
  group_by(cell_type, Category) %>%
  summarise(
    total = n(),
    detected = sum(!is.na(log2FoldChange)),
    significant = sum(significant, na.rm = TRUE),
    upregulated = sum(significant & log2FoldChange > 1, na.rm = TRUE),
    downregulated = sum(significant & log2FoldChange < -1, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(regen_summary, "liver_deg_results/regeneration_summary.csv", row.names = FALSE)
cat("Saved: liver_deg_results/regeneration_summary.csv\n")

print(regen_summary)

cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE!\n")
cat(rep("=", 70), "\n\n")

cat("Output files created in: liver_deg_results/\n")
cat("  - Individual DEG results per cell type\n")
cat("  - all_cell_types_deg.csv (combined)\n")
cat("  - regeneration_markers_screening.csv\n")
cat("  - deg_summary.csv\n")
cat("  - regeneration_summary.csv\n")
cat("  - RDS objects for further analysis\n")
