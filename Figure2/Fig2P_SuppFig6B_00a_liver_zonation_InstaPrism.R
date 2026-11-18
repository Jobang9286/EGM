# ============================================================================
# Liver EGM Project - Complete Analysis Workflow
# From Data Loading to Deconvolution and DEG Analysis
# Author: Analysis Pipeline
# Date: 2024
# ============================================================================

# ============================================================================
# Part 1: Library Loading and Setup
# ============================================================================

library(InstaPrism)
library(SingleCellExperiment)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Set working directory
setwd("EGM_project_for_deconvolution/liver")

# 출력 폴더 미리 생성 (없으면 write.csv/saveRDS 단계에서 에러)
for (d in c("fractions", "results/InstaPrism", "deg_results/zone_specific", "go_results", "new_result")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ============================================================================
# Part 2: Data Loading
# ============================================================================

cat("Step 1: Loading data...\n")

# Single cell reference data
# hepatocyte_zonation_sce.rds: GSE84498 간세포, colData$zone = Periportal/Midzonal/Pericentral
sc_data <- readRDS("liver_sc_data/hepatocyte_zonation_sce.rds")

# Bulk RNA-seq data
# NOTE: bulk count는 liver/counts_liver.csv (comma-separated)
bulk_data <- read.csv("counts_liver.csv", 
                      header = TRUE, 
                      row.names = 1, 
                      check.names = FALSE)

cat("  - SC data:", ncol(sc_data), "cells,", nrow(sc_data), "genes\n")
cat("  - Bulk data:", ncol(bulk_data), "samples,", nrow(bulk_data), "genes\n")

# ============================================================================
# Part 3: Data Preprocessing
# ============================================================================

cat("\nStep 2: Data preprocessing...\n")

# Extract common genes
common_genes <- intersect(rownames(sc_data), rownames(bulk_data))
cat("  - Common genes:", length(common_genes), "\n")

# Prepare single cell expression matrix
sc_expr <- as.matrix(assay(sc_data, "counts")[common_genes, ])

# Prepare bulk expression matrix
bulk_expr <- as.matrix(bulk_data[common_genes, ])

# Extract cell type labels (zonation info)
cell_type_labels <- as.character(colData(sc_data)$zone)  # colData 컬럼명은 "zone"

# Check cell type distribution
cat("  - Cell type distribution:\n")
print(table(cell_type_labels))

# ============================================================================
# Part 4: Create InstaPrism Reference
# ============================================================================

cat("\nStep 3: Creating InstaPrism reference...\n")

refPhi_liver <- refPrepare(
  sc_Expr = sc_expr,
  cell.type.labels = cell_type_labels,
  cell.state.labels = cell_type_labels,  # Using zonation as cell states
  pseudo.min = 1E-8
)

cat("  - Reference created with", length(refPhi_liver@map), "cell types\n")

# ============================================================================
# Part 5: Run InstaPrism Deconvolution
# ============================================================================

cat("\nStep 4: Running InstaPrism deconvolution...\n")

# Subset bulk data for HGPS samples (Control vs EGM)
hgps_samples <- grep("HGPS", colnames(bulk_expr), value = TRUE)
bulk_hgps <- bulk_expr[, hgps_samples]

cat("  - HGPS samples:", ncol(bulk_hgps), "\n")
cat("  - Sample names:\n")
print(colnames(bulk_hgps))

# Run deconvolution
deconv_hgps_zonation <- InstaPrism(
  bulk_Expr = bulk_hgps,
  refPhi_cs = refPhi_liver,
  n.iter = 100,
  filter = TRUE,
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
  pseudo.min = 1E-8,
  verbose = TRUE,
  convergence.plot = TRUE,
  n.core = 1
)

cat("\n  ✓ Deconvolution completed!\n")

# ============================================================================
# Part 6: Extract Cell Type Fractions
# ============================================================================

cat("\nStep 5: Extracting cell type fractions...\n")

# Extract cell type fractions (samples x cell types)
cell_fractions <- t(deconv_hgps_zonation@Post.ini.ct@theta)

# Add sample group information
sample_groups <- ifelse(grepl("control", colnames(bulk_hgps)), 
                        "HGPS_control", 
                        "HGPS_EGM")

cell_fractions_df <- data.frame(
  Sample = rownames(cell_fractions),
  Group = sample_groups,
  cell_fractions,
  check.names = FALSE
)

# Save fractions
write.csv(cell_fractions_df, 
          'fractions/liver_zonation_cell_fractions.csv',
          row.names = FALSE)

cat("  ✓ Cell fractions saved\n")
print(head(cell_fractions_df))

# ============================================================================
# Part 7: Extract Cell-Type Specific Expression (Z array)
# ============================================================================

cat("\nStep 6: Extracting cell-type specific expression...\n")

# Get 3D array: samples x genes x cell types
Z_array <- get_Z_array(deconv_hgps_zonation, resolution = 'ct', n.core = 1)

cat("  - Z array dimensions:", paste(dim(Z_array), collapse=" x "), "\n")
cat("  - (samples x genes x cell_types)\n")

# Save Z array
saveRDS(Z_array, 'results/InstaPrism/liver_Z_array.rds')

# saveRDS(deconv_hgps_zonation, "results/InstaPrism/deconv_hgps_zonation.rds")   # 00b reads this object
