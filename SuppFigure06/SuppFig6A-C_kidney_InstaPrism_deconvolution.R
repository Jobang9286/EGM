# 1. 라이브러리 로드
library(InstaPrism)
library(SingleCellExperiment)

# 2. 데이터 로드
sc_data <- readRDS("EGM_project_for_deconvolution/kidney/kidney_sc_data/Mousesub_sce.rds")
bulk_data <- read.table("EGM_project_for_deconvolution/kidney/kidney_bulk_data/counts_kidney.txt", 
                        header=TRUE, row.names=1)

# 3. 공통 유전자 추출
common_genes <- intersect(rownames(sc_data), rownames(bulk_data))
sc_expr <- as.matrix(assay(sc_data, "counts")[common_genes, ])
bulk_expr <- as.matrix(bulk_data[common_genes, ])
cell_type_labels <- as.character(colData(sc_data)$cellType)

# 4. Reference 생성
refPhi_kidney <- refPrepare(
  sc_Expr = sc_expr,
  cell.type.labels = cell_type_labels,
  cell.state.labels = cell_type_labels,
  pseudo.min = 1E-8
)

# 5. InstaPrism 실행
deconv_result <- InstaPrism(
  bulk_Expr = bulk_expr,
  refPhi_cs = refPhi_kidney,
  n.iter = 100,
  filter = TRUE,
  verbose = TRUE,
  convergence.plot = TRUE,
  n.core = 1
)

# 6. 결과 추출
cell_fractions <- t(deconv_result@Post.ini.ct@theta)
write.csv(cell_fractions, "EGM_project_for_deconvolution/kidney/fractions/kidney_cell_fractions1.csv")

cell_fractions

# 7. Cell-type specific 유전자 발현 (선택사항)
Z_array <- get_Z_array(deconv_result, resolution = 'ct', n.core = 1)
