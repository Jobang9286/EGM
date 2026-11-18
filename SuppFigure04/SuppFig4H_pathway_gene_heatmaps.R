library(pheatmap)

## ---- rlog (setup) -----------------------------------------------------------
rld <- rlogTransformation(dds)
rld$condition

rld

rld_mat <- assay(rld)
rld_mat

## ---- MAPK signaling pathway (n_MAPK.pdf) ------------------------------------
genes <- c("Flt4", "Hspb1", "Cacna1c", "Il1rap", "Dusp16", "Egfr", "Rap1b", "Rps6ka4", "Fgf7", "Pak1", "Dusp10", "Rac3", "Flna", "Map3k8", "Cd14", "Map3k6", "Map4k4", "Pdgfra", "Map2k1", "Tgfb1", "Il1r1", "Hgf", "Vegfc", "Pla2g4a", "Ngf", "Dusp9", "Mapk8ip3", "Dusp6", "Nfkb1", "Vegfa", "Mapk13", "Ereg", "Nfkb2", "Nr4a1", "Il1a", "Cacnb3", "Gdnf", "Rasa1", "Il1b", "Mapkapk2", "Fgfr1")
mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("condition")])
pheatmap(mat, annotation_col = anno, border_color = "white", cluster_cols = F)

## ---- glycolysis (n_Gly.pdf) -------------------------------------------------
genes <- c("Tpi1", "Pgam1", "Eno1", "Bpgm", "Hk2", "Hk1", "Gpi1", "Hk3", "Pfkl", "Ldha", "Pkm", "Pgk1", "Aldoa", "Gapdh", "Pgm1")
mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("condition")])
pheatmap(mat, annotation_col = anno, border_color = "white", cluster_cols = F)

## ---- PI3K-Akt signaling pathway (n_pi3k.pdf) --------------------------------
genes <- c("Csf3", "Lamc3", "Pkn3", "Itgb3", "Flt4", "Tnc", "Il4ra", "Thbs2", "Il2rg", "Thbs1", "Egfr", "Gys1", "Ccnd3", "Fgf7", "Ppp2r5e", "Itgav", "Jak2", "Ywhah", "Magi1", "Ifnar2", "Pdgfra", "Map2k1", "Hsp90aa1", "Hgf", "F2r", "Vegfc", "Ngf", "Osmr", "Gng11", "Nfkb1", "Vegfa", "Ereg", "Nr4a1", "Il6", "Gdnf", "Itga7", "Itga5", "Gnb5", "Tlr4", "Tlr2", "Fgfr1")
mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("condition")])
pheatmap(mat, annotation_col = anno, border_color = "white", cluster_cols = F)

## ---- lysosome (n_lyso.pdf) --------------------------------------------------
genes <- c("Scarb2", "Ctsa", "Hexb", "Gaa", "Fuca2", "Hexa", "Laptm5", "Naga", "Lipa", "Ctss", "Gga2", "Cln3", "Ap1g2", "Ids", "Man2b1", "Cd68", "Atp6v0d1", "Ctsd", "Lgmn")
mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("condition")])
pheatmap(mat, annotation_col = anno, border_color = "white", cluster_cols = F)

## ---- p53 & senescence (p53_senescence_heatmap.pdf) --------------------------
genes <- c("Tgfb2", "Cdkn1a", "Smad3", "Cdkn2a", "Pik3r3", "Itpr3", "Pik3r1", "Foxo3", "Zfp36l2", "Tgfbr1", "Tgfbr2", "Rbl2", "Ppp3ca", "Ccnb2", "Ccnd2", "Ccnd1", "Rbbp4", "Raf1", "Rrm2", "Aifm2", "Sesn1", "Sesn2", "Igf1")
mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("condition")])
pheatmap(mat, annotation_col = anno, border_color = "white", cluster_cols = F)

## ---- mTOR signaling pathway (n_mtor.pdf) ------------------------------------
genes <- c("Atp6v1g2", "Irs1", "Nprl3", "Pik3r3", "Igf1", "Pik3r1", "Slc7a5", "Stk11", "Sesn2", "Rps6ka1", "Atp6v1b2", "Grb10", "Lamtor1", "Raf1")
mat <- assay(rld)[genes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[c("condition")])
pheatmap(mat, annotation_col = anno, border_color = "white", cluster_cols = F)
