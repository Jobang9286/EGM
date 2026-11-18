# ------------------------------------------------------------------------------
# Venn diagrams (Figure E) + RRHO heatmaps (Figure G): EGM vs young fibroblast / young skin
# ------------------------------------------------------------------------------
library(RRHO2)

dir_root <- "reproduce"
dir_out  <- file.path(dir_root, "out")

rank_tbl <- function(f) {                       # DDE = log2FoldChange (NOT signed -log10 p)
  x <- read.csv(file.path(dir_root, "data", f))
  x <- x[!duplicated(x[, 1]), ]
  rownames(x) <- x[, 1]
  x
}

# in vitro fibroblast, EGM vs no-EGM
egm  <- rank_tbl("egm_vs_no_egm_for_RRHO.csv")
# Mahmoudi et al. 2018 ear fibroblast, 3 mo vs 29 mo
fb   <- rank_tbl("fb_young_vs_old_for_RRHO.csv")
# GSE132040 (Tabula Muris Senis bulk) skin, 1+3 mo vs 21+24+27 mo
skin <- rank_tbl("skin_young_vs_old_for_RRHO.csv")

rrho <- function(a, b, la, lb) {
  cg <- intersect(rownames(a), rownames(b))     # fb: 1803 genes, skin: 820 genes
  RRHO2_initialize(as.data.frame(a[cg, ]), as.data.frame(b[cg, ]),
                   labels = c(la, lb), log10.ind = TRUE)
}

o_fb <- rrho(fb,   egm, "Young fibroblast", "EGM")   # max -log10P = 11.62 -> key 12
o_sk <- rrho(skin, egm, "Young skin",       "EGM")   # max -log10P = 13.00 -> key 13

## ---- Figure G ----------------------------------------------------------------
pdf(file.path(dir_out, "figG_right_fibroblast_x_EGM.pdf"), 6, 6); RRHO2_heatmap(o_fb); dev.off()
pdf(file.path(dir_out, "figG_left_skin_x_EGM.pdf"),        6, 6); RRHO2_heatmap(o_sk); dev.off()

## ---- Figure E ----------------------------------------------------------------
pdf(file.path(dir_out, "figE_left_venn_EGM_x_fibroblast.pdf"), 6, 4); RRHO2_vennDiagram(o_fb, "uu"); dev.off()
pdf(file.path(dir_out, "figE_right_venn_EGM_x_skin.pdf"),      6, 4); RRHO2_vennDiagram(o_sk, "uu"); dev.off()

## ---- checks against the saved reference gene lists --------------------------
chk <- function(computed, saved) setequal(computed, saved)
stopifnot(
  chk(o_fb$genelist_uu$gene_list_overlap_uu, read.csv(file.path(dir_root,"ref","genelist_ear_uu.list.csv"))[,2]),  # 65
  chk(o_fb$genelist_dd$gene_list_overlap_dd, read.csv(file.path(dir_root,"ref","genelist_ear_dd.list.csv"))[,2]),  # 95
  chk(o_fb$genelist_du$gene_list_overlap_du, read.csv(file.path(dir_root,"ref","genelist_ear_du.list.csv"))[,2]),  # 57
  chk(o_sk$genelist_uu$gene_list_overlap_uu, readLines(file.path(dir_root,"ref","uu.list.csv"))),                  # 52
  length(o_fb$genelist_uu$gene_list2_uu) == 302,   # Figure E left:  EGM
  length(o_fb$genelist_uu$gene_list1_uu) == 173,   #                 Young fibroblast
  length(o_sk$genelist_uu$gene_list1_uu) == 204,   # Figure E right: Young skin
  length(o_sk$genelist_uu$gene_list2_uu) ==  88    #                 EGM
)

write.csv(data.frame(gene = o_fb$genelist_uu$gene_list_overlap_uu),
          file.path(dir_out, "overlap_EGM_youngFibroblast_65.csv"), row.names = FALSE)
write.csv(data.frame(gene = o_sk$genelist_uu$gene_list_overlap_uu),
          file.path(dir_out, "overlap_EGM_youngSkin_52.csv"), row.names = FALSE)
