# RRHO with P-BP vs C-BP on x, against two different y axes:
#   (a) fibroblast ageing, 24 vs 3 months (GSE304040 controls)
#   (b) EGM vs Control on x, P-BP vs C-BP on y  (axes swapped relative to (a))
#   (c) EGM vs Control on x, fibroblast ageing on y
# ---------------------------------------------------------------------------
# Quadrant convention: both rankings are NEGATED before RRHO2 sorts them, so
# up-up sits TOP-RIGHT (RRHO2's own default puts up-up bottom-left).
# With ageing on y, the *discordant* corners are the rejuvenation signal:
#   bottom-right = up in P-BP x down with age
#   top-left     = down in P-BP x up with age
#
# Expected:
#   P-BP x ageing (matched, 1780) : max 26.1 | up-up  56 | dn-dn  22 | disc 255 | rho -0.159
#   P-BP x EGM    (matched, 1780) : max 35.2 | up-up 310 | dn-dn 191 | disc  25 | rho +0.321
#   P-BP x EGM    (natural, 1948) : max 39.5 | up-up 317 | dn-dn 131 | disc  13 | rho +0.307
# The two y axes give mirror-image patterns: ageing loads the discordant corners,
# EGM loads the concordant ones.
suppressPackageStartupMessages({library(RRHO2); library(limma)})
OUT <- "RRHO_new/reproduced"

## ---- axis 1: P-BP vs CP -------------------------------------------------
## DDE = log2FC (r = 0.98 against our own DESeq2+apeglm lfc_raw); the file is a
## pre-filtered DE subset of 2,071 genes.
osk <- read.csv("RRHO_new/OSKM_vs_plasma_for_RRHO.csv"); names(osk)[1] <- "gene"

## ---- axis 2: EGM vs nEGM -------------------------------------------------
## Only needed to define the matched gene set, so the P-BP and EGM panels are
## comparable. Saved from the DESeq2 run on egm_blood/counts_symbol.csv
## (9 samples, ~group, contrast egm vs untreated, count>=10 in >=3 samples).
EGM <- readRDS(file.path(OUT, "DE_arms.rds"))$EGM

## ---- axis 3: ageing 24m vs 3m (GSE304040) --------------------------------
## Alle et al. OSK ear fibroblasts; untreated controls only, 4 vs 4, unpaired.
## The matrix is library-size-normalised counts, so limma on log2 rather than DESeq2.
X    <- read.csv(gzfile("GSE304040_3_mFB_OSK_onoff_30-_874248908_merged_featureCounts_norm.csv.gz"),
                 row.names = 1, check.names = FALSE)
s    <- colnames(X)
age  <- ifelse(grepl("-3m-", s), "m3", "m24")
cond <- ifelse(grepl("OSKonoff", s), "OSKonoff", ifelse(grepl("OSK", s), "OSK", "ctrl"))
L    <- log2(as.matrix(X[rowSums(X >= 10) >= 8, ]) + 1)
k    <- cond == "ctrl"
tt   <- topTable(eBayes(lmFit(L[, k], model.matrix(~ factor(age[k], c("m3","m24"))))),
                 coef = 2, number = Inf, sort.by = "none")
AG   <- data.frame(gene = rownames(tt), DDE = tt$logFC, padj = tt$adj.P.Val)

## ---- RRHO ---------------------------------------------------------------
gm <- Reduce(intersect, list(osk$gene, EGM$gene, AG$gene))     # matched set, 1780

build <- function(a, b, genes, labs) {
  g <- intersect(intersect(a$gene, b$gene), genes)
  n <- function(x) data.frame(gene = g, DDE = -x$DDE[match(g, x$gene)])  # up-up -> top-right
  set.seed(1)
  list(o = RRHO2_initialize(n(a), n(b), labels = labs, log10.ind = TRUE), g = g)
}
report <- function(r, a, b, lab) cat(sprintf(
  "%-30s common %4d | max %.1f | up-up %3d | dn-dn %3d | discord %3d | rho %+.3f\n", lab,
  length(r$g), max(r$o$hypermat, na.rm = TRUE),
  length(r$o$genelist_dd$gene_list_overlap_dd), length(r$o$genelist_uu$gene_list_overlap_uu),
  length(r$o$genelist_ud$gene_list_overlap_ud) + length(r$o$genelist_du$gene_list_overlap_du),
  cor(a$DDE[match(r$g, a$gene)], b$DDE[match(r$g, b$gene)], method = "spearman")))

## collapse RRHO2's 4-cell quadrant gap to a single thin one
thin <- function(H) {
  rNA <- which(apply(H, 1, function(r) all(is.na(r))))
  cNA <- which(apply(H, 2, function(c) all(is.na(c))))
  kR <- setdiff(seq_len(nrow(H)), rNA); kC <- setdiff(seq_len(ncol(H)), cNA)
  b1 <- sum(kR < min(rNA)); b2 <- sum(kC < min(cNA))
  ins <- function(M, at, m) {
    if (m == 1) rbind(M[1:at, , drop = FALSE], NA, M[(at+1):nrow(M), , drop = FALSE])
    else        cbind(M[, 1:at, drop = FALSE], NA, M[, (at+1):ncol(M), drop = FALSE])
  }
  list(M = ins(ins(H[kR, kC], b1, 1), b2, 2), b1 = b1, b2 = b2)
}

r_age  <- build(osk, AG,  gm,                             c("4FP vs CP", "ageing 24m vs 3m"))
# for the EGM panels EGM goes on x and P-BP on y, so the argument order flips
r_egm  <- build(EGM, osk, gm,                             c("EGM vs nEGM", "4FP vs CP"))
r_egmN <- build(EGM, osk, intersect(osk$gene, EGM$gene),  c("EGM vs nEGM", "4FP vs CP"))
report(r_age,  osk, AG,  "P-BP x ageing (matched)")
r_ega  <- build(EGM, AG, gm,                             c("EGM vs nEGM", "ageing 24m vs 3m"))
r_egaN <- build(EGM, AG, intersect(EGM$gene, AG$gene),   c("EGM vs nEGM", "ageing 24m vs 3m"))
report(r_egm,  EGM, osk, "EGM x P-BP (matched)")
report(r_egmN, EGM, osk, "EGM x P-BP (natural)")
report(r_ega,  EGM, AG,  "EGM x ageing (matched)")
report(r_egaN, EGM, AG,  "EGM x ageing (natural)")

## ---- figure --------------------------------------------------------------
draw_rrho <- function(M, b1, b2, xlab_l, xlab_r, ylab_b, ylab_t,
                      bar_title = expression("-log("*italic(P)*"-value)"),
                      cex_lab = 1.05, cex_bar = 0.95) {
  jet <- colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F",
                            "yellow","#FF7F00","red","#7F0000"))(101)
  nx <- nrow(M); ny <- ncol(M); xs <- seq_len(nx); ys <- seq_len(ny)
  par(fig = c(0.13, 0.99, 0.09, 0.80), mar = c(0,0,0,0), new = FALSE)
  image(xs, ys, M, col = jet, axes = FALSE, xlab = "", ylab = "",
        breaks = seq(min(M, na.rm = TRUE), max(M, na.rm = TRUE), length.out = 102))
  mtext(xlab_l, side = 1, at = mean(xs[1:b2]),      line = 0.9, cex = cex_lab)
  mtext(xlab_r, side = 1, at = mean(xs[(b2+2):nx]), line = 0.9, cex = cex_lab)
  mtext(ylab_b, side = 2, at = mean(ys[1:b1]),      line = 0.9, cex = cex_lab, las = 0)
  mtext(ylab_t, side = 2, at = mean(ys[(b1+2):ny]), line = 0.9, cex = cex_lab, las = 0)
  # horizontal colour bar, top right, high value on the LEFT
  par(fig = c(0.55, 0.99, 0.855, 0.925), mar = c(0,0,0,0), new = TRUE)
  plot.new(); plot.window(xlim = c(0,1), ylim = c(0,1))
  rect(seq(0, 1-1/101, length.out = 101), 0, seq(1/101, 1, length.out = 101), 1,
       col = rev(jet), border = NA)
  mtext(bar_title, side = 3, line = 0.25, cex = cex_bar)
}
save2 <- function(r, xl, xr, yb, yt, stem) {
  th <- thin(r$o$hypermat)
  png(file.path(OUT, paste0(stem, ".png")), width = 1150, height = 1350, res = 300)
  draw_rrho(th$M, th$b1, th$b2, xl, xr, yb, yt); dev.off()
  pdf(file.path(OUT, paste0(stem, ".pdf")), width = 3.8, height = 4.5)
  draw_rrho(th$M, th$b1, th$b2, xl, xr, yb, yt); dev.off()
  cat("wrote", stem, "\n")
}
save2(r_age,  "C-BP",    "P-BP", "3 month", "24 month", "RRHO_4FP_x_ageing_styled")
save2(r_egm,  "Control", "EGM",   "C-BP",    "P-BP",    "RRHO_EGM_x_4FP_styled_matched")
save2(r_egmN, "Control", "EGM",   "C-BP",    "P-BP",    "RRHO_EGM_x_4FP_styled")
save2(r_ega,  "Control", "EGM",   "3 month", "24 month", "RRHO_EGM_x_ageing_styled_matched")
save2(r_egaN, "Control", "EGM",   "3 month", "24 month", "RRHO_EGM_x_ageing_styled")
