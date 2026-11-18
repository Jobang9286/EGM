## ============ 1. EGM vs untreated (egm_blood, DESeq2) ============

suppressPackageStartupMessages(library(DESeq2))
setwd("4FP_plasma")

C  <- read.csv("egm_blood/counts_symbol.csv", row.names = 1, check.names = FALSE)
mp <- read.delim("egm_blood/mapfile.tsv")
mp$group <- factor(mp$TREATMENT, levels = c("untreated","egm","eno1","gclc","map1s"))
rownames(mp) <- mp$sample
C <- as.matrix(C[, mp$sample])
cat("counts:", nrow(C), "genes x", ncol(C), "samples\n")
print(table(mp$group))

keep <- rowSums(C >= 10) >= 3
C <- C[keep, ]
cat("after filter (count>=10 in >=3 samples):", nrow(C), "genes\n")

# one model over all 9 samples: the n=1 arms borrow dispersion from the
# replicated GFP and EGM groups, which is the only way to get usable stats there
dds <- DESeqDataSetFromMatrix(round(C), mp, ~ group)
dds <- DESeq(dds, quiet = TRUE)
cat("\ncoefficients:", paste(resultsNames(dds), collapse = ", "), "\n")

## per-arm DE vs untreated, saved for RRHO.
OUT <- "RRHO_new/reproduced"; dir.create(OUT, showWarnings = FALSE)

arms <- c(EGM = "egm", Eno1 = "eno1", Gclc = "gclc", Map1s = "map1s")
DE <- list()
for (nm in names(arms)) {
  r <- results(dds, contrast = c("group", arms[[nm]], "untreated"))
  d <- data.frame(gene = rownames(r), DDE = r$log2FoldChange,
                  pvalue = r$pvalue, padj = r$padj)
  d <- d[!is.na(d$DDE) & !is.na(d$pvalue), ]
  DE[[nm]] <- d
  write.csv(d[, c("gene","DDE")], file.path(OUT, sprintf("%s_vs_nEGM_for_RRHO.csv", nm)),
            row.names = FALSE)
  cat(sprintf("%-6s vs nEGM : %5d genes tested | padj<0.05 %5d | |lfc|>1 %5d | DDE %.2f .. %.2f\n",
              nm, nrow(d), sum(d$padj < .05, na.rm = TRUE),
              sum(abs(d$DDE) > 1), min(d$DDE), max(d$DDE)))
}
saveRDS(DE, file.path(OUT, "DE_arms.rds"))

## ============ 2. P-BP vs C-BP axis ============

osk  <- read.csv("RRHO_new/OSKM_vs_plasma_for_RRHO.csv"); names(osk)[1] <- "gene"

## ============ 3. partial OSKM, GSE231853 (WT and Ercc1 d/-; D2, D4 vs D0) ============

suppressPackageStartupMessages({library(DESeq2)})
cat("org.Mm.eg.db:", requireNamespace("org.Mm.eg.db", quietly = TRUE), "\n")

raw <- read.delim(gzfile("analysis/geo/GSE231853_counts_TOR.tsv.gz"), check.names = FALSE)
cm  <- as.matrix(raw[, 7:ncol(raw)]); rownames(cm) <- raw$Geneid
cat("counts:", nrow(cm), "genes x", ncol(cm), "samples\n")

# metadata from the series matrix titles: TOR_<n>_<genotype>_<treatment>
sm <- readLines("analysis/geo/GSE231853_sm.txt", warn = FALSE)
ttl <- unlist(strsplit(grep("^!Sample_title", sm, value = TRUE)[1], "\t"))[-1]
ttl <- gsub('"', '', ttl); ttl <- ttl[grepl("^TOR_", ttl)]
md <- data.frame(title = ttl,
                 id    = paste0("TOR_", sub("^TOR_(\\d+)_.*$", "\\1", ttl)),
                 geno  = ifelse(grepl("_wt", ttl), "WT", "ERCC1"),
                 vc    = grepl("_VC$", ttl),
                 day   = ifelse(grepl("control", ttl), "D0",
                         ifelse(grepl("day-2", ttl), "D2", "D4")))
md <- md[md$id %in% colnames(cm), ]
sel <- md[md$geno == "WT" & !md$vc, ]
sel$day <- factor(sel$day, c("D0","D2","D4"))
print(table(sel$day))
cm <- cm[, sel$id]
cat("\nlibrary sizes (M):", paste(round(colSums(cm)/1e6, 1), collapse = ", "), "\n")

## symbol mapping and DESeq2, WT D2 / D4 vs D0.
suppressPackageStartupMessages({library(org.Mm.eg.db); library(AnnotationDbi)})
sym <- AnnotationDbi::mapIds(org.Mm.eg.db, sub("\\..*$", "", rownames(cm)),
                             "SYMBOL", "ENSEMBL", multiVals = "first")
ok <- !is.na(sym)
cs <- rowsum(cm[ok, ], sym[ok])
cs <- cs[rowSums(cs >= 10) >= 3, ]
cat("after symbol mapping + filter:", nrow(cs), "genes\n")

dds2 <- DESeq(DESeqDataSetFromMatrix(round(cs), data.frame(day = sel$day), ~ day), quiet = TRUE)
POR <- list()
for (d in c("D2","D4")) {
  r <- results(dds2, contrast = c("day", d, "D0"))
  x <- data.frame(gene = rownames(r), DDE = r$log2FoldChange, padj = r$padj)
  x <- x[!is.na(x$DDE) & !is.na(r$pvalue), ]
  POR[[d]] <- x
  cat(sprintf("partial OSKM %s vs D0 : %5d genes | padj<0.05 %5d | up %4d down %4d\n",
              d, nrow(x), sum(x$padj < .05, na.rm = TRUE),
              sum(x$padj < .05 & x$DDE > 0, na.rm = TRUE),
              sum(x$padj < .05 & x$DDE < 0, na.rm = TRUE)))
}
cat("\nsanity - reprogramming factors and pluripotency markers (D4 log2FC):\n")
for (g in c("Pou5f1","Sox2","Klf4","Myc","Nanog","Lin28a","Dppa5a","Col1a1","Thy1"))
  if (g %in% POR$D4$gene) cat(sprintf("  %-8s %+6.2f\n", g, POR$D4$DDE[POR$D4$gene == g]))

## Ercc1 d/-: D2 / D4 vs D0.
suppressPackageStartupMessages({library(DESeq2)})
selE <- md[md$geno == "ERCC1" & !md$vc, ]
selE$day <- factor(selE$day, c("D0","D2","D4"))
print(table(selE$day))

cmE <- as.matrix(raw[, 7:ncol(raw)]); rownames(cmE) <- raw$Geneid
cmE <- cmE[, selE$id]
symE <- AnnotationDbi::mapIds(org.Mm.eg.db, sub("\\..*$", "", rownames(cmE)),
                              "SYMBOL", "ENSEMBL", multiVals = "first")
okE <- !is.na(symE)
csE <- rowsum(cmE[okE, ], symE[okE]); csE <- csE[rowSums(csE >= 10) >= 3, ]
cat("genes:", nrow(csE), "\n")

ddsE <- DESeq(DESeqDataSetFromMatrix(round(csE), data.frame(day = selE$day), ~ day), quiet = TRUE)
PORE <- list()
for (d in c("D2","D4")) {
  r <- results(ddsE, contrast = c("day", d, "D0"))
  x <- data.frame(gene = rownames(r), DDE = r$log2FoldChange, padj = r$padj)
  PORE[[d]] <- x[!is.na(x$DDE) & !is.na(r$pvalue), ]
  cat(sprintf("ERCC1 partial OSKM %s vs D0 : %5d genes | padj<0.05 %5d\n",
              d, nrow(PORE[[d]]), sum(PORE[[d]]$padj < .05, na.rm = TRUE)))
}
cat("\nsanity (ERCC1 D4 log2FC):\n")
for (g in c("Pou5f1","Sox2","Myc","Col1a1","Cdkn2a","Cdkn1a"))
  if (g %in% PORE$D4$gene) cat(sprintf("  %-8s %+6.2f\n", g, PORE$D4$DDE[PORE$D4$gene == g]))

## ============ 4. OSK / OSK on-off / ageing, GSE304040 (limma, animal-paired) ============

cat("limma:", requireNamespace("limma", quietly = TRUE), "\n")
X <- read.csv(gzfile("GSE304040_3_mFB_OSK_onoff_30-_874248908_merged_featureCounts_norm.csv.gz"),
              row.names = 1, check.names = FALSE)
cat("matrix:", nrow(X), "genes x", ncol(X), "samples\n")

s  <- colnames(X)
age <- ifelse(grepl("-3m-", s), "m3", "m24")
cond <- ifelse(grepl("OSKonoff", s), "OSKonoff", ifelse(grepl("OSK", s), "OSK", "ctrl"))
animal <- sub("^.*-", "", s)
meta <- data.frame(sample = s, age = age, cond = factor(cond, c("ctrl","OSK","OSKonoff")),
                   animal = animal)
print(table(meta$age, meta$cond))
cat("\nanimals per age:\n")
print(tapply(meta$animal, meta$age, function(x) paste(sort(unique(x)), collapse=", ")))
cat("\npaired? each animal appears in all three conditions:\n")
print(table(meta$animal, meta$cond))

## limma per age, animal-paired; ageing axis from the controls.
suppressPackageStartupMessages(library(limma))
keep <- rowSums(X >= 10) >= 8
L <- log2(as.matrix(X[keep, ]) + 1)
cat("genes after filter:", nrow(L), "\n\n")

OSKD <- list()
for (a in c("m3","m24")) {
  k <- meta$age == a
  d <- droplevels(meta[k, ])
  des <- model.matrix(~ animal + cond, data = d)
  fit <- eBayes(lmFit(L[, k], des))
  for (cn in c("condOSK","condOSKonoff")) {
    tt <- topTable(fit, coef = cn, number = Inf, sort.by = "none")
    lab <- paste0(a, "_", sub("^cond", "", cn))
    OSKD[[lab]] <- data.frame(gene = rownames(tt), DDE = tt$logFC, padj = tt$adj.P.Val)
    cat(sprintf("%-14s : adj.P<0.05 %5d   |logFC|>1 %4d\n", lab,
                sum(tt$adj.P.Val < .05), sum(abs(tt$logFC) > 1)))
  }
}
# ageing axis, unpaired 4 vs 4 controls
k <- meta$cond == "ctrl"
des <- model.matrix(~ factor(meta$age[k], c("m3","m24")))
tt <- topTable(eBayes(lmFit(L[, k], des)), coef = 2, number = Inf, sort.by = "none")
OSKD[["ageing_24m_vs_3m"]] <- data.frame(gene = rownames(tt), DDE = tt$logFC, padj = tt$adj.P.Val)
cat(sprintf("%-14s : adj.P<0.05 %5d\n", "ageing", sum(tt$adj.P.Val < .05)))

cat("\nsanity - OSK induction (m24 OSK logFC):\n")
for (g in c("Pou5f1","Sox2","Klf4","Myc","Gsta4","Col1a1"))
  if (g %in% OSKD$m24_OSK$gene) cat(sprintf("  %-8s %+6.2f\n", g, OSKD$m24_OSK$DDE[OSKD$m24_OSK$gene == g]))
cat("\n  ...and after 4-day withdrawal (m24 OSKonoff):\n")
for (g in c("Pou5f1","Sox2","Klf4","Gsta4"))
  if (g %in% OSKD$m24_OSKonoff$gene) cat(sprintf("  %-8s %+6.2f\n", g, OSKD$m24_OSKonoff$DDE[OSKD$m24_OSKonoff$gene == g]))

## ============ 5. gene sets: GO BP (gmt) + MSigDB Hallmark ============

## GO BP sets.
gmt <- readLines("chromatin/rank/go_bp_mouse.gmt")
sets <- lapply(strsplit(gmt, "\t"), function(x) x[-(1:2)])
names(sets) <- sapply(strsplit(gmt, "\t"), `[`, 1)
sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
cat(sprintf("GO BP sets: %d\n", length(sets)))

## GO id -> term name helper.
if (!requireNamespace("GO.db", quietly = TRUE))
  BiocManager::install("GO.db", ask = FALSE, update = FALSE)
suppressPackageStartupMessages(library(GO.db))
nm <- function(ids) sapply(ids, function(i) {
  x <- tryCatch(AnnotationDbi::Term(GOTERM[[i]]), error = function(e) NA); if (is.null(x)) NA else x })

## Hallmark sets.
suppressPackageStartupMessages(library(msigdbr))
mh <- msigdbr(species = "Mus musculus", collection = "H")
HS <- split(mh$gene_symbol, mh$gs_name)
cat("Hallmark sets:", length(HS), "\n")

## ============ 6. scan every set on every axis ============

ALLS <- c(sets, setNames(HS, base::sub("^HALLMARK_", "H:", names(HS))))
ALLS <- ALLS[lengths(ALLS) >= 12 & lengths(ALLS) <= 400]
cat("gene sets scanned:", length(ALLS), "\n")

AXS <- list(FP4 = osk, EGM = DE$EGM,
            OSKM_WT = POR$D4, OSKM_ERCC1 = PORE$D4,
            OSK24 = OSKD$m24_OSK, OSKonoff24 = OSKD$m24_OSKonoff,
            ageing = OSKD$ageing_24m_vs_3m)
V   <- lapply(AXS, function(x) setNames(x$DDE, x$gene))
BGm <- sapply(V, median)

# step 1: centred medians for every set on every axis (fast)
MED <- sapply(names(V), function(k) {
  v <- V[[k]]
  sapply(ALLS, function(S) { i <- v[intersect(S, names(v))]
    if (length(i) < 8) NA_real_ else median(i) - BGm[[k]] })
})
MED <- MED[complete.cases(MED), , drop = FALSE]
cat("sets with >=8 genes on every axis:", nrow(MED), "\n")

# step 2: keep sets where treatments AND reprogramming agree, and ageing opposes
tx   <- sign(MED[,"FP4"]) == sign(MED[,"EGM"])
rp   <- rowSums(sign(MED[, c("OSKM_WT","OSKM_ERCC1","OSK24","OSKonoff24")]) == sign(MED[,"EGM"]))
opp  <- sign(MED[,"ageing"]) == -sign(MED[,"EGM"])
cand <- which(tx & opp & rp >= 3)
cat("candidate sets (P-BP=EGM, >=3/4 reprogramming agree, ageing opposite):", length(cand), "\n")

## set-name helper.
nmx <- function(s) if (startsWith(s,"H:")) base::sub("^H:","Hallmark ",s) else nm(s)

## sign filter + Wilcoxon P for every candidate set (M2 / P2).
RPB <- c("OSKM_WT","OSKM_ERCC1","OSK24")     # active reprogramming
# permissive: ageing opposite & significant; both treatments same sign, >=1 significant;
# >=2 of the active reprogramming axes same sign, >=1 significant
tx  <- sign(MED[,"FP4"]) == sign(MED[,"EGM"])
opp <- sign(MED[,"ageing"]) == -sign(MED[,"EGM"])
c2  <- rownames(MED)[tx & opp]
P2 <- t(sapply(c2, function(s) sapply(names(V), function(k) {
  v <- V[[k]]; suppressWarnings(wilcox.test(v[intersect(ALLS[[s]], names(v))], v)$p.value) })))
colnames(P2) <- names(V); M2 <- MED[c2, , drop = FALSE]

## axis definitions (extra, core, lab6).
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
extra <- list(OSK3 = OSKD$m3_OSK, ONOFF3 = OSKD$m3_OSKonoff)
core <- c("FP4","EGM","OSK24","OSKonoff24","ageing")
lab6 <- c(FP4="P-BP", EGM="EGM",
          OSK3="OSK 3m", ONOFF3="OSK on-off 3m",
          OSK24="OSK 24m", OSKonoff24="OSK on-off 24m", ageing="Ageing")

## dot-size range.
RNG <- c(1.0, 4.6)

## ============ 7. relaxed scan and curation of the programmes ============

## relaxed rule (ageing P < 0.10); down-direction candidates, ranked by ageing P.
rel <- sign(M2[,"FP4"]) == sign(M2[,"EGM"]) &
       sign(M2[,"OSK24"]) == sign(M2[,"EGM"]) &
       sign(M2[,"ageing"]) == -sign(M2[,"EGM"]) &
       P2[,"ageing"] < .10 &
       (P2[,"FP4"] < .05 | P2[,"EGM"] < .05 | P2[,"OSK24"] < .05)
KR <- rownames(M2)[rel]
DNr <- KR[M2[KR,"EGM"] < 0]; DNr <- DNr[order(P2[DNr,"ageing"])]
kept <- list(); selDr <- character()
for (s in DNr) { g <- ALLS[[s]]
  if (!any(sapply(kept, function(k) length(intersect(g,k))/min(length(g),length(k)) > .5))) {
    kept[[length(kept)+1]] <- g; selDr <- c(selDr, s) }
  if (length(selDr) >= 22) break }

## up-direction candidates, ranked by OSK 24m P.
UPr <- KR[M2[KR,"EGM"] > 0]
UPr <- UPr[order(P2[UPr,"OSK24"] , -abs(M2[UPr,"OSK24"]))]
kept <- list(); selU <- character()
for (s in UPr) { g <- ALLS[[s]]
  if (!any(sapply(kept, function(k) length(intersect(g,k))/min(length(g),length(k)) > .45))) {
    kept[[length(kept)+1]] <- g; selU <- c(selU, s) }
  if (length(selU) >= 20) break }

## pick the programmes and their short labels.
fid <- function(pool, pat) { i <- pool[grepl(pat, sapply(pool, nmx), ignore.case=TRUE)]; i[1] }
upP <- c(fid(selU,"^Hallmark G2M"), fid(selU,"^chromosome segregation"),
         fid(selU,"^positive regulation of cell cycle"), fid(selU,"MYC_TARGETS_V1"),
         fid(selU,"^Hallmark MITOTIC_SPINDLE"), "GO:0061621")
dnP <- c(fid(selDr,"^homophilic cell adhesion"), fid(selDr,"^extracellular matrix assembly"),
         fid(selDr,"^mesenchyme development"), fid(selDr,"canonical Wnt"),
         fid(selDr,"epithelial to mesenchymal"), fid(selDr,"^positive regulation of ossification"))
upP <- upP[!is.na(upP) & upP %in% rownames(M2)]; dnP <- dnP[!is.na(dnP) & dnP %in% rownames(M2)]
picks <- c(upP, dnP)
short2 <- setNames(c("G2M checkpoint","Chromosome segregation","Positive reg. of cell cycle",
                     "MYC targets","Mitotic spindle","Glycolysis (via F6P)",
                     "Homophilic cell adhesion","ECM assembly","Mesenchyme development",
                     "Canonical Wnt signalling","Neg. reg. of EMT","Ossification")[seq_along(picks)],
                   picks)

## plot data D8 and plot function mk().
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
EX3 <- do.call(rbind, lapply(names(extra), function(k) {
  v <- setNames(extra[[k]]$DDE, extra[[k]]$gene); b <- median(v)
  do.call(rbind, lapply(picks, function(s) { i <- v[intersect(ALLS[[s]], names(v))]
    data.frame(id=s, axis=k, med=median(i)-b, p=suppressWarnings(wilcox.test(i,v)$p.value)) })) }))
D8 <- bind_rows(
  as.data.frame(M2[picks, core]) |> tibble::rownames_to_column("id") |>
    pivot_longer(-id, names_to="axis", values_to="med") |>
    left_join(as.data.frame(P2[picks, core]) |> tibble::rownames_to_column("id") |>
                pivot_longer(-id, names_to="axis", values_to="p"), by=c("id","axis")), EX3) |>
  mutate(prog = factor(short2[id], short2[picks]),
         axis = factor(lab6[axis], rev(lab6)),
         v = pmax(pmin(med, .5), -.5), lp = pmin(-log10(p), 12),
         dir = factor(ifelse(id %in% upP, "restored by treatment", "suppressed by treatment"),
                      c("restored by treatment","suppressed by treatment")))
mk <- function(bare) {
  p <- ggplot(D8, aes(prog, axis)) +
    geom_hline(yintercept=1.5, colour="grey85", linewidth=.5) +
    geom_hline(yintercept=5.5, colour="grey85", linewidth=.5) +
    geom_point(aes(size=lp, fill=v), shape=21, colour="grey25", stroke=.35) +
    facet_grid(. ~ dir, scales="free_x", space="free_x") +
    scale_fill_gradient2(low="#4A7DB5", mid="white", high="#C4443C", midpoint=0,
                         limits=c(-.5,.5), name=expression(median~log[2]*FC)) +
    scale_size_continuous(range=RNG, limits=c(0,12), name=expression(-log[10]*P), breaks=c(2,5,10)) +
    scale_y_discrete(expand=expansion(add=.7)) + scale_x_discrete(expand=expansion(add=.7)) +
    labs(x=NULL, y=NULL) + theme_bw(base_size=10) +
    theme(panel.grid=element_blank(),
          panel.border=element_rect(colour="black", fill=NA, linewidth=.5),
          strip.background=element_blank())
  if (bare) p + theme(strip.text=element_blank(), axis.text=element_blank(),
                      axis.ticks=element_blank(), legend.position="none",
                      plot.margin=margin(4,4,4,4))
  else p + theme(strip.text=element_text(face="bold", size=9),
                 axis.text.x=element_text(angle=45, hjust=1, size=9),
                 axis.text.y=element_text(face="bold", size=10),
                 legend.key.size=unit(11,"pt"))
}

## drop Ossification and Homophilic cell adhesion.
drop <- picks[short2[picks] %in% c("Ossification","Homophilic cell adhesion")]
picks <- setdiff(picks, drop); upP <- intersect(upP, picks); dnP <- intersect(dnP, picks)
D8 <- D8 |> filter(!id %in% drop) |> mutate(prog = droplevels(prog))

## add Mesenchymal cell proliferation.
add <- fid(selDr, "^regulation of mesenchymal cell proliferation")
short2[add] <- "Mesenchymal cell proliferation"
picks <- c(picks, add); dnP <- c(dnP, add)
newrow <- bind_rows(
  data.frame(id = add, axis = core, med = M2[add, core], p = P2[add, core]),
  do.call(rbind, lapply(names(extra), function(k) {
    v <- setNames(extra[[k]]$DDE, extra[[k]]$gene); i <- v[intersect(ALLS[[add]], names(v))]
    data.frame(id = add, axis = k, med = median(i) - median(v),
               p = suppressWarnings(wilcox.test(i, v)$p.value)) })))
D8 <- bind_rows(D8 |> mutate(prog = as.character(prog)),
                newrow |> mutate(prog = short2[add],
                                 axis = as.character(factor(lab6[axis], rev(lab6))),
                                 v = pmax(pmin(med,.5),-.5), lp = pmin(-log10(p),12),
                                 dir = "suppressed by treatment")) |>
  mutate(prog = factor(prog, short2[picks]),
         axis = factor(axis, rev(lab6)),
         dir = factor(dir, c("restored by treatment","suppressed by treatment")))

## non-cell-cycle up candidates.
ccpat <- "cell cycle|mitot|meiot|chromosom|spindle|cytokinesis|centrosome|nuclear division|segregation|kinetochore|cell division|DNA replication|G2M|E2F"
UPr2 <- KR[M2[KR,"EGM"] > 0]
UPr2 <- UPr2[!grepl(ccpat, sapply(UPr2, nmx), ignore.case = TRUE)]
UPr2 <- UPr2[order(P2[UPr2,"OSK24"], -abs(M2[UPr2,"OSK24"]))]
kept <- list(); selU2 <- character()
for (s in UPr2) { g <- ALLS[[s]]
  if (!any(sapply(kept, function(k) length(intersect(g,k))/min(length(g),length(k)) > .45))) {
    kept[[length(kept)+1]] <- g; selU2 <- c(selU2, s) }
  if (length(selU2) >= 22) break }

## halve the cell-cycle rows; add Energy metabolism, Deubiquitination, Pol II transcription.
keepU <- picks[short2[picks] %in% c("G2M checkpoint","Chromosome segregation","MYC targets","Glycolysis (via F6P)")]
newU  <- c(fid(selU2,"^generation of precursor metabolites"),
           fid(selU2,"^protein modification by small protein removal"),
           fid(selU2,"^transcription initiation at RNA polymerase II"))
newU  <- newU[!is.na(newU)]
short2[newU] <- c("Energy metabolism","Deubiquitination","Pol II transcription initiation")[seq_along(newU)]
upP   <- c(keepU, newU)
picks <- c(upP, dnP)

rows <- do.call(rbind, lapply(newU, function(s) bind_rows(
  data.frame(id=s, axis=core, med=M2[s,core], p=P2[s,core]),
  do.call(rbind, lapply(names(extra), function(k) {
    v <- setNames(extra[[k]]$DDE, extra[[k]]$gene); i <- v[intersect(ALLS[[s]], names(v))]
    data.frame(id=s, axis=k, med=median(i)-median(v),
               p=suppressWarnings(wilcox.test(i,v)$p.value)) })))))
D8 <- bind_rows(
  D8 |> filter(id %in% picks) |> mutate(prog=as.character(prog), axis=as.character(axis)),
  rows |> mutate(prog=short2[id], axis=as.character(factor(lab6[axis], rev(lab6))),
                 v=pmax(pmin(med,.5),-.5), lp=pmin(-log10(p),12),
                 dir="restored by treatment")) |>
  mutate(prog=factor(prog, short2[picks]), axis=factor(axis, rev(lab6)),
         dir=factor(dir, c("restored by treatment","suppressed by treatment")))

## lipid-related sets.
lip <- rownames(MED)[grepl("lipid|fatty acid|cholesterol|sterol|lipoprotein|triglyceride|phospholipid|adipo|lipogen|beta-oxidation",
                           sapply(rownames(MED), nmx), ignore.case = TRUE)]
PL <- t(sapply(lip, function(s) sapply(names(V), function(k) {
  v <- V[[k]]; suppressWarnings(wilcox.test(v[intersect(ALLS[[s]], names(v))], v)$p.value) })))
colnames(PL) <- names(V); ML <- MED[lip, , drop = FALSE]
# rank by evidence in the treatments
sc <- -log10(PL[,"FP4"]) - log10(PL[,"EGM"])
ordL <- lip[order(-sc)]
kept <- list(); selL <- character()
for (s in ordL) { g <- ALLS[[s]]
  if (!any(sapply(kept, function(k) length(intersect(g,k))/min(length(g),length(k)) > .5))) {
    kept[[length(kept)+1]] <- g; selL <- c(selL, s) }
  if (length(selL) >= 16) break }

## add Fatty acid transport.
fat <- selL[grepl("^fatty acid transport$", sapply(selL, nmx), ignore.case = TRUE)][1]
short2[fat] <- "Fatty acid transport"
upP <- c(upP, fat); picks <- c(upP, dnP)
rowF <- bind_rows(
  data.frame(id = fat, axis = core, med = ML[fat, core], p = PL[fat, core]),
  do.call(rbind, lapply(names(extra), function(k) {
    v <- setNames(extra[[k]]$DDE, extra[[k]]$gene); i <- v[intersect(ALLS[[fat]], names(v))]
    data.frame(id = fat, axis = k, med = median(i) - median(v),
               p = suppressWarnings(wilcox.test(i, v)$p.value)) })))
D8 <- bind_rows(D8 |> mutate(prog = as.character(prog), axis = as.character(axis)),
                rowF |> mutate(prog = short2[fat],
                               axis = as.character(factor(lab6[axis], rev(lab6))),
                               v = pmax(pmin(med,.5),-.5), lp = pmin(-log10(p),12),
                               dir = "restored by treatment")) |>
  mutate(prog = factor(prog, short2[picks]), axis = factor(axis, rev(lab6)),
         dir = factor(dir, c("restored by treatment","suppressed by treatment")))

## drop Mesenchyme development; save the labelled version.
drop2 <- picks[short2[picks] == "Mesenchyme development"]
picks <- setdiff(picks, drop2); dnP <- setdiff(dnP, drop2)
D8 <- D8 |> filter(!id %in% drop2) |> mutate(prog = droplevels(prog))
ggsave("RRHO_new/reproduced/shared_programmes_dot.png", mk(FALSE), width=9.4, height=4.5, dpi=300)
ggsave("RRHO_new/reproduced/shared_programmes_dot.pdf", mk(FALSE), width=9.4, height=4.5)

## ============ 8. final text-free panel ============

cat("up  :", paste(short2[upP], collapse = " | "), "\n")
cat("down:", paste(short2[dnP], collapse = " | "), "\n")
ggsave("RRHO_new/reproduced/shared_programmes_dot_bare.png", mk(TRUE), width = 7.4, height = 3.0, dpi = 300)
ggsave("RRHO_new/reproduced/shared_programmes_dot_bare.pdf", mk(TRUE), width = 7.4, height = 3.0)
print(mk(TRUE))

