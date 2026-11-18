## ---- output directory -------------------------------------------------------
fig_o <- "out"

## ---- read DESeq2 table; DEG counts per cut-off ------------------------------
fig_v <- read.csv("data/egm_DESeq2_all.csv", row.names = 1)
cat("EGM vs no-EGM DESeq2:", nrow(fig_v), "genes | +LFC = EGM에서 up (Mmp3 =",
    round(fig_v["Mmp3","log2FoldChange"],2), ")\n")
cat("max -log10 p =", round(max(-log10(fig_v$pvalue), na.rm=TRUE), 1), "\n\n")
cat("컷별 DEG 수 (원본 그림: 왼쪽 1,097 / 오른쪽 1,250)\n")
for (lf in c(0.5, 0.585, 0.6, 1)) for (pc in c("pvalue","padj")) {
  s <- fig_v[!is.na(fig_v[[pc]]) & fig_v[[pc]] < 0.05, ]
  cat(sprintf("  |LFC|>%.3f & %s<0.05  ->  down %4d / up %4d\n", lf, pc,
              sum(s$log2FoldChange < -lf), sum(s$log2FoldChange > lf))) }

## ---- -log10 p (capped) and groups -------------------------------------------
cat("pvalue == 0 인 유전자:", sum(fig_v$pvalue == 0, na.rm = TRUE), "개\n")
fig_v <- fig_v[!is.na(fig_v$pvalue) & !is.na(fig_v$log2FoldChange), ]
fig_v$nlp <- -log10(fig_v$pvalue)
fig_cap <- 250
cat("cap", fig_cap, "초과(또는 Inf):", sum(fig_v$nlp > fig_cap), "개 -> 250 으로 고정\n")
fig_v$nlp_c <- pmin(fig_v$nlp, fig_cap)
fig_v$grp <- with(fig_v, ifelse(pvalue < .05 & log2FoldChange >  .5, "EGM",
                         ifelse(pvalue < .05 & log2FoldChange < -.5, "Control", "n.s.")))
print(table(fig_v$grp))

## ---- final volcano ----------------------------------------------------------
fig_volc <- function() {
  col <- c(`n.s.` = "#E3E3E3", Control = "#8C8C8C", EGM = "#5B7CFA")
  par(mar = c(4.2, 4.4, 3.6, 1.2))
  o <- order(match(fig_v$grp, c("n.s.","Control","EGM")))
  plot(fig_v$log2FoldChange[o], fig_v$nlp_c[o], pch = 16, cex = .55,
       col = col[fig_v$grp[o]], xlab = expression(Log[2]*FoldChange),
       ylab = "-log(Pvalue)", las = 1, bty = "o",          # 사방 테두리
       xlim = c(-9, 13), ylim = c(0, 255))
  box(lwd = 1.2)
  abline(v = c(-.5, .5), lty = 2, lwd = 1.6)
  abline(h = -log10(.05), lty = 2, lwd = 1.6)
  legend("top", legend = c("Control","EGM"), fill = col[c("Control","EGM")], border = NA,
         horiz = TRUE, bty = "n", inset = c(0, -0.14), xpd = TRUE, cex = 1.05)
}
pdf(file.path(fig_o, "volcano_EGM_vs_noEGM.pdf"), width = 5.2, height = 5.2); fig_volc(); dev.off()
