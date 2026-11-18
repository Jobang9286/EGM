# ==============================================================================
# EGM-Young module along OSKM reprogramming pseudotime — figure code
#
# Source data (both already contain everything needed; nothing is recomputed):
#   SOT/sce_seurat.RData
#       sce.seurat     Seurat, 20,273 genes x 912 cells, assay "originalexp"
#                      meta.data already carries the AddModuleScore columns
#                      young_ear_signature1 (65 genes) and young_skin_signature1 (52)
#       pseudotimeED   912 x 1 slingshot pseudotime, Lineage1, MEF 1.0 -> ESC 31.6
#   Upstream analysis: code_SOT1.R
#
# Gene lists: <root>/ref/genelist_ear_uu.list.csv  (65, fibroblast)
#             <root>/ref/uu.list.csv               (52, skin, no header)
# ==============================================================================
suppressMessages({library(Seurat); library(ggplot2); library(ggh4x); library(Matrix)})

root   <- "reproduce"
outdir <- file.path(root, "out")

## ---- 1. load -----------------------------------------------------------------
e <- new.env(); load("SOT/sce_seurat.RData", envir = e)   # ~10 s, ~3.5 GB
so   <- e$sce.seurat
expr <- GetAssayData(so, assay = "originalexp", layer = "data")
pt   <- e$pseudotimeED[colnames(so), 1]

## ---- 2. timepoints ------------------------------------------------------------
# osk_d4 has 0 cells and is dropped; iPSC (31.2) and ESC (31.6) are nearly identical
lab <- c(mef="MEF", osk_d0="Day0", osk_d1="Day1", osk_d2="Day2", osk_d3="Day3",
         osk_d5="Day5", osk_d7="Day7", osk_d8="Day8", ips="iPSC", esc="ESC")
tp <- factor(unname(lab[as.character(so$Day)]), levels = unname(lab))
names(tp) <- colnames(so); so$Timepoint <- droplevels(tp)

## ---- 3. even-spaced pseudotime axis -------------------------------------------
# Raw pseudotime is very unevenly spaced (MEF 1.0, Day1 15.0, Day2 19.9, Day7 24.8,
# iPSC 31.2). Piecewise-linear rescale so each timepoint sits on an integer, keeping
# every cell's position *within* its interval. Curve shape is preserved locally.
anchor_pt <- tapply(pt, so$Timepoint, median)
anchor_pt <- c(anchor_pt[c("MEF","Day0","Day1","Day2","Day3","Day5","Day7")],
               `iPSC/ESC` = unname(mean(anchor_pt[c("iPSC","ESC")])))
ax <- c(min(pt), unname(anchor_pt), max(pt))
ay <- c(1 - (anchor_pt[1] - min(pt)) / (anchor_pt[2] - anchor_pt[1]),
        seq_along(anchor_pt),
        length(anchor_pt) + (max(pt) - tail(anchor_pt,1)) /
          (tail(anchor_pt,1) - tail(anchor_pt,2)[1]))
remap <- function(p) approx(ax, ay, xout = p, rule = 2)$y

## ---- 4. gene groups -----------------------------------------------------------
# Selected for an inverted-U along pseudotime (rise mid, fall before pluripotency).
# Dropped after inspection: Pfkfb3 (quadratic p = 0.18, flat), Angptl4 (quadratic
# +7.3, rises INTO iPSC), Bst1 / Mt2 / Upp1 (same problem), Mmp3 / Mmp13 / Hp / Lcn2
# (absent or all-zero in this object). Sod3 is kept but is detected in only 2.2% of
# cells — its post-Day3 upturn is not trustworthy.
groups <- list(`1. Early reprogramming` = c("Cd24a","Nt5e"),
               `2. Metabolic`           = c("Gch1","Lipg"),
               `3. Redox (decrease)`    = c("Sod3","Grem1"))
pal <- c(Cd24a="#E74A48", Nt5e="#C0392B", Gch1="#4E7BAB", Lipg="#1F4E79",
         Sod3="#7FBF70", Grem1="#4E8C3F")

d <- do.call(rbind, lapply(names(groups), function(k)
  do.call(rbind, lapply(groups[[k]], function(g)
    data.frame(group = k, gene = g, xe = remap(pt), z = as.numeric(scale(expr[g, ])))))))
d$gene  <- factor(d$gene,  levels = unlist(groups))
d$group <- factor(d$group, levels = names(groups))

## ---- 5. inverted-U test --------------------------------------------------------
# z ~ poly(pseudotime, 2); a significant NEGATIVE quadratic term = rises then falls.
stats <- t(sapply(unlist(groups), function(g) {
  s <- summary(lm(as.numeric(scale(expr[g, ])) ~ poly(pt, 2)))$coefficients[3, c(1,4)]
  c(pct_cells = round(100 * mean(expr[g, ] > 0), 1), quad = round(s[1], 1), p = signif(s[2], 2)) }))
print(stats)
#           pct_cells   quad         p
#   Cd24a        91.1  -22.3  3.1e-160
#   Nt5e         67.5  -15.1   1.8e-59
#   Gch1         42.0   -8.0   4.7e-16
#   Lipg         25.3   -4.2   2.3e-05
#   Sod3          2.2   +4.0   5.7e-05   <- U-shaped, and only 20 cells
#   Grem1        26.9   +1.9   3.7e-02   <- declines from MEF

## ---- 6. figure ----------------------------------------------------------------
p <- ggplot(d, aes(xe, z, colour = gene, fill = gene)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = .3) +
  geom_smooth(method = "loess", span = .75, se = TRUE, alpha = .3, linewidth = .95) +
  facet_wrap(~ group, nrow = 1, scales = "free_y") +
  # only the Metabolic panel gets 3 y ticks; its range is ~1/4 of the others
  facetted_pos_scales(y = list(
    `1. Early reprogramming` = scale_y_continuous(),
    `2. Metabolic`           = scale_y_continuous(breaks = scales::breaks_pretty(n = 3)),
    `3. Redox (decrease)`    = scale_y_continuous())) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Reprogramming pseudotime", y = "Scaled expression",
       title = "Transient expression of EGM-Young module genes during reprogramming") +
  theme_bw(base_size = 10) +
  theme(panel.grid  = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.65),
        plot.title  = element_text(hjust = .5, size = 10),
        strip.background = element_rect(fill = "grey93"),
        strip.text  = element_text(colour = NA),      # strip boxes kept, labels hidden
        legend.position = "right", legend.text = element_text(size = 8.5),
        legend.key.height = unit(12, "pt"))

ggsave(file.path(outdir, "panelN_grouped_3panel.pdf"), p, width = 9, height = 3.2)
p
