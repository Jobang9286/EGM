# =============================================================================
# ATAC tracks at Ccnb1 and Itga5, P-BP vs CP (chromatin/figures/locus_tracks_ATAC_Ccnb1_Itga5.*)
# Input: chromatin/tracks/atac_locus_signal2.tsv (100-bp bins of the CPM bigWigs),
#        chromatin/genes_de_long.tsv (RNA log2FC / padj for the subtitles)
# =============================================================================
suppressMessages({library(tidyverse); library(patchwork); library(ggrepel)})
D   <- "chromatin/tracks"
out <- "chromatin/figures"

PAL  <- c(`P-BP` = "#5FA8CE", CP = "#ADADAD")
OUTL <- c(`P-BP` = "#2E7CA6", CP = "#7A7A7A")

gm <- tribble(
  ~locus, ~gene, ~start, ~end,
  "Ccnb1","Cenph",  100896181,100912407,
  "Ccnb1","Ccnb1",  100915157,100923078,
  "Itga5","Zfp385a",103222321,103248520,
  "Itga5","Itga5",  103252712,103275190)

# ytext: size of the y tick labels; ymarg_r: their right margin in pt (2.4 = ggplot default),
# raised for a smaller label so the panel stays where it was; hide: draw title, subtitle and strip labels
# in transparent ink (they still take their space, so the layout does not move); a hidden strip label is
# padded by STRIP_PAD pt on the left so the track sits where the Supp Fig. 1J crop box expects it
STRIP_PAD <- 5.0187
mkpan <- function(d, loc, ttl, sub = NULL, ytext = 8.5, ymarg_r = 2.4, hide = FALSE,
                  pal = PAL, outl = OUTL){
  dd <- filter(d, locus==loc); rng <- range(c(dd$start, dd$end))/1e6
  gg <- filter(gm, locus==loc)
  ink <- if (hide) NA else "black"
  trk <- ggplot(dd, aes(mid/1e6, v, fill=arm, colour=arm)) +
    geom_area(linewidth=.2) + geom_blank(aes(y=ymax)) +
    facet_grid(arm ~ ., switch="y") +
    scale_fill_manual(values=pal, guide="none") +
    scale_colour_manual(values=outl, guide="none") +
    scale_x_continuous(limits=rng, expand=expansion(mult=.012)) +
    scale_y_continuous(n.breaks=2, expand=expansion(c(0,.05))) +
    labs(x=NULL, y=NULL, title=ttl, subtitle=sub) +
    theme_classic(base_size=12) +
    theme(strip.placement="outside", strip.background=element_blank(),
          strip.text.y.left=element_text(angle=0, size=10, face="bold", hjust=1,
                                         margin=margin(r=5, l=if (hide) STRIP_PAD else 0), colour=ink),
          axis.line.y=element_line(colour="black", linewidth=.4), axis.line.x=element_blank(),
          axis.ticks.y=element_line(colour="black", linewidth=.4),
          axis.ticks.length=unit(3,"pt"),
          axis.text.y=element_text(size=ytext, colour="black", margin=margin(r=ymarg_r)),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          plot.title=element_text(size=13, face="bold", hjust=.5, margin=margin(b=2), colour=ink),
          plot.subtitle=element_text(size=8.5, colour=if (hide) NA else "grey35", hjust=.5, margin=margin(b=3)),
          panel.spacing.y=unit(1.6,"mm"), plot.margin=margin(3,18,1,4))
  ann <- ggplot(gg) +
    geom_segment(aes(x=start/1e6, xend=end/1e6, y=0, yend=0), linewidth=3.2, colour="grey25") +
    geom_text_repel(aes(x=(start+end)/2e6, y=0, label=gene), size=3.1, colour="grey10",
                    nudge_y=-1, direction="x", segment.size=.28, segment.colour="grey60",
                    min.segment.length=0, box.padding=.13, max.overlaps=Inf, seed=1) +
    scale_x_continuous(limits=rng, expand=expansion(mult=.012)) +
    scale_y_continuous(limits=c(-2,.5), breaks=NULL, expand=c(0,0)) +
    labs(x=paste0(dd$chrom[1]," (Mb)"), y=NULL) +
    theme_classic(base_size=12) +
    theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          axis.line.x=element_line(colour="black", linewidth=.4),
          axis.ticks.x=element_line(colour="black", linewidth=.4),
          axis.ticks.length=unit(3,"pt"),
          axis.text.x=element_text(size=8.5, colour="black"),
          axis.title.x=element_text(size=9.5, margin=margin(t=2)),
          plot.margin=margin(1,18,3,4))
  wrap_elements(trk / ann + plot_layout(heights=c(2.6,1.6)))
}

prep <- function(f){
  read_tsv(file.path(D,f), show_col_types=FALSE, progress=FALSE) |>
    pivot_longer(c(sig_4FP,sig_CP), names_to="arm", values_to="v") |>
    mutate(arm = factor(if_else(arm=="sig_4FP","P-BP","CP"), c("P-BP","CP")),
           mid = (start+end)/2) |>
    group_by(locus) |> mutate(ymax = max(v)) |> ungroup()
}
de <- read_tsv("chromatin/genes_de_long.tsv",
               show_col_types=FALSE, progress=FALSE) |>
  filter(mark=="ATAC") |> distinct(gene, RNA=logFC, padj=adj.P.Val)
d2 <- prep("atac_locus_signal2.tsv")
sb <- function(g){ r <- filter(de, gene==g); sprintf("RNA  log2FC %+.2f   padj %.1e", r$RNA, r$padj) }
FIG_W <- 6.0; FIG_H <- 4.4

# ---- 1. both loci, with strip labels and titles ------------------------------------------------------------
p_orig <- mkpan(d2,"Ccnb1","Ccnb1", sb("Ccnb1")) / mkpan(d2,"Itga5","Itga5", sb("Itga5"))
if (exists("WRITE_ORIG") && WRITE_ORIG) {   # off by default
  ggsave(file.path(out,"locus_tracks_ATAC_Ccnb1_Itga5.pdf"), p_orig, width=FIG_W, height=FIG_H)
  ggsave(file.path(out,"locus_tracks_ATAC_Ccnb1_Itga5.png"), p_orig, width=FIG_W, height=FIG_H, dpi=300)
}

# ---- 2. Supp Fig. 1J: the Itga5 half, cropped, with smaller y tick labels so the 0 of P-BP and the 3 of
#         C-BP no longer touch.
# Crop box, in inches of the 6.0 x 4.4 in figure.
CROP <- c(x0 = 0.7635, top = 2.6844, w = 4.9005, h = 1.4549)
# Fill / outline colours for 1J (#70A6CA in place of #5FA8CE). Grey is unchanged.
PAL_PPT  <- c(`P-BP` = "#70A6CA", CP = "#ADADAD")
OUTL_PPT <- c(`P-BP` = "#447AA2", CP = "#7A7A7A")
YTEXT <- 7                                    # y tick label size (mkpan default 8.5)

draw_crop <- function(p){
  vp <- grid::viewport(x = unit(-CROP[["x0"]], "in"), y = unit(-(FIG_H - CROP[["top"]] - CROP[["h"]]), "in"),
                       width = unit(FIG_W, "in"), height = unit(FIG_H, "in"), just = c("left", "bottom"))
  grid::grid.newpage(); grid::pushViewport(vp); grid::grid.draw(patchwork::patchworkGrob(p)); grid::popViewport()
}
# keep the panel where it was: widen the label margin by what the smaller digit saves
digit_w <- function(sz) grid::convertWidth(grid::grobWidth(grid::textGrob("0", gp = grid::gpar(fontsize = sz))),
                                           "pt", valueOnly = TRUE)
make_1J <- function(ytext, pal, outl)
  mkpan(d2,"Ccnb1","Ccnb1", sb("Ccnb1")) /
  mkpan(d2,"Itga5","Itga5", sb("Itga5"), ytext = ytext, ymarg_r = 2.4 + digit_w(8.5) - digit_w(ytext),
        hide = TRUE, pal = pal, outl = outl)

f1J <- file.path(out, "SuppFig1J_ATAC_Itga5")
PX_W <- round(3.1597 * 1200)                  # 3.16 x 0.94 in at 1200 dpi
ragg::agg_png(paste0(f1J, "_1200dpi.png"), width = PX_W, height = round(PX_W * CROP[["h"]] / CROP[["w"]]),
              units = "px", res = PX_W / CROP[["w"]], background = "white")
draw_crop(make_1J(YTEXT, PAL_PPT, OUTL_PPT)); invisible(dev.off())
# ragg stamps the drawing scale (774 px per figure inch) as dpi; restamp 1200 so the file opens 3.16 in wide
png::writePNG(png::readPNG(paste0(f1J, "_1200dpi.png")), paste0(f1J, "_1200dpi.png"), dpi = 1200)
cairo_pdf(paste0(f1J, ".pdf"), width = CROP[["w"]], height = CROP[["h"]])
draw_crop(make_1J(YTEXT, PAL_PPT, OUTL_PPT)); invisible(dev.off())
cat(sprintf("Supp Fig 1J: %s_1200dpi.png (%d px wide = 3.16 in at 1200 dpi) and .pdf; y labels %.1f pt\n",
            basename(f1J), PX_W, YTEXT))
