suppressMessages({library(tidyverse); library(patchwork); library(ggrepel)})
D <- "chromatin/tracks"
out <- "chromatin/figures"

PAL  <- c(`P-BP` = "#5FA8CE", CP = "#ADADAD")     # deeper blue / deeper grey
OUTL <- c(`P-BP` = "#2E7CA6", CP = "#7A7A7A")

gm <- tribble(
  ~locus, ~gene, ~start, ~end,
  "Cdkn2a_locus","Cdkn2a", 89192707, 89212890,
  "Cdkn2a_locus","Cdkn2b", 89224535, 89229276,
  "HoxA_cluster","Hoxa1",  52132569, 52135297,
  "HoxA_cluster","Hoxa5",  52178733, 52181567,
  "HoxA_cluster","Hoxa9",  52200079, 52208069,
  "HoxA_cluster","Hoxa10", 52208176, 52217834,
  "HoxA_cluster","Hoxa13", 52234673, 52237788,
  "Ccnb1","Cenph",  100896181,100912407,
  "Ccnb1","Ccnb1",  100915157,100923078,
  "Itga5","Zfp385a",103222321,103248520,
  "Itga5","Itga5",  103252712,103275190)

mkpan <- function(d, loc, ttl, sub = NULL){
  dd <- filter(d, locus==loc); rng <- range(c(dd$start, dd$end))/1e6
  gg <- filter(gm, locus==loc)
  trk <- ggplot(dd, aes(mid/1e6, v, fill=arm, colour=arm)) +
    geom_area(linewidth=.2) + geom_blank(aes(y=ymax)) +
    facet_grid(arm ~ ., switch="y") +
    scale_fill_manual(values=PAL, guide="none") +
    scale_colour_manual(values=OUTL, guide="none") +
    scale_x_continuous(limits=rng, expand=expansion(mult=.012)) +
    scale_y_continuous(n.breaks=2, expand=expansion(c(0,.05))) +
    labs(x=NULL, y=NULL, title=ttl, subtitle=sub) +
    theme_classic(base_size=12) +
    theme(strip.placement="outside", strip.background=element_blank(),
          strip.text.y.left=element_text(angle=0, size=10, face="bold", hjust=1,
                                         margin=margin(r=5)),
          axis.line.y=element_line(colour="black", linewidth=.4), axis.line.x=element_blank(),
          axis.ticks.y=element_line(colour="black", linewidth=.4),
          axis.ticks.length=unit(3,"pt"),
          axis.text.y=element_text(size=8.5, colour="black"),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          plot.title=element_text(size=13, face="bold", hjust=.5, margin=margin(b=2)),
          plot.subtitle=element_text(size=8.5, colour="grey35", hjust=.5, margin=margin(b=3)),
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

prep <- function(f, keep=NULL){
  x <- read_tsv(file.path(D,f), show_col_types=FALSE, progress=FALSE)
  if (!is.null(keep)) x <- filter(x, mark==keep)
  x |> pivot_longer(c(sig_4FP,sig_CP), names_to="arm", values_to="v") |>
    mutate(arm = factor(if_else(arm=="sig_4FP","P-BP","CP"), c("P-BP","CP")),
           mid = (start+end)/2) |>
    group_by(locus) |> mutate(ymax = max(v)) |> ungroup()
}

# H3K27me3
d1 <- prep("locus_signal.tsv", "H3K27me3")
p1 <- mkpan(d1,"Cdkn2a_locus","Cdkn2a") / mkpan(d1,"HoxA_cluster","HoxA cluster")
ggsave(file.path(out,"locus_tracks_H3K27me3_Cdkn2a_HoxA.pdf"), p1, width=6.0, height=4.2)
ggsave(file.path(out,"locus_tracks_H3K27me3_Cdkn2a_HoxA.png"), p1, width=6.0, height=4.2, dpi=300)
