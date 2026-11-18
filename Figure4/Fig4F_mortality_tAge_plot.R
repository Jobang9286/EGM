suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
d   <- fread("outputs/tage_mortality_persample.csv")
tis <- c("Brain", "Liver", "Kidney")   # column order
HX  <- 30                               # off-axis x for HGPS (chronological age not applicable)

## colour scheme — Y-Ctrl grey, M-Ctrl light blue, O-Ctrl dark blue, Saline red, EGM green
gcol <- c("Y-Ctrl"="#C9C9C9","M-Ctrl"="#7FB2E5","O-Ctrl"="#2166AC","Saline"="#C00000","EGM"="#70AD47")
gof  <- function(a) ifelse(a %in% c(1,3), "Y-Ctrl", ifelse(a %in% c(6,12,18), "M-Ctrl", "O-Ctrl"))

mk <- function(t){
  dd <- d[tissue==t]; wt <- dd[!is.na(age)]
  hc <- dd[treatment=="HGPS_ctrl"]; he <- dd[treatment=="HGPS_EGM"]
  Rr  <- cor(wt$tAge, wt$age)
  dlt <- mean(he$tAge) - mean(hc$tAge)

  ## --- top: scatter (WT regression + WT R; HGPS at x=HX) ---
  sc <- rbind(data.table(x=wt$age, y=wt$tAge, grp=gof(wt$age)),
              data.table(x=HX, y=hc$tAge, grp="Saline"),
              data.table(x=HX, y=he$tAge, grp="EGM"))
  sc$grp <- factor(sc$grp, levels=names(gcol))
  top <- ggplot(sc, aes(x,y)) +
    geom_smooth(data=data.table(x=wt$age,y=wt$tAge), aes(x,y), method="lm", se=TRUE,
                color="black", fill="grey80", linewidth=.6) +
    geom_vline(xintercept=28.5, linetype=3, color="grey65") +
    geom_point(aes(color=grp), size=2) + scale_color_manual(values=gcol, drop=FALSE, name=NULL) +
    scale_x_continuous(limits=c(0,32), breaks=c(0,6,12,18,24,30),
                       labels=c("0","6","12","18","24","HGPS")) +
    labs(x="Chronological age (month)", y="Mortality tAge (log10HR)", title=NULL) +
    annotate("text", x=0.5, y=max(sc$y), hjust=0, vjust=1,
             label=sprintf("WT R = %.3f", Rr), size=3.5) +
    theme_classic(base_size=12) + theme(axis.text=element_text(color="black"))

  ## --- bottom: boxplot (Y/M/O-Ctrl + Saline + EGM) ---
  bx <- rbind(data.table(grp=gof(wt$age), y=wt$tAge),
              data.table(grp="Saline", y=hc$tAge),
              data.table(grp="EGM",    y=he$tAge))
  bx$grp <- factor(bx$grp, levels=names(gcol))
  bot <- ggplot(bx, aes(grp,y,color=grp)) +
    geom_hline(yintercept=0, linetype=2, color="grey60") +
    geom_boxplot(fill=NA, outlier.shape=NA, width=.6, linewidth=.6) +
    geom_jitter(width=.1, size=1.5, show.legend=FALSE) +
    scale_color_manual(values=gcol, drop=FALSE, name=NULL) +
    geom_vline(xintercept=3.5, linetype=3, color="grey65") +
    labs(x=NULL, y="Mortality tAge", title=sprintf("delta = %.2f", dlt)) +
    theme_classic(base_size=12) +
    theme(axis.text=element_text(color="black"), plot.title=element_text(size=10.5, hjust=.5),
          axis.text.x=element_text(angle=35, hjust=1, size=11.5))
  ## NOTE: p-value labels are not drawn. deltas / p summary:
  ##   Brain -0.52 (p 0.063), Liver -0.47 (p 0.011), Kidney -0.87 (p 0.069)  [center = WT_1M]
  list(top=top, bot=bot)
}

P   <- lapply(tis, mk)
leg <- theme(legend.position="right", legend.title=element_blank(),
             legend.text=element_text(size=11), legend.key.size=unit(5,"mm"))
row_top <- (P[[1]]$top | P[[2]]$top | P[[3]]$top) + plot_layout(guides="collect") & leg
row_bot <- (P[[1]]$bot | P[[2]]$bot | P[[3]]$bot) + plot_layout(guides="collect") & leg
fig <- row_top / row_bot
ggsave("outputs/mouse_mortality_fig4C.png", fig, width=12.2, height=7.6, dpi=160)
cat("saved outputs/mouse_mortality_fig4C.png\n")
