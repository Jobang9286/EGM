## --- setup ----------------------------------------------------------------
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

# palette carried over from the rest of the project
BP <- "#5FA8CE"; BP_D <- "#2E7CA6"; CP <- "#ADADAD"; CP_D <- "#7A7A7A"

FS_BASE <- 11 * 1.2                       # 20% up from the previous base_size

SC <- 0.75                                  # dial this: smaller = smaller panel

LW <- 1.75

## --- read the P-BP vs C-BP signature -------------------------------------
sig <- read.csv("analysis/4FP_vs_CP_signature.csv")
sig <- sig[!is.na(sig$pvalue) & !is.na(sig$lfc), ]
cat(sprintf("genes: %d   max -log10(p) = %.0f   lfc range %.2f to %.2f\n",
            nrow(sig), max(-log10(sig$pvalue)), min(sig$lfc), max(sig$lfc)))

## --- DE classes and the padj = 0.05 line ----------------------------------
LFC <- 0.5
pcut <- max(sig$pvalue[!is.na(sig$padj) & sig$padj < 0.05])   # p at which padj crosses 0.05
sig$y  <- -log10(sig$pvalue)
sig$de <- ifelse(!is.na(sig$padj) & sig$padj < .05 & sig$lfc >  LFC, "up",
          ifelse(!is.na(sig$padj) & sig$padj < .05 & sig$lfc < -LFC, "down", "ns"))
cat(sprintf("up (P-BP) %d   down (CP) %d   ns %d\n", sum(sig$de=="up"),
            sum(sig$de=="down"), sum(sig$de=="ns")))
cat(sprintf("horizontal line at -log10(p) = %.2f  (padj = 0.05)\n", -log10(pcut)))

WV <- 3.0 * SC; HV <- 2.6 * SC

PT <- 2.2                                    # point-size dial

## --- final panel ----------------------------------------------------------
bare_theme <- function(sc, lw, box = TRUE) {
  th <- theme_classic(base_size = FS_BASE * sc) +
    theme(axis.ticks = element_line(colour = "black", linewidth = .4 * sc * lw),
          axis.ticks.length = unit(2.6 * sc, "pt"),
          axis.text = element_blank(), axis.title = element_blank(),
          plot.title = element_blank(), plot.subtitle = element_blank(),
          legend.position = "none", plot.margin = margin(6, 8, 6, 6) * sc)
  if (box) th + theme(axis.line = element_blank(),
                      panel.border = element_rect(colour = "black", fill = NA,
                                                  linewidth = .4 * sc * lw))
  else     th + theme(axis.line = element_line(colour = "black", linewidth = .4 * sc * lw))
}

v_bare <- ggplot(sig[order(sig$de != "ns"), ], aes(lfc, y)) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dotdash",
             linewidth = .35 * SC * LW, colour = "grey25") +
  geom_hline(yintercept = -log10(pcut), linetype = "dotdash",
             linewidth = .35 * SC * LW, colour = "grey25") +
  geom_point(aes(colour = de), size = PT * SC, stroke = 0) +
  scale_colour_manual(values = c(ns = "black", down = CP, up = BP)) +
  scale_x_continuous(breaks = seq(-4, 4, 2), limits = c(-4.4, 4.4)) +
  labs(x = NULL, y = NULL) + bare_theme(SC, LW, box = TRUE)

ggsave("analysis/figures/volcano_4FBP_bare.pdf", v_bare, width = WV, height = HV)
ggsave("analysis/figures/volcano_4FBP_bare.png", v_bare, width = WV, height = HV, dpi = 300)
print(v_bare)
