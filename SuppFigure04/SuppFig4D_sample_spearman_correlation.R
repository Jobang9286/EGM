## ---- counts, sample sheet, expression filter --------------------------------
suppressPackageStartupMessages(library(DESeq2))
C9 <- read.csv("egm_blood/counts_symbol.csv", row.names = 1, check.names = FALSE)
mp <- read.delim("egm_blood/mapfile.tsv")
cat("samples:", paste(mp$sample, collapse=", "), "\n\n")
C9 <- as.matrix(C9[, mp$sample])
keep <- rowSums(C9 >= 10) >= 3

## ---- VST, top-2,000 variable genes, Spearman --------------------------------
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
ord <- c("GFP-1","GFP-2","GFP-3","Eno1-1","Gclc-1","Map1s-1","EGM-1","EGM-2","EGM-3")
grpc <- c(GFP="#9E9E9E", single="#C4A46A", EGM="#2E7CA6")
lab  <- setNames(ifelse(grepl("^GFP", ord), "GFP",
                 ifelse(grepl("^EGM", ord), "EGM", "single")), ord)

dd <- DESeqDataSetFromMatrix(round(C9[keep, ord]), data.frame(g = lab[ord]), ~ 1)
vs <- assay(vst(dd, blind = TRUE))
tv <- vs[order(-apply(vs, 1, var))[1:2000], ]
CR <- cor(tv, method = "spearman")
cat("mean Spearman vs the three GFP samples:\n")
for (s in ord[4:9]) cat(sprintf("  %-8s %.3f\n", s, mean(CR[s, c("GFP-1","GFP-2","GFP-3")])))

## ---- long-format table ------------------------------------------------------
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
hm <- as.data.frame(CR) |> tibble::rownames_to_column("a") |>
  pivot_longer(-a, names_to = "b", values_to = "r") |>
  mutate(a = factor(a, ord), b = factor(b, rev(ord)))

## ---- final heatmap ----------------------------------------------------------
g <- ggplot(hm, aes(a, b, fill = r)) +
  geom_tile(colour = "white", linewidth = .5) +
  annotate("segment", x = c(3.5, 6.5), xend = c(3.5, 6.5), y = .5, yend = 9.5,
           colour = "black", linewidth = .5) +
  annotate("segment", x = .5, xend = 9.5, y = c(3.5, 6.5), yend = c(3.5, 6.5),
           colour = "black", linewidth = .5) +
  scale_fill_gradientn(colours = c("#F7F7F7","#BBD7E8","#5FA8CE","#2E7CA6","#17394F"),
                       limits = c(0.6, 1), breaks = c(.7,.8,.9,1), name = "Spearman") +
  scale_x_discrete(position = "top", expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(x = NULL, y = NULL) + coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_text(angle = 45, hjust = 0, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.key.width = unit(8,"pt"), legend.key.height = unit(26,"pt"))
ggsave("RRHO_new/reproduced/egm_sample_correlation.png", g, width = 5.4, height = 4.4, dpi = 300)
ggsave("RRHO_new/reproduced/egm_sample_correlation.pdf", g, width = 5.4, height = 4.4)
print(g)
