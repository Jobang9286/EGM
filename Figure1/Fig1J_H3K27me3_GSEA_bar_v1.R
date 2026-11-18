# H3K27me3 GSEA bar: GO terms regaining H3K27me3 in P-BP.
#
# Ranking metric is promoter log2(P-BP/CP) H3K27me3 signal, TSS +/- 2 kb, all
# 21,757 protein-coding genes. GSEA is rank-based, so the global offset in this
# metric (median +0.29 between the arms) does not
# affect it; only the ordering matters.
#
# Terms are hand-picked from the significant set into two categories. The
# selection is in `want` below and is the only editorial step; everything else
# is computed.
suppressMessages({library(tidyverse); library(fgsea)})
# project root: the directory holding data/ and figures/. Works when the script
# is run from the project root, from scripts/, or via Rscript with a full path.
find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  cands <- c(if (length(f)) dirname(dirname(normalizePath(f[1]))), getwd(), dirname(getwd()))
  hit <- cands[dir.exists(file.path(cands, "data"))]
  if (!length(hit)) stop("cannot locate the data/ directory; run from the project root")
  hit[1]
}
here <- find_root()
DATA <- file.path(here, "data"); OUT <- file.path(here, "figures")
dir.create(OUT, showWarnings = FALSE)

rds <- file.path(DATA, "fgsea_h3k27me3.rds")
if (file.exists(rds)) {
  fg <- readRDS(rds)                       # cached run, identical to the block below
} else {
  gmt <- gmtPathways(file.path(DATA, "go_bp_mouse.gmt"))
  gl  <- read_lines(file.path(DATA, "go_bp_mouse.gmt"))
  nm  <- tibble(pathway = str_extract(gl, "^[^\t]+"), term = str_split_i(gl, "\t", 2))
  s <- read_tsv(file.path(DATA, "h3k27me3_promoter_gene.tsv"), show_col_types = FALSE)
  r <- setNames(s$log2FC, s$gene); r <- r[is.finite(r)]
  r <- sort(r - median(r), decreasing = TRUE)
  set.seed(1)
  fg <- as_tibble(fgsea(gmt, r, minSize = 15, maxSize = 400, nPermSimple = 10000)) |>
    left_join(nm, by = "pathway")
  saveRDS(fg, rds)
}

want <- tribble(
  ~cat, ~term,
  "Development / patterning", "embryonic organ morphogenesis",
  "Development / patterning", "cell fate commitment",
  "Development / patterning", "mesenchymal cell differentiation",
  "Development / patterning", "stem cell differentiation",
  "SASP / senescence",        "regulation of response to oxidative stress",
  "SASP / senescence",        "regulation of cellular senescence",
  "SASP / senescence",        "interleukin-6-mediated signaling pathway",
  "SASP / senescence",        "response to interleukin-1")

d <- want |> left_join(fg, by = "term")
stopifnot(!any(is.na(d$NES)), all(d$pval < 0.05))
d <- d |>
  mutate(cat  = factor(cat, c("SASP / senescence", "Development / patterning")),
         term = str_to_sentence(term)) |>
  arrange(cat, NES) |> mutate(term = fct_inorder(term))

pal <- c("Development / patterning" = "#7B5EA7", "SASP / senescence" = "#3F8F84")
p <- ggplot(d, aes(NES, term, fill = cat)) +
  geom_col(width = .68) +
  scale_fill_manual(values = pal, name = NULL,
                    breaks = c("Development / patterning", "SASP / senescence")) +
  scale_x_continuous(limits = c(0, 2.5), breaks = 0:2, expand = c(0, 0)) +
  labs(y = NULL, x = "NES", title = "H3K27me3",
       subtitle = "GSEA on promoter H3K27me3 log2(P-BP/CP), GO BP, p < 0.05") +
  theme_classic(base_size = 11) +
  theme(axis.line  = element_line(colour = "black", linewidth = .45),
        axis.ticks = element_line(colour = "black", linewidth = .45),
        axis.ticks.length = unit(3.5, "pt"),
        axis.text  = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 10, margin = margin(t = 6)),
        plot.title = element_text(size = 13, face = "bold", hjust = .5),
        plot.subtitle = element_text(size = 8, colour = "grey40", hjust = .5),
        legend.position = "top", legend.text = element_text(size = 8.5),
        legend.key.size = unit(3.5, "mm"), plot.margin = margin(6, 16, 6, 6))

ggsave(file.path(OUT, "H3K27me3_4FPspecific_NES.pdf"), p, width = 6.4, height = 3.0)
ggsave(file.path(OUT, "H3K27me3_4FPspecific_NES.png"), p, width = 6.4, height = 3.0, dpi = 300)
write_tsv(d |> select(cat, term, NES, pval, padj, size),
          file.path(OUT, "H3K27me3_4FPspecific_NES.tsv"))
message("wrote H3K27me3_4FPspecific_NES.{pdf,png,tsv}")
