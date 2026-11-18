REPO_SHARED <- "../shared/enrichment_plotting"
source(file.path(REPO_SHARED, "plot_module_enrichment.R"))
# Set TDIR / PDIR here if the tables are not under the helper's default path.

mods <- c("M1","M3","M5","M6","M9","M2","M11")
figs <- list()
for (m in mods) {
  nk <- if (is.null(REACT_SEL[[m]])) 10 else 7
  figs[[m]] <- make_fig(m, half = TRUE, n_kegg = nk, save = TRUE)   # 저장 + 객체 보관
}
for (m in mods) print(figs[[m]])    # RStudio Plots: M1 -> M3 -> M5 -> M6 -> M9 -> M2 -> M11
