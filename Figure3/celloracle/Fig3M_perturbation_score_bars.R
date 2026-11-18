# ============================================================
# lipid TFs (reversing / no effect) + non-lipid controls
#
# Block 1 (red)  : 17 lipid-metabolism TFs whose KO reverses the rejuvenation flow
# Block 2 (grey) : lipid-metabolism TFs without effect (|PS| < 0.01)  -> not every lipid TF acts
# Block 3 (grey) : non-lipid control TFs, mean PS closest to 0 among 23+ screened
#                  (data/PS_negctrl_all.csv)
# y = perturbation score = -mean PS over the Q1-Q3 grid
# Output: figs/FigS_KO_PS_lipid_ns_controls.{pdf,png}
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
ROOT <- "EGM_oracle"; PS_THR <- 0.01

LIPID    <- c("Srebf2","Nr4a1","Klf4","Tfap2b","Klf5","Atf2","Tfe3","Foxo1","Klf9","Nr2f2",
              "Srebf1","Atf4","Nfe2l1","Usf2","Nfe2l2","Atf6","E2f1")
# Block 2 now = near-zero TFs from the 84-TF screen whose primary function is NOT lipid metabolism
# (they entered that screen only through broad GO terms such as "fat cell differentiation")
LIPID_NS <- c("Egr2","Zbtb7a","Mitf","Pbx1","Mecom","Gata2")
CTRL     <- c("Mbd2","Fli1")          # dedicated non-lipid control screen (data/PS_negctrl_all.csv)

a <- read.csv(file.path(ROOT, "data/PS_TF_KO_all.csv"))
b <- read.csv(file.path(ROOT, "data/PS_negctrl_all.csv"))
d <- bind_rows(
  a %>% filter(TF %in% LIPID)    %>% mutate(group = "lipid")    %>% arrange(mean_ps_early),
  bind_rows(a %>% filter(TF %in% LIPID_NS), b %>% filter(TF %in% CTRL)) %>%
    mutate(group = "control") %>% arrange(mean_ps_early)) %>%
  select(TF, mean_ps_early, group) %>%
  mutate(score = -mean_ps_early, TF = factor(TF, levels = TF))
stopifnot(nrow(d) == length(c(LIPID, LIPID_NS, CTRL)))
print(d %>% mutate(score = round(score, 4)), row.names = FALSE)

n1 <- sum(d$group == "lipid"); n3 <- sum(d$group == "control")
ytop <- max(d$score) * 1.08
p <- ggplot(d, aes(TF, score, fill = group)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -PS_THR, ymax = PS_THR, fill = "grey85", alpha = 0.8) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = n1 + 0.5, linetype = "22", linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values = c(lipid = "#B2182B", lipid_ns = "grey65", control = "grey65"), guide = "none") +
  scale_y_continuous(breaks = c(0, 0.10), labels = c("0", "0.1"),
                     expand = expansion(c(0.04, 0.10))) +
  labs(title = "In silico perturbation", x = NULL, y = "Perturbation score") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6), axis.line = element_blank())

ggsave(file.path(ROOT, "figs/FigS_KO_PS_lipid_ns_controls.pdf"), p, width = 6.8, height = 3.6, device = "pdf")
ggsave(file.path(ROOT, "figs/FigS_KO_PS_lipid_ns_controls.png"), p, width = 6.8, height = 3.6, dpi = 300, bg = "white")
cat("saved: figs/FigS_KO_PS_lipid_ns_controls.{pdf,png}\n")
