#!/usr/bin/env Rscript
SCRIPT_DIR <- "."   # folder holding Fig5M-N_03_plot_organ_enrichment.R
# 원숭이 장기별 GO/KEGG 그림 전체 재생성
#   Rscript Fig5M-N_04_run_all_organs.R   (또는 RStudio 에서 source)
# 라벨 칸 폭(lab_go / lab_kegg)은 5장 모두 같은 값이어야 패널 위치가 맞는다.
source(file.path(SCRIPT_DIR, "Fig5M-N_03_plot_organ_enrichment.R"))

KEYS <- c("Kidney", "Liver", "Brain", "Spleen_up", "Spleen_down")

for (k in KEYS) {
  organ_fig(k, save = TRUE)                      # 라벨 버전 + 글자 없는 버전 + 사용 term CSV
  message(sprintf("[%s] 저장 완료", k))
}

# 라벨 칸 실측치 점검: lab_go(1.9) / lab_kegg(2.6) 를 넘으면 그 그림만 패널이 밀린다
chk <- t(sapply(KEYS, function(o) c(
  GO   = label_width_in(plot_go(organ_go(o), panel_w = 1.0, txt = 8.5, wrap = 34)),
  KEGG = label_width_in(plot_kegg(organ_kegg(o), panel_w = 0.85, txt = 8.5,
                                  n_size = 3, group = "full", wrap = 48)))))
print(round(chk, 2))
if (any(chk[, "GO"] > 1.9) || any(chk[, "KEGG"] > 2.6))
  warning("라벨 칸이 기준을 넘었다. lab_go / lab_kegg 를 키우고 5장 모두 다시 저장할 것")
