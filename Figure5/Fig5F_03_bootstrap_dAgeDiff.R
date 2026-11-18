# Bootstrap interval for ΔAgeDiff (clock design of Fig5F_01_native_clock_ymctrl.R).
# Each iteration resamples, with replacement and within group, the public Y and M animals, the three saline libraries
# (training) and the three EGM libraries (projected), then reruns the whole pipeline: aging-gene selection (limma,
# p < 0.05), elastic net with alpha fixed at the value in native_clock_ymctrl_summary.csv and lambda by leave-one-sample-out CV
# (copies of a resampled sample share a fold), saline as leave-one-library-out and EGM projected. Saline draws with a
# single distinct library are redrawn (about 1 in 9), since leaving that library out would remove the 22-yr anchor.
# Output: bootstrap_dAgeDiff.csv (per iteration) and bootstrap_dAgeDiff_summary.csv (point estimate + 95% percentile).
suppressMessages({library(glmnet); library(limma); library(data.table); library(parallel)})
B <- "monkey/reproduce data/repro_native_clock_egm/data"
OUT <- "monkey/reproduce data/repro_native_clock_egm/reviewer_checks_2026-11"
ALPHA <- c(Brain = 0.1, Liver = 0.55, Kidney = 0.9)   # outputs/native_clock_ymctrl_summary.csv
ages  <- c(Y = 5, M = 11, CTRL = 22)
NB <- 1000

one <- function(cpm, meta, a, resample = TRUE) {
  pick <- function(g) { s <- rownames(meta)[meta$group == g]; if (!resample) return(s)
    repeat { r <- sample(s, length(s), replace = TRUE)          # the saline animal is the only 22-yr anchor: a draw with
      if (g != "CTRL" || length(unique(r)) >= 2) return(r) } }  # one distinct library would lose it under leave-one-out
  tr <- unlist(lapply(names(ages), pick)); te <- pick("EGM")
  age <- unname(ages[meta[tr, "group"]])
  Xtr <- cpm[, tr, drop = FALSE]; keep <- rowMeans(Xtr > 1) >= 0.5
  tt <- topTable(eBayes(lmFit(Xtr[keep, ], model.matrix(~age))), coef = "age", number = Inf, sort.by = "none")
  ag <- rownames(tt)[tt$P.Value < 0.05]
  fit_cv <- function(idx) { f <- as.integer(factor(tr[idx]))
    cv.glmnet(t(Xtr[ag, idx, drop = FALSE]), age[idx], alpha = a, foldid = f, grouped = FALSE, type.measure = "mae") }
  full <- fit_cv(seq_along(tr))
  egm <- mean(predict(full, t(cpm[ag, te, drop = FALSE]), s = "lambda.min"))
  sal_ids <- unique(tr[meta[tr, "group"] == "CTRL"])
  sal <- mean(sapply(sal_ids, function(s) {           # leave this library (and its copies) out
    idx <- which(tr != s); predict(fit_cv(idx), t(cpm[ag, s, drop = FALSE]), s = "lambda.min") }))
  egm - sal
}

res <- rbindlist(lapply(names(ALPHA), function(t) {
  cnt <- as.matrix(read.csv(file.path(B, sprintf("%s_counts.csv", t)), row.names = 1, check.names = FALSE))
  meta <- read.csv(file.path(B, sprintf("%s_meta.csv", t)), row.names = 1)
  cpm <- log2(sweep(cnt, 2, colSums(cnt), "/") * 1e6 + 1)
  set.seed(1); pt <- one(cpm, meta, ALPHA[[t]], resample = FALSE)
  RNGkind("L'Ecuyer-CMRG"); set.seed(2026)
  bs <- unlist(mclapply(seq_len(NB), function(i) tryCatch(one(cpm, meta, ALPHA[[t]]), error = function(e) NA_real_),
                        mc.cores = max(1, detectCores() - 2), mc.set.seed = TRUE))
  cat(t, "point", round(pt, 2), "boot median", round(median(bs, na.rm = TRUE), 2),
      "95%", round(quantile(bs, c(.025, .975), na.rm = TRUE), 2), "failed", sum(is.na(bs)), "\n")
  data.table(tissue = t, iter = c(0, seq_len(NB)), dAgeDiff = c(pt, bs))
}))
fwrite(res, file.path(OUT, "bootstrap_dAgeDiff.csv"))
sm <- res[iter > 0, .(boot_median = median(dAgeDiff, na.rm = TRUE), lo95 = quantile(dAgeDiff, .025, na.rm = TRUE),
                      hi95 = quantile(dAgeDiff, .975, na.rm = TRUE), frac_below0 = mean(dAgeDiff < 0, na.rm = TRUE),
                      n_ok = sum(!is.na(dAgeDiff))), by = tissue]
sm <- merge(res[iter == 0, .(tissue, point = dAgeDiff)], sm, by = "tissue")
fwrite(sm, file.path(OUT, "bootstrap_dAgeDiff_summary.csv")); print(sm)
