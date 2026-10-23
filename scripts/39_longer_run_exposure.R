# 39: FULL PANEL: longer-run within-pair exposure with calendar-year windows (dyad fixed-effects logit). Run with ALLIANCE_DEF=defense.
#      Predictors are measured in year t (or earlier) and the outcome is the dispute beginning in t+1.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(fixest) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
d <- d %>% arrange(dyad, year) %>% mutate(onset_m = as.integer(m1 == 1))
win <- function(yrs, z, lo, hi, minn = 3) vapply(seq_along(yrs), function(i) { k <- yrs >= yrs[i] - lo & yrs <= yrs[i] - hi; if (sum(k) >= minn) mean(z[k]) else NA_real_ }, numeric(1))
d <- d %>% group_by(dyad) %>% mutate(Wa = win(year, z_mx, 4, 0), Wb = win(year, z_mx, 5, 1), Pa = win(year, z_mx, 9, 5), Pb = win(year, z_mx, 10, 6), cha = Wa - Pa, chb = Wb - Pb) %>% ungroup()
sink("results/longer_run_exposure_output.txt", split = TRUE)
cat("Windows (calendar years, >= 3 observed years): Wa = years t-4..t; Wb = t-5..t-1; change_a = Wa minus mean over t-9..t-5; change_b = Wb minus mean over t-10..t-6. Outcome: dispute begins in t+1.\n")
for (cg in c(0, 1)) { n5 <- d %>% filter(m == 0, contig == cg)
  cat(sprintf("\n%s pairs\n", ifelse(cg == 0, "NON-CONTIGUOUS", "CONTIGUOUS")))
  for (v in c("z_mx", "Wa", "Wb", "cha", "chb")) for (sp in c("", "+ z_vol")) {
    x <- n5 %>% filter(!is.na(Wa), if (v %in% c("cha")) !is.na(cha) else if (v == "chb") !is.na(chb) else if (v == "Wb") !is.na(Wb) else TRUE); x$l_dur <- log(x$dur_m); x$v <- x[[v]]
    m <- tryCatch(feglm(as.formula(paste("onset_m ~ v + z_mn + z_cap + l_dur + l_same + l_diff", sp, "| dyad")), x, family = "binomial", cluster = ~dyad, notes = FALSE), error = function(e) NULL)
    if (is.null(m)) { cat(v, "failed\n"); next }
    lab <- c(z_mx = "current-year share (same sample as Wa)", Wa = "mean, years t-4..t", Wb = "mean, years t-5..t-1", cha = "change: Wa minus mean t-9..t-5", chb = "change: Wb minus mean t-10..t-6")[v]
    cat(sprintf("%-40s %-9s %s | dyads with onset %d, obs %d\n", lab, ifelse(sp == "", "baseline", "+volume"), ci(coef(m)["v"], se(m)["v"]), length(unique(x$dyad[x$onset_m == 1])), nobs(m))) } }
sink()
