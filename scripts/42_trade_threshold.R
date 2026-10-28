# 42: FULL PANEL, dispute onset: sensitivity of the geographical interaction to the minimum-trade rule for the denominator of the trade shares.
#     Rules: (a) total trade >= US$50 million nominal (the rule used in earlier versions), (b) no minimum (positive total trade), (c) total trade >= 0.01% (the article's rule, script 02) and (d) >= 0.05% of recorded world trade that year.
#     Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
nat <- read_csv("raw/trade/COW_Trade_4.0/National_COW_4.0.csv", show_col_types = FALSE) %>%
  mutate(tt = imports + exports, tt = ifelse(tt <= 0, NA, tt)) %>% select(ccode, year, tt) %>% filter(!is.na(tt))
nat <- nat %>% group_by(year) %>% mutate(world = sum(tt), rel = tt / world) %>% ungroup()
rules <- list("US$50m nominal (earlier rule)" = function(x) x$tt >= 50, "no minimum" = function(x) x$tt > 0,
              "0.01% of recorded world trade (article)" = function(x) x$rel >= 1e-4, "0.05% of recorded world trade" = function(x) x$rel >= 5e-4)
era <- function(y) cut(y, c(1869, 1913, 1945, 1989, 2011), labels = c("1870-1913", "1914-45", "1946-89", "1990-2010"))
sink("results/trade_threshold_output.txt", split = TRUE)
cat("State-years with recorded total trade, and share meeting each rule, by era\n")
for (nm in names(rules)) { k <- nat %>% mutate(ok = rules[[nm]](cur_data()), era = era(year)) %>% group_by(era) %>% summarise(n = n(), kept = mean(ok)); cat(sprintf("%-32s %s\n", nm, paste(sprintf("%s: %d (%.0f%%)", k$era, k$n, 100 * k$kept), collapse = " | "))) }
tie0 <- read_csv("data/tie_trade.csv", show_col_types = FALSE)
dy <- read_csv("raw/trade/COW_Trade_4.0/Dyadic_COW_4.0.csv", show_col_types = FALSE) %>% mutate(across(c(flow1, flow2, smoothtotrade), ~ na_if(., -9))) %>%
  mutate(ccode_low = pmin(ccode1, ccode2), ccode_high = pmax(ccode1, ccode2)) %>% group_by(year, ccode_low, ccode_high) %>% summarise(dyad = sum(smoothtotrade, na.rm = TRUE), .groups = "drop")
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
for (nm in names(rules)) {
  den <- nat %>% filter(rules[[nm]](.)) %>% select(ccode, year, tt)
  alt <- dy %>% left_join(den, by = c("ccode_low" = "ccode", "year")) %>% rename(tl = tt) %>% left_join(den, by = c("ccode_high" = "ccode", "year")) %>% rename(th = tt) %>%
    mutate(ov = dyad > pmin(tl, th, na.rm = TRUE), sl = ifelse(ov, NA, dyad / tl), sh = ifelse(ov, NA, dyad / th)) %>% select(year, ccode_low, ccode_high, sl, sh)
  t1 <- tie0 %>% select(-trade_share_low, -trade_share_high) %>% inner_join(alt, by = c("year", "ccode_low", "ccode_high")) %>% rename(trade_share_low = sl, trade_share_high = sh)
  tf <- tempfile(fileext = ".csv"); write_csv(t1, tf)
  e <- new.env(); eval(parse(text = gsub("data/tie_trade.csv", tf, src[1:(k - 1)], fixed = TRUE)), envir = e)
  x <- e$d %>% filter(m == 0) %>% mutate(onset_m = as.integer(m1 == 1), l_dur = log(dur_m))
  cat(sprintf("\n== %s: %d dyad-years with valid shares in the onset sample, %d onsets; by era (dyad-years/onsets): %s\n", nm, nrow(x), sum(x$onset_m),
      paste(sprintf("%s %d/%d", levels(era(x$year)), as.integer(table(era(x$year))), as.integer(tapply(x$onset_m, era(x$year), sum))), collapse = ", ")))
  for (sp in c("", "+ z_vol + decade")) {
    m <- glm(as.formula(paste("onset_m ~ z_mx * contig + z_mn + z_cap + l_dur + l_same + l_diff", sp)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m)
    cat(sprintf("  %-18s non-contiguous %s | interaction %s | contiguous %s\n", ifelse(sp == "", "baseline", "+ volume, decade"), ci(b["z_mx"], sqrt(V["z_mx", "z_mx"])), ci(b["z_mx:contig"], sqrt(V["z_mx:contig", "z_mx:contig"])),
      ci(b["z_mx"] + b["z_mx:contig"], sqrt(sum(V[c("z_mx", "z_mx:contig"), c("z_mx", "z_mx:contig")]))))) } }
sink()
