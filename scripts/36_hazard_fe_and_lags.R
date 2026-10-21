# 36: Robustness of the dyad-year hazard (script 35): lagged predictors, bloc/decade FE, dyad FE.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(fixest) })

src <- readLines("scripts/35_dyad_year_hazard.R"); cut <- grep("^z <- function", src)[1] - 1
eval(parse(text = src[1:cut]))          # rebuilds panel `d` (no output sink yet)

# --- lagged trade/capability within dyad (consecutive years only) ---
tt2 <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), dep = pmax(trade_share_low, trade_share_high),
            mut = pmin(trade_share_low, trade_share_high), cap = cap_asymmetry)
lagv <- function(k) d %>% select(year, lo, hi) %>% mutate(year_k = year - k) %>%
  left_join(tt2, by = c("year_k" = "year", "lo", "hi")) %>% select(dep, mut, cap) %>% rename_with(~ paste0(., "_l", k))
d <- bind_cols(d, lagv(1), lagv(3), lagv(5))

# --- bloc: both members in same large (>=8 member) alliance at t; id of first such alliance ---
mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE)
big <- mem %>% group_by(version4id) %>% summarise(n = n_distinct(ccode)) %>% filter(n >= 8) %>% pull(version4id)
mb <- mem %>% filter(version4id %in% big) %>% transmute(id = version4id, ccode, y0 = mem_st_year,
        y1 = ifelse(is.na(mem_end_year) | mem_end_year == 0, 9999, mem_end_year))
d$bloc <- "none"
for (id in big) {
  m <- mb %>% filter(id == !!id)
  inb <- function(cc, yr) sapply(seq_along(cc), function(i) any(m$ccode == cc[i] & m$y0 <= yr[i] & m$y1 >= yr[i]))
  hit <- d$bloc == "none" & inb(d$lo, d$year) & inb(d$hi, d$year)
  d$bloc[hit] <- as.character(id)
}
cat("share of dyad-years in a large-bloc dyad:", mean(d$bloc != "none"), "\n")
d$decade <- floor(d$year / 10) * 10
d <- d %>% mutate(flip_i = as.integer(flip), gone_i = as.integer(gone), end_i = as.integer(end_))

prep <- function(s, sfx) {
  s <- s %>% mutate(dep_ = .data[[paste0("dep", sfx)]], mut_ = .data[[paste0("mut", sfx)]], cap_ = .data[[paste0("cap", sfx)]]) %>%
    filter(!is.na(dep_), !is.na(mut_), !is.na(cap_))
  s %>% mutate(z_dep = as.numeric(scale(dep_)), z_mut = as.numeric(scale(mut_)), z_cap = as.numeric(scale(cap_)),
               l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
}
sink("results/hazard_fe_and_lags_output.txt", split = TRUE)
base <- "z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal"
rep1 <- function(s, y, sfx, fe, label) {
  x <- prep(s, sfx)
  f <- as.formula(paste(y, "~", base, if (nzchar(fe)) paste("|", fe) else NULL))
  m <- tryCatch(feglm(f, data = x, family = "binomial", cluster = ~dyad), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(m)) { cat(sprintf("%-55s FAILED\n", label)); return(invisible()) }
  ct <- coeftable(m)
  cat(sprintf("%-55s b_dep=%6.3f se=%5.3f p=%6.4f | unbal=%6.3f bal=%6.3f | n=%d ev=%d\n", label,
      ct["z_dep",1], ct["z_dep",2], ct["z_dep",4], ct["l_unbal",1], ct["l_bal",1], nobs(m), sum(x[[y]][x$dyad %in% unique(x$dyad)] )))
}
for (sg in c(-1, 1)) {
  s <- d %>% filter(sign == sg); nm <- ifelse(sg == -1, "DISPUTES", "ALLIANCES")
  for (y in c("flip_i", "gone_i")) {
    cat("\n####", nm, "-", y, "####\n")
    rep1(s, y, "",    "",             "baseline (t)")
    rep1(s, y, "",    "decade",       "+ decade FE")
    rep1(s, y, "",    "bloc",         "+ bloc FE")
    rep1(s, y, "",    "bloc + decade","+ bloc + decade FE")
    rep1(s, y, "_l1", "",             "asymmetry lagged 1y")
    rep1(s, y, "_l3", "",             "asymmetry lagged 3y")
    rep1(s, y, "_l5", "",             "asymmetry lagged 5y")
    rep1(s, y, "_l1", "bloc + decade","lag 1y + bloc + decade FE")
    rep1(s, y, "",    "dyad",         "DYAD FE (within-dyad)")
    rep1(s, y, "",    "dyad + decade","DYAD + decade FE")
  }
}
sink()
