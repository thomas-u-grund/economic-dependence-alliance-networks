# 37: On the dyad-year hazard (script 35): (a) out-of-sample predictive gain from asymmetry,
#     folds grouped by dyad; (b) heterogeneity by trade diversification of the more exposed state.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
set.seed(20261021)
src <- readLines("scripts/35_dyad_year_hazard.R"); eval(parse(text = src[1:(grep("^z <- function", src)[1] - 1)]))
d <- d %>% select(-any_of(c("shl", "shh"))) %>% mutate(flip_i = as.integer(flip), gone_i = as.integer(gone))

# --- diversification of the more exposed state ---
tie <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), shl = trade_share_low, shh = trade_share_high)
hhi <- bind_rows(tie %>% transmute(year, cc = lo, sh = shl), tie %>% transmute(year, cc = hi, sh = shh)) %>%
  filter(!is.na(sh)) %>% group_by(year, cc) %>% summarise(hhi = sum(sh^2), .groups = "drop")
d <- d %>% left_join(tie, by = c("year", "lo", "hi")) %>%
  mutate(exp_cc = ifelse(!is.na(shl) & !is.na(shh) & shl >= shh, lo, hi), sh_exp = pmax(shl, shh)) %>%
  left_join(hhi, by = c("year", "exp_cc" = "cc")) %>% mutate(hhio = pmax(hhi - sh_exp^2, 0))

prep <- function(s) s %>% filter(!is.na(dep), !is.na(mut), !is.na(cap), !is.na(hhi)) %>%
  mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)),
         l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal),
         z_hhi = as.numeric(scale(hhi)), z_hhio = as.numeric(scale(hhio)))
auc <- function(p, y) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
base <- "z_mut + z_cap + l_dur + l_unbal + l_bal"
sink("results/hazard_predictive_and_diversification_output.txt", split = TRUE)
for (sg in c(-1, 1)) for (y in c("flip_i", "gone_i")) {
  x <- prep(d %>% filter(sign == sg)); nm <- paste(ifelse(sg == -1, "DISPUTES", "ALLIANCES"), y)
  cat("\n########", nm, " N =", nrow(x), " events =", sum(x[[y]]), " dyads =", n_distinct(x$dyad), "########\n")
  # (a) predictive gain, 10-fold CV grouped by dyad, repeated 5 times
  res <- replicate(5, {
    g <- unique(x$dyad); f <- setNames(sample(rep(1:10, length.out = length(g))), g); fo <- f[x$dyad]
    p1 <- p0 <- rep(NA_real_, nrow(x))
    for (k in 1:10) {
      tr <- x[fo != k, ]; te <- x[fo == k, ]
      m1 <- glm(as.formula(paste(y, "~ z_dep +", base)), tr, family = binomial)
      m0 <- glm(as.formula(paste(y, "~", base)), tr, family = binomial)
      p1[fo == k] <- predict(m1, te, type = "response"); p0[fo == k] <- predict(m0, te, type = "response")
    }
    yy <- x[[y]]; ll <- function(p) mean(yy * log(p) + (1 - yy) * log(1 - p))
    c(ll_full = ll(p1), ll_null = ll(p0), auc_full = auc(p1, yy), auc_null = auc(p0, yy))
  })
  r <- rowMeans(res)
  cat(sprintf("(a) held-out log-lik: full %.5f vs without asymmetry %.5f (gain %.5f, %.2f%% of null)\n    AUC: full %.4f vs without %.4f (gain %.4f)\n",
      r[1], r[2], r[1] - r[2], 100 * (r[1] - r[2]) / abs(r[2]), r[3], r[4], r[3] - r[4]))
  # in-sample LR test for comparison
  m1 <- glm(as.formula(paste(y, "~ z_dep +", base)), x, family = binomial); m0 <- update(m1, . ~ . - z_dep)
  cat(sprintf("    in-sample LR chi2(1) = %.1f\n", 2 * (logLik(m1) - logLik(m0))))
  # (b) diversification
  for (v in c("hhi", "hhio")) {
    xx <- x; xx$hi_ <- as.numeric(xx[[v]] > median(xx[[v]])); xx$zv <- xx[[paste0("z_", v)]]
    sub <- function(h) { s <- xx[xx$hi_ == h, ]; m <- glm(as.formula(paste(y, "~ z_dep +", base)), s, family = binomial)
      ct <- coeftest(m, vcov = vcovCL(m, cluster = s$dyad)); sprintf("b=%6.3f (se %.3f, p=%.4f, ev=%d)", ct["z_dep",1], ct["z_dep",2], ct["z_dep",4], sum(s[[y]])) }
    mi <- glm(as.formula(paste(y, "~ z_dep * zv +", base)), xx, family = binomial)
    ci <- coeftest(mi, vcov = vcovCL(mi, cluster = xx$dyad))
    cat(sprintf("(b) %-5s concentrated: %s\n          diversified:  %s\n          interaction z_dep x %s: b=%6.3f (p=%.4f)\n",
        v, sub(1), sub(0), v, ci["z_dep:zv",1], ci["z_dep:zv",4]))
  }
}
sink()
