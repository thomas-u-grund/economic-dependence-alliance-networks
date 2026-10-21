# 62: Standardised coefficients and single-predictor held-out AUC within allied dispute-years and within alliance-years.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
set.seed(20261026)
d0 <- build_panel(base_signed(), "baseline")
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>% transmute(year, lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2), allied = 1L) %>% distinct()
auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
prep <- function(x) x %>% filter(!is.na(dep), !is.na(mut), !is.na(cap)) %>% mutate(z_larger = as.numeric(scale(dep)), z_smaller = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), z_dur = as.numeric(scale(log(dur))), z_U = as.numeric(scale(log1p(n_unbal))), z_B = as.numeric(scale(log1p(n_bal))))
preds <- c("z_larger", "z_smaller", "z_cap", "z_dur", "z_U", "z_B")
run <- function(x, label) {
  cat(sprintf("\n######## %s: N=%d, events=%d, dyads=%d\n", label, nrow(x), sum(x$flip_i), n_distinct(x$dyad)))
  m <- glm(as.formula(paste("flip_i ~", paste(preds, collapse = "+"))), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  dy <- unique(x$dyad)
  for (v in preds) { a <- numeric(5); for (r in 1:5) { f <- sample(rep(1:10, length.out = length(dy))); names(f) <- dy; fo <- f[x$dyad]; p <- rep(NA_real_, nrow(x))
      for (k in 1:10) { mm <- glm(as.formula(paste("flip_i ~", v)), x[fo != k, ], family = binomial); p[fo == k] <- predict(mm, x[fo == k, ], type = "response") }; a[r] <- auc(x$flip_i, p) }
    cat(sprintf("%-10s joint coef (per SD) %6.3f (se %.3f, p=%.4f) | single-predictor held-out AUC %.3f\n", v, ct[v, 1], ct[v, 2], ct[v, 4], mean(a))) } }
x <- prep(d0 %>% filter(sign == -1) %>% left_join(al, by = c("year", "lo", "hi")) %>% mutate(allied = replace_na(allied, 0L)) %>% filter(allied == 1))
sink("results/strongest_predictor_output.txt", split = TRUE)
run(x, "ALLIED DISPUTE-YEARS (dispute ends inside the alliance)")
cat(sprintf("corr(log(1+U), larger share) = %.3f; share of U that is common allies is not separately measured here\n", cor(x$z_U, x$z_larger)))
run(prep(d0 %>% filter(sign == 1)), "ALLIANCE-YEARS (dispute begins inside the alliance)")
sink()
