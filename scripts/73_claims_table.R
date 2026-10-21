# 73: Claims table for 'what economic dependence does and does not do'. One row per outcome x specification; larger trade share (per SD). Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(fixest) })

s72 <- readLines("scripts/72_negative_tie_resolution.R"); k72 <- grep("^sink\\(", s72)[1]
eval(parse(text = s72[1:(k72 - 1)]))      # builds d (full panel incl. l_same, l_diff, decade...), dm, S (script-72 sample)
BOUND <- 0.15
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), tri = (n_same + n_diff) >= 1, form_a = as.integer(a1 == 1), end_a = as.integer(a1 == 0))
rows <- list()
add <- function(outcome, spec, est, se, n, ev, dyads, sc = NA, sd_ = NA) {
  lo <- est - 1.96 * se; hi <- est + 1.96 * se; lo90 <- est - 1.645 * se; hi90 <- est + 1.645 * se
  eqv <- ifelse(lo > 0 | hi < 0, "different from 0", ifelse(lo90 > -BOUND & hi90 < BOUND, "no effect beyond +-0.15", "inconclusive"))
  rows[[length(rows) + 1]] <<- data.frame(outcome, spec, N = n, events = ev, dyads, est = est, se = se, ci_lo = lo, ci_hi = hi, MDE = 2.8 * se, verdict = eqv, split_std = sc, common_std = sd_)
}
fitT <- function(x, y, outcome, spec, extra = "", dur = "dur_m") {
  x$l_dur <- log(x[[dur]]); x$zs_diff <- as.numeric(scale(x$l_diff)); x$zs_same <- as.numeric(scale(x$l_same))
  m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur + zs_same + zs_diff", extra)), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  add(outcome, spec, ct["z_mx", 1], ct["z_mx", 2], nobs(m), sum(x[[y]]), n_distinct(x$dyad[x[[y]] == 1]), ct["zs_diff", 1], ct["zs_same", 1])
}
pops <- list(
  list("dispute onset", "onset_m", d %>% filter(tri, m == 0), "dur_m"), list("dispute ending", "end_m", d %>% filter(tri, m == 1), "dur_m"),
  list("alliance formation", "form_a", d %>% filter(tri, a == 0, m == 0), "dur_a"), list("alliance ending", "end_a", d %>% filter(tri, a == 1, m == 0), "dur_a"))
for (p in pops) { fitT(p[[3]], p[[2]], p[[1]], "baseline", dur = p[[4]]); fitT(p[[3]], p[[2]], p[[1]], "+ volume + decade FE", "+ z_vol + decade", p[[4]])
  cx <- p[[3]] %>% filter(contig == 1); if (sum(cx[[p[[2]]]]) >= 30) fitT(cx, p[[2]], p[[1]], "contiguous only", "", p[[4]]) }
# dispute ending inside unbalanced (+,+,-) triads (script 72, spec D and with pact FE)
S$zs_diff <- as.numeric(scale(S$l_diff)); S$zs_pp <- as.numeric(scale(S$l_pp))
for (fe in c("", "| pact")) {
  m <- feglm(as.formula(paste("end_m ~ l_pp + l_diff + l_nn + l_dur + z_mx + z_mn + z_cap + z_vol + contig + z_lcap + allied", fe)), data = S, family = "binomial", cluster = ~dyad); ct <- coeftest(m)
  add("dispute ending, unbalanced (+,+,-) triads", ifelse(fe == "", "spec D (controls, allied)", "spec D + pact FE"), ct["z_mx", 1], ct["z_mx", 2], nobs(m), sum(S$end_m), n_distinct(S$dyad), NA, NA) }
tab <- bind_rows(rows) %>% mutate(across(c(est, se, ci_lo, ci_hi, MDE, split_std, common_std), ~ round(.x, 3)))
sink("results/claims_table_output.txt", split = TRUE)
cat(sprintf("Claims table: larger trade share, log-odds per SD; equivalence bound +-%.2f; dyad-clustered SEs. split_std/common_std = standardised coefficient of log1p split/common partners in the same model.\n\n", BOUND))
print(as.data.frame(tab %>% select(outcome, spec, N, events, dyads, est, ci_lo, ci_hi, MDE, verdict)), row.names = FALSE)
cat("\nStructure benchmark (standardised, same models):\n"); print(as.data.frame(tab %>% filter(!is.na(split_std)) %>% select(outcome, spec, split_std, common_std)), row.names = FALSE)
sink(); write_csv(tab, "results/claims_table.csv")
