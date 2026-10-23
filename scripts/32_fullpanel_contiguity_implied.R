# 32: FULL PANEL: larger share x contiguity with implied effects (non-contiguous, interaction, contiguous) for onset, ending, alliance formation. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), form_a = as.integer(a1 == 1))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
lin <- function(b, V, w) { est <- sum(b[names(w)] * w); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); c(est, se) }
sink("results/fullpanel_contiguity_implied_output.txt", split = TRUE)
for (p in list(list("dispute onset", "onset_m", "dur_m", d %>% filter(m == 0)), list("dispute ending", "end_m", "dur_m", d %>% filter(m == 1)), list("alliance formation", "form_a", "dur_a", d %>% filter(a == 0, m == 0)))) for (sp in c("", "+ z_vol + decade")) {
  x <- p[[4]]; x$l_dur <- log(x[[p[[3]]]]); m <- glm(as.formula(paste(p[[2]], "~ z_mx * contig + z_mn + z_cap + l_dur + l_same + l_diff", sp)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  a <- lin(b, V, c(z_mx = 1)); i <- lin(b, V, c(`z_mx:contig` = 1)); c1 <- lin(b, V, c(z_mx = 1, `z_mx:contig` = 1))
  cat(sprintf("%-20s %-18s non-contiguous %s | interaction %s | contiguous %s | N=%d events=%d (contiguous rows %d)\n", p[[1]], ifelse(sp == "", "baseline", "+ vol + decade FE"), ci(a[1], a[2]), ci(i[1], i[2]), ci(c1[1], c1[2]), nrow(x), sum(x[[p[[2]]]]), sum(x$contig))) }
sink()
