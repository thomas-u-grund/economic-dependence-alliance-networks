# 95: FULL PANEL: larger share x contiguity, dyad-clustered vs actor-overlap (dyadic-robust) SEs for the implied effects. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), form_a = as.integer(a1 == 1))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f] (p=%.3f)", est, est - 1.96 * se, est + 1.96 * se, 2 * pnorm(-abs(est / se)))
fit <- function(x, y, dur, extra) {
  x$l_dur <- log(x[[dur]]); m <- glm(as.formula(paste(y, "~ z_mx * contig + z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial)
  X <- model.matrix(m); keep <- !is.na(coef(m)); X <- X[, keep, drop = FALSE]; fit <- glm.fit(X, x[[y]], family = binomial()); mu <- fit$fitted.values; U <- X * (x[[y]] - mu)
  A <- solve(crossprod(X, X * (mu * (1 - mu)))); n <- nrow(X); b <- coef(fit); names(b) <- colnames(X)
  Sd <- rowsum(U, x$dyad); Gd <- nrow(Sd); Vd <- A %*% crossprod(Sd) %*% A * (Gd / (Gd - 1))
  st <- rbind(data.frame(row = seq_len(n), s = x$lo), data.frame(row = seq_len(n), s = x$hi)); Ss <- rowsum(U[st$row, , drop = FALSE], group = st$s); Gs <- length(unique(st$s))
  Vs <- A %*% (crossprod(Ss) - crossprod(Sd)) %*% A * (Gs / (Gs - 1)); list(b = b, Vd = Vd, Vs = Vs, n = n) }
lin <- function(r, w) { est <- sum(r$b[names(w)] * w); c(est, sqrt(drop(t(w) %*% r$Vd[names(w), names(w)] %*% w)), sqrt(max(drop(t(w) %*% r$Vs[names(w), names(w)] %*% w), 0))) }
sink("results/geography_actor_overlap_output.txt", split = TRUE)
cat("FULL PANEL, larger share x contiguity: dyad-clustered SE vs actor-overlap SE\n")
for (p in list(list("dispute onset", "onset_m", "dur_m", d %>% filter(m == 0)), list("alliance formation", "form_a", "dur_a", d %>% filter(a == 0, m == 0)))) for (sp in c("", "+ z_vol + decade")) {
  r <- fit(p[[4]], p[[2]], p[[3]], sp); cat(sprintf("\n%s, %s (N=%d, events=%d)\n", p[[1]], ifelse(sp == "", "baseline", "+ volume + decade FE"), r$n, sum(p[[4]][[p[[2]]]])))
  for (w in list(list("non-contiguous", c(z_mx = 1)), list("interaction", c(`z_mx:contig` = 1)), list("contiguous", c(z_mx = 1, `z_mx:contig` = 1)))) { e <- lin(r, w[[2]])
    cat(sprintf("  %-15s dyad-clustered %s | actor-overlap %s | SE ratio %.2f\n", w[[1]], ci(e[1], e[2]), ci(e[1], e[3]), e[3] / e[2])) } }
sink()
