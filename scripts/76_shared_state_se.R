# 76: Dyad-clustered vs actor-overlap (dyadic-robust) SEs for the larger-share coefficient in the embedded-pair panel. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })

src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), tri = (n_same + n_diff) >= 1, form_a = as.integer(a1 == 1), end_a = as.integer(a1 == 0))
robust_se <- function(x, y, extra, dur) {
  x$l_dur <- log(x[[dur]]); f <- as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra))
  m <- glm(f, x, family = binomial); X <- model.matrix(m); keep <- !is.na(coef(m)); X <- X[, keep, drop = FALSE]
  fit <- glm.fit(X, x[[y]], family = binomial()); mu <- fit$fitted.values; U <- X * (x[[y]] - mu); n <- nrow(X)
  A <- solve(crossprod(X, X * (mu * (1 - mu))))                       # inverse information
  Sd <- rowsum(U, x$dyad); Gd <- nrow(Sd); Vd <- A %*% crossprod(Sd) %*% A * (Gd / (Gd - 1))
  chk <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad)); j <- which(colnames(X) == "z_mx")
  stopifnot(abs(sqrt(Vd[j, j]) - chk["z_mx", 2]) / chk["z_mx", 2] < 1e-3)      # manual dyad SE reproduces sandwich
  st <- rbind(data.frame(row = seq_len(n), s = x$lo), data.frame(row = seq_len(n), s = x$hi))
  Ss <- rowsum(U[st$row, , drop = FALSE], group = st$s); Gs <- length(unique(st$s))
  meatD <- crossprod(Ss) - crossprod(Sd); Vs <- A %*% meatD %*% A * (Gs / (Gs - 1))
  list(est = unname(coef(fit)["z_mx"]), se_dyad = sqrt(Vd[j, j]), se_dyadic = sqrt(max(Vs[j, j], 0)), neg_var = Vs[j, j] < 0, N = n, ev = sum(x[[y]]), states = Gs)
}
pops <- list(list("dispute onset", "onset_m", d %>% filter(tri, m == 0), "dur_m"), list("dispute ending", "end_m", d %>% filter(tri, m == 1), "dur_m"),
             list("alliance formation", "form_a", d %>% filter(tri, a == 0, m == 0), "dur_a"), list("alliance ending", "end_a", d %>% filter(tri, a == 1, m == 0), "dur_a"))
res <- list()
for (p in pops) for (sp in c("baseline", "+ volume + decade FE")) {
  r <- robust_se(p[[3]], p[[2]], if (sp == "baseline") "" else "+ z_vol + decade", p[[4]])
  res[[length(res) + 1]] <- data.frame(outcome = p[[1]], spec = sp, N = r$N, events = r$ev, est = r$est, se_dyad = r$se_dyad, se_dyadic = r$se_dyadic, ratio = r$se_dyadic / r$se_dyad,
    p_dyad = 2 * pnorm(-abs(r$est / r$se_dyad)), p_dyadic = 2 * pnorm(-abs(r$est / r$se_dyadic)), ci_lo_dyadic = r$est - 1.96 * r$se_dyadic, ci_hi_dyadic = r$est + 1.96 * r$se_dyadic, negvar = r$neg_var) }
tab <- bind_rows(res) %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))
sink("results/shared_state_se_output.txt", split = TRUE); cat("Larger share, embedded pairs: dyad-clustered vs actor-overlap (dyadic-robust) SE\n"); print(as.data.frame(tab), row.names = FALSE); sink()
write_csv(tab, "results/shared_state_se.csv")
