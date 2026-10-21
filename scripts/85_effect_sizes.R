# 85: Substantive effect sizes (predicted annual onset probabilities at percentiles of the larger share, by contiguity) and influence checks. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
MU <- mean(d$mx); SDV <- sd(d$mx)                                   # global scale used for z_mx
d <- d %>% mutate(onset_m = as.integer(m1 == 1), tri = (n_same + n_diff) >= 1)
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
ctrl <- "z_mn + z_cap + l_dur + l_same + l_diff + z_vol"
sink("results/effect_sizes_output.txt", split = TRUE)
for (pop in c("embedded pairs (T)", "full panel")) {
  x <- d %>% filter(m == 0); if (pop != "full panel") x <- x %>% filter(tri); x$l_dur <- log(x$dur_m)
  cat(sprintf("\n######## %s: %d dyad-years, %d onsets; larger share percentiles (raw share): p50 %.4f, p90 %.4f, p99 %.4f; SD %.4f ########\n", toupper(pop), nrow(x), sum(x$onset_m), quantile(x$mx, .5), quantile(x$mx, .9), quantile(x$mx, .99), SDV))
  m <- glm(as.formula(paste("onset_m ~ z_mx * contig +", ctrl)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  base <- c(z_mn = mean(x$z_mn), z_cap = mean(x$z_cap), l_dur = mean(x$l_dur), l_same = mean(x$l_same), l_diff = mean(x$l_diff), z_vol = mean(x$z_vol))
  pr <- function(share, ct) { z <- unname((share - MU) / SDV); v <- c(`(Intercept)` = 1, z_mx = z, contig = ct, base, `z_mx:contig` = z * ct)[names(b)]; eta <- sum(b * v); se <- sqrt(drop(t(v) %*% V %*% v))
    c(1000 * plogis(eta), 1000 * plogis(eta - 1.96 * se), 1000 * plogis(eta + 1.96 * se)) }
  cat("Predicted onsets per 1,000 dyad-years (other covariates at sample means; 95% CI, dyad-clustered)\n")
  for (ct in c(0, 1)) for (q in c(.5, .9, .99)) { r <- pr(quantile(x$mx, q), ct); cat(sprintf("  %-14s share at p%02d = %.4f : %6.2f [%6.2f, %6.2f]\n", ifelse(ct == 0, "non-contiguous", "contiguous"), q * 100, quantile(x$mx, q), r[1], r[2], r[3])) }
  cat("Odds ratios (non-contiguous | contiguous):\n")
  for (mv in list(c(.5, .9), c(.9, .99))) { for (ct in c(0, 1)) { z1 <- unname((quantile(x$mx, mv[1]) - MU) / SDV); z2 <- unname((quantile(x$mx, mv[2]) - MU) / SDV); w <- c(z_mx = z2 - z1); if (ct == 1) w <- c(w, `z_mx:contig` = z2 - z1)
      est <- sum(b[names(w)] * w); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); cat(sprintf("  p%d -> p%d, %s: OR %.2f [%.2f, %.2f]\n", mv[1] * 100, mv[2] * 100, ifelse(ct == 0, "non-contiguous", "contiguous"), exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se))) } }
  cat("Influence checks, larger-share coefficient (non-contiguous | interaction with contiguity), controls as above:\n")
  chk <- function(xx, label, var = "z_mx") { mm <- glm(as.formula(paste("onset_m ~", var, "* contig +", ctrl)), xx, family = binomial); VV <- vcovCL(mm, cluster = xx$dyad); bb <- coef(mm); ok2 <- names(bb)[!is.na(bb)]; VV <- VV[ok2, ok2]; bb <- bb[ok2]
    a <- c(bb[var], sqrt(VV[var, var])); it <- paste0(var, ":contig"); cat(sprintf("  %-44s %s | interaction %s | onsets %d\n", label, ci(a[1], a[2]), ci(bb[it], sqrt(VV[it, it])), sum(xx$onset_m))) }
  chk(x, "all observations")
  chk(x %>% filter(mx <= quantile(mx, .99)), "drop dyad-years above the 99th percentile")
  top20 <- x %>% group_by(dyad) %>% summarise(mm = max(mx), .groups = "drop") %>% arrange(desc(mm)) %>% slice(1:20) %>% pull(dyad); chk(x %>% filter(!dyad %in% top20), "drop the 20 most exposed dyads")
  x$z_lmx <- as.numeric(scale(log(x$mx + 1e-4))); chk(x, "log(larger share + 0.0001), standardised", "z_lmx")
  x$z_rk <- as.numeric(scale(rank(x$mx) / nrow(x))); chk(x, "percentile rank of the larger share", "z_rk") }
sink()
