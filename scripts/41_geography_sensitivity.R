# 41: FULL PANEL, dispute onset: larger share x contiguity without the top 1% of exposure, with calendar-year effects, and with flexible forms. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(splines) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
x <- d %>% filter(m == 0) %>% mutate(onset_m = as.integer(m1 == 1), l_dur = log(dur_m), yr = factor(year), lmx = log(mx + 0.001), lmn = log(mn + 0.001))
cut99 <- quantile(x$mx, 0.99)
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f] (p=%.3f)", est, est - 1.96 * se, est + 1.96 * se, 2 * pnorm(-abs(est / se)))
rep_lin <- function(m, dat, label) {
  b <- coef(m); V <- vcovCL(m, cluster = dat$dyad); n1 <- "z_mx"; n2 <- "z_mx:contig"
  w <- c(1, 1); est <- b[n1] + b[n2]; se <- sqrt(sum(V[c(n1, n2), c(n1, n2)]))
  cat(sprintf("%-46s N=%d ev=%d | non-contig %s | interaction %s | contig %s\n", label, nrow(dat), sum(dat$onset_m), ci(b[n1], sqrt(V[n1, n1])), ci(b[n2], sqrt(V[n2, n2])), ci(est, se))) }
f <- "onset_m ~ z_mx * contig + z_mn + z_cap + l_dur + l_same + l_diff"
sink("results/geography_sensitivity_output.txt", split = TRUE)
cat(sprintf("Onset sample: %d dyad-years, %d onsets. 99th percentile of larger share = %.4f (%d observations above)\n", nrow(x), sum(x$onset_m), cut99, sum(x$mx > cut99)))
for (sp in c("", "+ z_vol + decade")) {
  tag <- ifelse(sp == "", "baseline", "+ volume + decade FE")
  rep_lin(glm(as.formula(paste(f, sp)), x, family = binomial), x, paste("all observations,", tag))
  xt <- x %>% filter(mx <= cut99); rep_lin(glm(as.formula(paste(f, sp)), xt, family = binomial), xt, paste("excluding top 1% of larger share,", tag))
  xt5 <- x %>% filter(mx <= quantile(x$mx, 0.95)); rep_lin(glm(as.formula(paste(f, sp)), xt5, family = binomial), xt5, paste("excluding top 5% of larger share,", tag)) }
cat("\nCalendar-year fixed effects instead of decade effects (+ volume)\n")
m <- glm(as.formula(paste(f, "+ z_vol + yr")), x, family = binomial); rep_lin(m, x, "calendar-year FE + volume")
cat("\nFlexible form (+ volume + decade FE): log share, and natural spline in the larger share by contiguity\n")
m <- glm(onset_m ~ lmx * contig + lmn + z_cap + l_dur + l_same + l_diff + z_vol + decade, x, family = binomial); ct <- coeftest(m, vcovCL(m, cluster = x$dyad))
cat(sprintf("log larger share: non-contig %s | interaction %s\n", ci(ct["lmx", 1], ct["lmx", 2]), ci(ct["lmx:contig", 1], ct["lmx:contig", 2])))
ms <- glm(onset_m ~ ns(lmx, 3) * contig + z_mn + z_cap + l_dur + l_same + l_diff + z_vol + decade, x, family = binomial)
mb <- update(ms, . ~ . - ns(lmx, 3):contig); lr <- anova(mb, ms, test = "LRT"); cat(sprintf("spline x contiguity LR test: chi2 = %.2f, df = %d, p = %.4g\n", lr$Deviance[2], lr$Df[2], lr$`Pr(>Chi)`[2]))
base <- x %>% summarise(across(c(z_mn, z_cap, l_dur, l_same, l_diff, z_vol), median))
for (cg in 0:1) { cat(sprintf("Predicted onsets per 1,000 dyad-years, contig=%d (covariates at medians, decade 1980s):\n", cg))
  for (q in c(.5, .9, .95, .99)) { v <- quantile(x$mx[x$contig == cg], q); nd <- cbind(base, lmx = log(v + 0.001), lmn = median(x$lmn), contig = cg, decade = factor(198, levels = levels(x$decade)))
    cat(sprintf("  larger share at p%02d (%.4f): %.2f\n", q * 100, v, 1000 * predict(ms, nd, type = "response"))) } }
sink()
