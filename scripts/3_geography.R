# 3_geography: Section 5.3 and Table 1: the geographical interaction and its robustness (Appendix B, Tables 4 to 7, and the calendar-year, spline and minimum-trade checks).
# Run from the repository root: Rscript scripts/3_geography.R (or run everything with Rscript run_all.R).
# Each block below is self-contained and writes the results/ file named in its sink() call.
source("scripts/setup.R")

# ---- FULL PANEL: larger share x contiguity with implied effects (non-contiguous, interaction, contiguous) for onset, ending, alliance formation.
local({
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), form_a = as.integer(a1 == 1))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
lin <- function(b, V, w) { est <- sum(b[names(w)] * w); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); c(est, se) }
sink("results/fullpanel_contiguity_implied_output.txt", split = TRUE)
for (p in list(list("dispute onset", "onset_m", "dur_m", d %>% filter(m == 0)), list("dispute ending", "end_m", "dur_m", d %>% filter(m == 1)), list("alliance formation", "form_a", "dur_a", d %>% filter(a == 0, m == 0)))) for (sp in c("", "+ z_vol + decade")) {
  x <- p[[4]]; x$l_dur <- log(x[[p[[3]]]]); m <- glm(as.formula(paste(p[[2]], "~ z_mx * contig + z_mn + z_cap + l_dur + l_same + l_diff", sp)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  a <- lin(b, V, c(z_mx = 1)); i <- lin(b, V, c(`z_mx:contig` = 1)); c1 <- lin(b, V, c(z_mx = 1, `z_mx:contig` = 1))
  cat(sprintf("%-20s %-18s non-contiguous %s | interaction %s | contiguous %s | N=%d events=%d (contiguous rows %d)\n", p[[1]], ifelse(sp == "", "baseline", "+ vol + decade FE"), ci(a[1], a[2]), ci(i[1], i[2]), ci(c1[1], c1[2]), nrow(x), sum(x[[p[[2]]]]), sum(x$contig))) }
sink()
})

# ---- FULL PANEL: larger share x contiguity, dyad-clustered vs actor-overlap (dyadic-robust) SEs for the implied effects.
local({
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
})

# ---- Exposure x land/water contiguity, exposure x embedding (full panel), and exposure x log capital distance, each with explicit interaction CIs.
local({
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), tri = (n_same + n_diff) >= 1, T = as.integer(tri))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
lin <- function(b, V, w) { est <- sum(b[names(w)] * w); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); c(est, se) }
fit <- function(x, y, dur, mod, label, extra = "") {
  x$l_dur <- log(x[[dur]]); f <- as.formula(paste(y, "~ z_mx *", mod, "+ z_mn + z_cap + l_dur + l_same + l_diff", if (mod != "contig") "+ contig" else "", extra))
  m <- glm(f, x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  inter <- paste0("z_mx:", mod); main <- lin(b, V, c(z_mx = 1)); it <- lin(b, V, setNames(1, inter)); both <- lin(b, V, setNames(c(1, 1), c("z_mx", inter)))
  cat(sprintf("%-52s base (%s=0) %s | interaction %s | implied at %s=1 %s | N=%d events=%d\n", label, mod, ci(main[1], main[2]), ci(it[1], it[2]), mod, ci(both[1], both[2]), nrow(x), sum(x[[y]])))
}
# capital distance for dyad-years of the embedded onset/ending samples (1886+)
yrs <- sort(unique(d$year[d$tri & d$year >= 1886]))
cap <- bind_rows(lapply(yrs, function(y) { mm <- distmatrix(as.Date(paste0(y, "-01-01")), type = "capdist", useGW = FALSE); ids <- as.numeric(rownames(mm)); ix <- which(upper.tri(mm), arr.ind = TRUE)
  tibble(year = y, lo = pmin(ids[ix[, 1]], ids[ix[, 2]]), hi = pmax(ids[ix[, 1]], ids[ix[, 2]]), capdist = mm[ix]) }))
d <- d %>% left_join(cap, by = c("year", "lo", "hi")) %>% mutate(z_lcap = as.numeric(scale(log(capdist + 1))))
sink("results/heterogeneity_output.txt", split = TRUE)
cat("(a) Embedded pairs: larger share x contiguity (land or water, any type)\n")
for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE")
  fit(d %>% filter(tri, m == 0), "onset_m", "dur_m", "contig", paste("onset,", s), sp); fit(d %>% filter(tri, m == 1), "end_m", "dur_m", "contig", paste("ending,", s), sp) }
cat("\n(b) Full panel: larger share x embedded pair (T)\n")
for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE")
  fit(d %>% filter(m == 0), "onset_m", "dur_m", "T", paste("onset,", s), sp); fit(d %>% filter(m == 1), "end_m", "dur_m", "T", paste("ending,", s), sp) }
cat("\n(c) Embedded pairs, 1886+: larger share x standardised log capital distance (controls include contiguity)\n")
x <- d %>% filter(tri, !is.na(capdist))
for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE"); fit(x %>% filter(m == 0), "onset_m", "dur_m", "z_lcap", paste("onset,", s), sp); fit(x %>% filter(m == 1), "end_m", "dur_m", "z_lcap", paste("ending,", s), sp) }
sink()
})

# ---- Who are the distant exposed pairs? Interactions of the larger share with major-power involvement and same-region status among non-contiguous pairs.
local({
nmc <- read_csv("raw/nmc/NMCv7/abridged/NMC-70-abridged.csv", show_col_types = FALSE) %>% select(ccode, year, cinc) %>% group_by(year) %>% mutate(rk = rank(-cinc, ties.method = "min")) %>% ungroup()
MAJ <- c(2, 200, 220, 255, 300, 325, 365, 710, 740)
region <- function(cc) cut(cc, c(0, 199, 399, 626, 699, 999), labels = c("Americas", "Europe", "Africa", "MiddleEast", "Asia"))
d <- d %>% left_join(nmc %>% transmute(year, lo = ccode, rk_lo = rk), by = c("year", "lo")) %>% left_join(nmc %>% transmute(year, hi = ccode, rk_hi = rk), by = c("year", "hi")) %>%
  mutate(top5 = as.integer(coalesce(rk_lo <= 5, FALSE) | coalesce(rk_hi <= 5, FALSE)), ever_major = as.integer(lo %in% MAJ | hi %in% MAJ), same_region = as.integer(region(lo) == region(hi)),
         onset_m = as.integer(m1 == 1), form_a = as.integer(a1 == 1), tri = (n_same + n_diff) >= 1) %>% filter(contig == 0)
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
lin <- function(b, V, w) { est <- sum(b[names(w)] * w); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); c(est, se) }
fit <- function(x, y, dur, mod, label, extra = "") {
  x$l_dur <- log(x[[dur]]); m <- glm(as.formula(paste(y, "~ z_mx *", mod, "+ z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial)
  V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]; inter <- paste0("z_mx:", mod)
  a <- lin(b, V, c(z_mx = 1)); i <- lin(b, V, setNames(1, inter)); c1 <- lin(b, V, setNames(c(1, 1), c("z_mx", inter)))
  cat(sprintf("%-44s mod=0 %s | interaction %s | mod=1 %s | N=%d events=%d (mod=1: %d rows)\n", label, ci(a[1], a[2]), ci(i[1], i[2]), ci(c1[1], c1[2]), nrow(x), sum(x[[y]]), sum(x[[mod]])))
}
sub <- function(x, y, dur, cond, label, extra = "") {
  x <- x %>% filter(eval(parse(text = cond))); x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%-44s larger share %s | N=%d events=%d\n", label, ci(ct["z_mx", 1], ct["z_mx", 2]), nrow(x), sum(x[[y]])))
}
sink("results/distant_pairs_output.txt", split = TRUE)
for (pop in c("embedded pairs (T)", "full panel")) {
  cat(sprintf("\n######## NON-CONTIGUOUS %s ########\n", toupper(pop)))
  dd <- if (pop == "full panel") d else d %>% filter(tri)
  for (oc in list(list("onset", "onset_m", "dur_m", dd %>% filter(m == 0)), list("alliance formation", "form_a", "dur_a", dd %>% filter(a == 0, m == 0)))) {
    cat(sprintf("\n== %s ==\n", oc[[1]]))
    for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE")
      for (mod in c("top5", "ever_major", "same_region")) fit(oc[[4]], oc[[2]], oc[[3]], mod, paste(mod, s), sp)
      sub(oc[[4]], oc[[2]], oc[[3]], "top5 == 0", paste("no top-5 state,", s), sp); sub(oc[[4]], oc[[2]], oc[[3]], "ever_major == 0", paste("no major power,", s), sp); sub(oc[[4]], oc[[2]], oc[[3]], "TRUE", paste("all non-contiguous,", s), sp) } } }
cat("\nShare of non-contiguous embedded onset risk-set rows with a top-5 state:", round(mean(d$top5[d$tri & d$m == 0]), 3), "; ever_major:", round(mean(d$ever_major[d$tri & d$m == 0]), 3), "; same region:", round(mean(d$same_region[d$tri & d$m == 0]), 3), "\n")
sink()
})

# ---- Substantive effect sizes (predicted annual onset probabilities at percentiles of the larger share, by contiguity) and influence checks.
local({
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
})

# ---- FULL PANEL, dispute onset: larger share x contiguity without the top 1% of exposure, with calendar-year effects, and with flexible forms.
local({
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
})

# ---- FULL PANEL, dispute onset: sensitivity of the geographical interaction to the minimum-trade rule for the denominator of the trade shares.
#     Rules: (a) total trade >= US$50 million nominal (the rule used in earlier versions), (b) no minimum (positive total trade), (c) total trade >= 0.01% (the article's rule, scripts/build_data.R) and (d) >= 0.05% of recorded world trade that year.
local({
nat <- read_csv("raw/trade/COW_Trade_4.0/National_COW_4.0.csv", show_col_types = FALSE) %>%
  mutate(tt = imports + exports, tt = ifelse(tt <= 0, NA, tt)) %>% select(ccode, year, tt) %>% filter(!is.na(tt))
nat <- nat %>% group_by(year) %>% mutate(world = sum(tt), rel = tt / world) %>% ungroup()
rules <- list("US$50m nominal (earlier rule)" = function(x) x$tt >= 50, "no minimum" = function(x) x$tt > 0,
              "0.01% of recorded world trade (article)" = function(x) x$rel >= 1e-4, "0.05% of recorded world trade" = function(x) x$rel >= 5e-4)
era <- function(y) cut(y, c(1869, 1913, 1945, 1989, 2011), labels = c("1870-1913", "1914-45", "1946-89", "1990-2010"))
sink("results/trade_threshold_output.txt", split = TRUE)
cat("State-years with recorded total trade, and share meeting each rule, by era\n")
for (nm in names(rules)) { k <- nat %>% mutate(ok = rules[[nm]](cur_data()), era = era(year)) %>% group_by(era) %>% summarise(n = n(), kept = mean(ok)); cat(sprintf("%-32s %s\n", nm, paste(sprintf("%s: %d (%.0f%%)", k$era, k$n, 100 * k$kept), collapse = " | "))) }
tie0 <- read_csv("data/tie_trade.csv.gz", show_col_types = FALSE)
dy <- read_csv("raw/trade/COW_Trade_4.0/Dyadic_COW_4.0.csv", show_col_types = FALSE) %>% mutate(across(c(flow1, flow2, smoothtotrade), ~ na_if(., -9))) %>%
  mutate(ccode_low = pmin(ccode1, ccode2), ccode_high = pmax(ccode1, ccode2)) %>% group_by(year, ccode_low, ccode_high) %>% summarise(dyad = sum(smoothtotrade, na.rm = TRUE), .groups = "drop")
src <- readLines("scripts/panel.R")
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
for (nm in names(rules)) {
  den <- nat %>% filter(rules[[nm]](.)) %>% select(ccode, year, tt)
  alt <- dy %>% left_join(den, by = c("ccode_low" = "ccode", "year")) %>% rename(tl = tt) %>% left_join(den, by = c("ccode_high" = "ccode", "year")) %>% rename(th = tt) %>%
    mutate(ov = dyad > pmin(tl, th, na.rm = TRUE), sl = ifelse(ov, NA, dyad / tl), sh = ifelse(ov, NA, dyad / th)) %>% select(year, ccode_low, ccode_high, sl, sh)
  t1 <- tie0 %>% select(-trade_share_low, -trade_share_high) %>% inner_join(alt, by = c("year", "ccode_low", "ccode_high")) %>% rename(trade_share_low = sl, trade_share_high = sh)
  tf <- tempfile(fileext = ".csv"); write_csv(t1, tf)
  e <- new.env(); eval(parse(text = gsub("data/tie_trade.csv.gz", tf, src, fixed = TRUE)), envir = e)
  x <- e$d %>% filter(m == 0) %>% mutate(onset_m = as.integer(m1 == 1), l_dur = log(dur_m))
  cat(sprintf("\n== %s: %d dyad-years with valid shares in the onset sample, %d onsets; by era (dyad-years/onsets): %s\n", nm, nrow(x), sum(x$onset_m),
      paste(sprintf("%s %d/%d", levels(era(x$year)), as.integer(table(era(x$year))), as.integer(tapply(x$onset_m, era(x$year), sum))), collapse = ", ")))
  for (sp in c("", "+ z_vol + decade")) {
    m <- glm(as.formula(paste("onset_m ~ z_mx * contig + z_mn + z_cap + l_dur + l_same + l_diff", sp)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m)
    cat(sprintf("  %-18s non-contiguous %s | interaction %s | contiguous %s\n", ifelse(sp == "", "baseline", "+ volume, decade"), ci(b["z_mx"], sqrt(V["z_mx", "z_mx"])), ci(b["z_mx:contig"], sqrt(V["z_mx:contig", "z_mx:contig"])),
      ci(b["z_mx"] + b["z_mx:contig"], sqrt(sum(V[c("z_mx", "z_mx:contig"), c("z_mx", "z_mx:contig")]))))) } }
sink()
})
