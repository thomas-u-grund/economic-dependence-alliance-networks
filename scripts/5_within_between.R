# 5_within_between: Section 5.4 and Table 2: between- and within-pair variation, dyad fixed effects, longer-run exposure, embedded pairs (Appendix A, Table 3; Appendix B, Table 4; Appendix D, Table 14).
# Run from the repository root: Rscript scripts/5_within_between.R (or run everything with Rscript run_all.R).
# Each block below is self-contained and writes the results/ file named in its sink() call.
source("scripts/setup.R")

# ---- FULL PANEL: within-between decomposition, larger share x log capital distance, escalation given onset.
local({
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), form_a = as.integer(a1 == 1))
yrs <- sort(unique(d$year[d$year >= 1886]))
cap <- bind_rows(lapply(yrs, function(y) { mm <- distmatrix(as.Date(paste0(y, "-01-01")), type = "capdist", useGW = FALSE); ids <- as.numeric(rownames(mm)); ix <- which(upper.tri(mm), arr.ind = TRUE)
  tibble(year = y, lo = pmin(ids[ix[, 1]], ids[ix[, 2]]), hi = pmax(ids[ix[, 1]], ids[ix[, 2]]), capdist = mm[ix]) }))
d <- d %>% left_join(cap, by = c("year", "lo", "hi")) %>% mutate(z_lcap = as.numeric(scale(log(capdist + 1))))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
wb <- function(x, y, dur, label, extra = "") {
  x <- x %>% group_by(dyad) %>% mutate(mx_b = mean(z_mx), mx_w = z_mx - mx_b, n_dyad_rows = n()) %>% ungroup(); x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ mx_b + mx_w + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  dif <- b["mx_b"] - b["mx_w"]; sed <- sqrt(V["mx_b", "mx_b"] + V["mx_w", "mx_w"] - 2 * V["mx_b", "mx_w"])
  cat(sprintf("%-40s between %s | within %s | difference %s | N=%d events=%d dyads=%d (with outcome variation: %d)\n", label, ci(b["mx_b"], sqrt(V["mx_b", "mx_b"])), ci(b["mx_w"], sqrt(V["mx_w", "mx_w"])), ci(dif, sed), nrow(x), sum(x[[y]]), n_distinct(x$dyad),
      x %>% group_by(dyad) %>% filter(n_distinct(.data[[y]]) > 1) %>% pull(dyad) %>% n_distinct()))
}
lin <- function(b, V, w) { est <- sum(b[names(w)] * w); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); c(est, se) }
dist_int <- function(x, y, dur, label, extra = "") {
  x <- x %>% filter(!is.na(capdist)); x$l_dur <- log(x[[dur]]); m <- glm(as.formula(paste(y, "~ z_mx * z_lcap + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra)), x, family = binomial)
  V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]; a <- lin(b, V, c(z_mx = 1)); i <- lin(b, V, c(`z_mx:z_lcap` = 1)); c1 <- lin(b, V, c(z_mx = 1, `z_mx:z_lcap` = 1))
  cat(sprintf("%-40s at average distance %s | interaction %s | at +1 SD distance %s | N=%d events=%d\n", label, ci(a[1], a[2]), ci(i[1], i[2]), ci(c1[1], c1[2]), nrow(x), sum(x[[y]])))
}
sink("results/fullpanel_within_het_output.txt", split = TRUE)
cat("######## Within-between, full panel (larger share per SD, global scale) ########\n")
pops <- list(list("dispute onset", "onset_m", "dur_m", d %>% filter(m == 0)), list("dispute ending", "end_m", "dur_m", d %>% filter(m == 1)), list("alliance formation", "form_a", "dur_a", d %>% filter(a == 0, m == 0)))
for (p in pops) { wb(p[[4]], p[[2]], p[[3]], paste(p[[1]], "baseline")); wb(p[[4]], p[[2]], p[[3]], paste(p[[1]], "+ vol + decade FE"), "+ z_vol + decade") }
cat("\n######## Larger share x standardised log capital distance, full panel 1886+ (contiguity as control) ########\n")
for (p in pops[1:2]) for (sp in c("", "+ z_vol + decade")) dist_int(p[[4]], p[[2]], p[[3]], paste(p[[1]], ifelse(sp == "", "baseline", "+ vol + decade FE")), sp)
for (sp in c("", "+ z_vol + decade")) dist_int(pops[[3]][[4]], "form_a", "dur_a", paste("alliance formation", ifelse(sp == "", "baseline", "+ vol + decade FE")), sp)
cat("\n######## Escalation given onset, full panel ########\n")
mid <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb), fatlev = ifelse(fatlev < 0, NA, fatlev)) %>%
  group_by(year, lo, hi) %>% summarise(h1 = max(hihost), fat1 = suppressWarnings(max(fatlev, na.rm = TRUE)), .groups = "drop") %>% mutate(fat1 = ifelse(is.infinite(fat1), NA, fat1), year = year - 1)
o <- d %>% filter(m == 0, m1 == 1) %>% left_join(mid, by = c("year", "lo", "hi")) %>% filter(!is.na(h1)) %>% mutate(l_dur = log(dur_m), y_force = as.integer(h1 >= 4), y_fat = as.integer(!is.na(fat1) & fat1 >= 1))
cat(sprintf("Onsets with hostility information: %d; use of force or war: %d (%.1f%%); any fatalities: %d (%.1f%%)\n", nrow(o), sum(o$y_force), 100 * mean(o$y_force), sum(o$y_fat), 100 * mean(o$y_fat)))
for (sp in c("", "+ z_vol + decade")) { cat(sprintf("  spec: %s\n", ifelse(sp == "", "baseline", sp)))
  for (y in c("y_force", "y_fat")) { m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff", sp)), o, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = o$dyad)); cat(sprintf("    %-8s %s p=%.4f events %d of %d\n", y, ci(ct["z_mx", 1], ct["z_mx", 2]), ct["z_mx", 4], sum(o[[y]]), nrow(o))) } }
sink()
})

# ---- FULL PANEL: dyad fixed-effects logit for dispute onset and ending (uses only dyads with variation in the outcome).
local({
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), decade = as.factor(decade))
sink("results/fullpanel_dyad_fe_output.txt", split = TRUE)
for (spec in list(list("DISPUTE ONSET", "onset_m", d %>% filter(m == 0), "dur_m"), list("DISPUTE ENDING", "end_m", d %>% filter(m == 1), "dur_m"))) {
  x <- spec[[3]]; x$l_dur <- log(x[[spec[[4]]]]); cat(sprintf("\n## %s, full panel\n", spec[[1]]))
  m0 <- feglm(as.formula(paste(spec[[2]], "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff")), x, family = "binomial", cluster = ~dyad); c0 <- coeftable(m0)
  m1 <- feglm(as.formula(paste(spec[[2]], "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff | dyad")), x, family = "binomial", cluster = ~dyad); c1 <- coeftable(m1)
  m2 <- feglm(as.formula(paste(spec[[2]], "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff | dyad + decade")), x, family = "binomial", cluster = ~dyad); c2 <- coeftable(m2)
  cat(sprintf("pooled              larger %6.3f (se %.3f, p=%.4f)\n", c0["z_mx", 1], c0["z_mx", 2], c0["z_mx", 4]))
  cat(sprintf("FE: dyad            larger %6.3f (se %.3f, p=%.4f) | N used=%d\n", c1["z_mx", 1], c1["z_mx", 2], c1["z_mx", 4], m1$nobs))
  cat(sprintf("FE: dyad + decade   larger %6.3f (se %.3f, p=%.4f) | N used=%d\n", c2["z_mx", 1], c2["z_mx", 2], c2["z_mx", 4], m2$nobs)) }
sink()
})

# ---- FULL PANEL: within-between decomposition of the larger share separately for contiguous and non-contiguous dyads.
local({
d <- d %>% mutate(onset_m = as.integer(m1 == 1), form_a = as.integer(a1 == 1))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
wb <- function(x, y, dur, label, extra = "") {
  x <- x %>% group_by(dyad) %>% mutate(mx_b = mean(z_mx), mx_w = z_mx - mx_b) %>% ungroup(); x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ mx_b + mx_w + z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  dif <- b["mx_b"] - b["mx_w"]; sed <- sqrt(V["mx_b", "mx_b"] + V["mx_w", "mx_w"] - 2 * V["mx_b", "mx_w"])
  cat(sprintf("%-46s between %s | within %s | difference %s | N=%d events=%d dyads=%d\n", label, ci(b["mx_b"], sqrt(V["mx_b", "mx_b"])), ci(b["mx_w"], sqrt(V["mx_w", "mx_w"])), ci(dif, sed), nrow(x), sum(x[[y]]), n_distinct(x$dyad))) }
wbi <- function(x, y, dur, label, extra = "") {
  x <- x %>% group_by(dyad) %>% mutate(mx_b = mean(z_mx), mx_w = z_mx - mx_b) %>% ungroup(); x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ (mx_b + mx_w) * contig + z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  cat(sprintf("%-46s contiguity difference in between %s | in within %s\n", label, ci(b["mx_b:contig"], sqrt(V["mx_b:contig", "mx_b:contig"])), ci(b["mx_w:contig"], sqrt(V["mx_w:contig", "mx_w:contig"])))) }
sink("results/within_between_by_contiguity_output.txt", split = TRUE)
pops <- list(list("dispute onset", "onset_m", "dur_m", d %>% filter(m == 0)), list("alliance formation", "form_a", "dur_a", d %>% filter(a == 0, m == 0)))
for (p in pops) for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE")
  cat("\n")
  wb(p[[4]] %>% filter(contig == 0), p[[2]], p[[3]], paste(p[[1]], "non-contiguous,", s), sp)
  wb(p[[4]] %>% filter(contig == 1), p[[2]], p[[3]], paste(p[[1]], "contiguous,", s), sp)
  wbi(p[[4]], p[[2]], p[[3]], paste(p[[1]], "interaction model,", s), sp) }
cat("\nDyad fixed-effects logit among non-contiguous dyads, larger share (per SD)\n")
for (p in pops) for (sp in c("", "+ z_vol")) { x <- p[[4]] %>% filter(contig == 0); x$l_dur <- log(x[[p[[3]]]])
  m <- feglm(as.formula(paste(p[[2]], "~ z_mx + z_mn + z_cap + l_dur + l_same + l_diff", sp, "| dyad")), x, family = "binomial", cluster = ~dyad)
  cat(sprintf("%-20s %-10s %s (p=%.3f) | dyads used %d, obs %d\n", p[[1]], ifelse(sp == "", "baseline", sp), ci(coef(m)["z_mx"], se(m)["z_mx"]), pvalue(m)["z_mx"], length(unique(x$dyad[x[[p[[2]]]] == 1])), nobs(m))) }
sink()
})

# ---- FULL PANEL: longer-run within-pair exposure with calendar-year windows (dyad fixed-effects logit).
#      Predictors are measured in year t (or earlier) and the outcome is the dispute beginning in t+1.
local({
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
d <- d %>% arrange(dyad, year) %>% mutate(onset_m = as.integer(m1 == 1))
win <- function(yrs, z, lo, hi, minn = 3) vapply(seq_along(yrs), function(i) { k <- yrs >= yrs[i] - lo & yrs <= yrs[i] - hi; if (sum(k) >= minn) mean(z[k]) else NA_real_ }, numeric(1))
d <- d %>% group_by(dyad) %>% mutate(Wa = win(year, z_mx, 4, 0), Wb = win(year, z_mx, 5, 1), Pa = win(year, z_mx, 9, 5), Pb = win(year, z_mx, 10, 6), cha = Wa - Pa, chb = Wb - Pb) %>% ungroup()
sink("results/longer_run_exposure_output.txt", split = TRUE)
cat("Windows (calendar years, >= 3 observed years): Wa = years t-4..t; Wb = t-5..t-1; change_a = Wa minus mean over t-9..t-5; change_b = Wb minus mean over t-10..t-6. Outcome: dispute begins in t+1.\n")
for (cg in c(0, 1)) { n5 <- d %>% filter(m == 0, contig == cg)
  cat(sprintf("\n%s pairs\n", ifelse(cg == 0, "NON-CONTIGUOUS", "CONTIGUOUS")))
  for (v in c("z_mx", "Wa", "Wb", "cha", "chb")) for (sp in c("", "+ z_vol")) {
    x <- n5 %>% filter(!is.na(Wa), if (v %in% c("cha")) !is.na(cha) else if (v == "chb") !is.na(chb) else if (v == "Wb") !is.na(Wb) else TRUE); x$l_dur <- log(x$dur_m); x$v <- x[[v]]
    m <- tryCatch(feglm(as.formula(paste("onset_m ~ v + z_mn + z_cap + l_dur + l_same + l_diff", sp, "| dyad")), x, family = "binomial", cluster = ~dyad, notes = FALSE), error = function(e) NULL)
    if (is.null(m)) { cat(v, "failed\n"); next }
    lab <- c(z_mx = "current-year share (same sample as Wa)", Wa = "mean, years t-4..t", Wb = "mean, years t-5..t-1", cha = "change: Wa minus mean t-9..t-5", chb = "change: Wb minus mean t-10..t-6")[v]
    cat(sprintf("%-40s %-9s %s | dyads with onset %d, obs %d\n", lab, ifelse(sp == "", "baseline", "+volume"), ci(coef(m)["v"], se(m)["v"]), length(unique(x$dyad[x$onset_m == 1])), nobs(m))) } }
sink()
})

# ---- Exposure and dispute onset/ending within the triad-embedded population (pairs with >=1 common tie partner).
local({
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), tri = (n_same + n_diff) >= 1)
sink("results/triad_population_output.txt", split = TRUE)
cat(sprintf("Dyad-years in panel: %d; with >=1 common tie partner (T): %d (%.1f%%)\n", nrow(d), sum(d$tri), 100 * mean(d$tri)))
cat(sprintf("Allied dyad-years: %d; of which in T: %d (%.1f%%); share of T that is allied: %.1f%%\n", sum(d$a == 1), sum(d$a == 1 & d$tri), 100 * mean(d$tri[d$a == 1]), 100 * mean(d$a[d$tri] == 1)))
cat("\n== Composition: embedded (T) x allied ==\n")
comp <- d %>% mutate(grp = case_when(tri & a == 1 ~ "Embedded, allied", tri & a == 0 ~ "Embedded, not allied", TRUE ~ "Not embedded")) %>% group_by(grp) %>%
  summarise(dyad_years = n(), no_dispute_years = sum(m == 0), onsets = sum(m == 0 & onset_m == 1), dispute_years = sum(m == 1), endings = sum(m == 1 & end_m == 1), .groups = "drop") %>%
  mutate(onset_rate_pct = round(100 * onsets / no_dispute_years, 2), ending_rate_pct = round(100 * endings / dispute_years, 1))
print(as.data.frame(comp)); cat(sprintf("Share of all onsets in T: %.1f%%; share of dispute-years in T: %.1f%%\n", 100 * mean(d$tri[d$m == 0 & d$onset_m == 1]), 100 * mean(d$tri[d$m == 1])))
fit <- function(x, y, label, extra = "", dur = "dur_m") {
  x$l_dur <- log(x[[dur]]); if (sum(x[[y]]) < 30) { cat(sprintf("%-58s too few events (%d)\n", label, sum(x[[y]]))); return(invisible()) }
  m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra)), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%-58s larger %6.3f (se %.3f, p=%.4f) | events %d of %d | dyads %d\n", label, ct["z_mx",1], ct["z_mx",2], ct["z_mx",4], sum(x[[y]]), nrow(x), n_distinct(x$dyad[x[[y]] == 1])))
}
d <- d %>% mutate(form_a = as.integer(a1 == 1), end_a = as.integer(a1 == 0))
for (pop in c("T (embedded)", "complement of T")) {
  s <- if (grepl("^T", pop)) d %>% filter(tri) else d %>% filter(!tri)
  cat("\n=====", pop, "=====\n")
  fit(s %>% filter(m == 0), "onset_m", "dispute onset: baseline")
  fit(s %>% filter(m == 0), "onset_m", "dispute onset: + volume + decade FE", "+ z_vol + decade")
  fit(s %>% filter(m == 0, contig == 1), "onset_m", "dispute onset: contiguous only")
  fit(s %>% filter(m == 1), "end_m", "dispute ending: baseline")
  fit(s %>% filter(m == 1), "end_m", "dispute ending: + volume + decade FE", "+ z_vol + decade")
  fit(s %>% filter(m == 1, contig == 1), "end_m", "dispute ending: contiguous only")
  fit(s %>% filter(a == 0, m == 0), "form_a", "alliance formation (no dispute): baseline", dur = "dur_a")
  fit(s %>% filter(a == 0, m == 0), "form_a", "alliance formation: + volume + decade FE", "+ z_vol + decade", dur = "dur_a")
  fit(s %>% filter(a == 0, m == 0, contig == 1), "form_a", "alliance formation: contiguous only", dur = "dur_a")
  fit(s %>% filter(a == 1, m == 0), "end_a", "alliance ending (no dispute): baseline", dur = "dur_a")
  fit(s %>% filter(a == 1, m == 0), "end_a", "alliance ending: + volume + decade FE", "+ z_vol + decade", dur = "dur_a")
  fit(s %>% filter(a == 1, m == 0, contig == 1), "end_a", "alliance ending: contiguous only", dur = "dur_a")
}
cat("\n===== pooled: larger share x T (embedded) =====\n")
for (spec in list(list("onset_m", d %>% filter(m == 0), "dur_m", "dispute onset"), list("end_m", d %>% filter(m == 1), "dur_m", "dispute ending"))) {
  x <- spec[[2]] %>% mutate(l_dur = log(.data[[spec[[3]]]]), T = as.integer(tri)); m <- glm(as.formula(paste(spec[[1]], "~ T * (z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff)")), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%-16s larger (not embedded) %6.3f (p=%.4f) | x embedded %6.3f (se %.3f, p=%.4f) | N=%d events=%d\n", spec[[4]], ct["z_mx",1], ct["z_mx",4], ct["T:z_mx",1], ct["T:z_mx",2], ct["T:z_mx",4], nobs(m), sum(x[[spec[[1]]]])))
}
cat("\n===== within T: larger share x alliance in force =====\n")
for (spec in list(list("onset_m", d %>% filter(tri, m == 0), "dispute onset"), list("end_m", d %>% filter(tri, m == 1), "dispute ending"))) {
  x <- spec[[2]] %>% mutate(l_dur = log(dur_m)); m <- glm(as.formula(paste(spec[[1]], "~ z_mx * a + z_mn + z_cap + contig + l_dur + l_same + l_diff")), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%-16s larger (no alliance) %6.3f (p=%.4f) | x alliance %6.3f (se %.3f, p=%.4f) | N=%d events=%d\n", spec[[3]], ct["z_mx",1], ct["z_mx",4], ct["z_mx:a",1], ct["z_mx:a",2], ct["z_mx:a",4], nobs(m), sum(x[[spec[[1]]]])))
}
sink()
})

# ---- Within-between (Mundlak) decomposition of the larger share for dispute onset, dispute ending, and alliance formation in the embedded population T.
local({
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), form_a = as.integer(a1 == 1), tri = (n_same + n_diff) >= 1)
fit <- function(x, y, dur, label, extra = "") {
  x <- x %>% group_by(dyad) %>% mutate(mx_b = mean(z_mx), mx_w = z_mx - mx_b, n_dyad_rows = n()) %>% ungroup(); x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ mx_b + mx_w + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra)), x, family = binomial)
  V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
  dif <- b["mx_b"] - b["mx_w"]; sed <- sqrt(V["mx_b", "mx_b"] + V["mx_w", "mx_w"] - 2 * V["mx_b", "mx_w"])
  cat(sprintf("%-46s between %s | within %s | difference %s | N=%d events=%d dyads=%d (dyads with >1 row: %d, with variation in the outcome: %d)\n", label, ci(b["mx_b"], sqrt(V["mx_b", "mx_b"])), ci(b["mx_w"], sqrt(V["mx_w", "mx_w"])), ci(dif, sed),
      nrow(x), sum(x[[y]]), n_distinct(x$dyad), n_distinct(x$dyad[x$n_dyad_rows > 1]), x %>% group_by(dyad) %>% filter(n_distinct(.data[[y]]) > 1) %>% pull(dyad) %>% n_distinct()))
}
sink("results/within_between_output.txt", split = TRUE)
cat("Larger share, per SD (global scale): between-dyad mean (mx_b) and within-dyad deviation (mx_w); embedded pairs.\n\n")
pops <- list(list("dispute onset", "onset_m", d %>% filter(tri, m == 0), "dur_m"), list("dispute ending", "end_m", d %>% filter(tri, m == 1), "dur_m"), list("alliance formation", "form_a", d %>% filter(tri, a == 0, m == 0), "dur_a"))
for (p in pops) { fit(p[[3]], p[[2]], p[[4]], paste(p[[1]], "baseline")); fit(p[[3]], p[[2]], p[[4]], paste(p[[1]], "+ volume + decade FE"), "+ z_vol + decade") }
sink()
})
