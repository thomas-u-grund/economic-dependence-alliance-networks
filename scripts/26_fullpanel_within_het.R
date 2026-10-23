# 26: FULL PANEL: within-between decomposition, larger share x log capital distance, escalation given onset. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(cshapes) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
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
