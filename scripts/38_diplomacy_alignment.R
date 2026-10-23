# 38: FULL PANEL: diplomatic exchange (COW) and UN voting alignment (Bailey-Strezhnev-Voeten) as indicators of political engagement among the distant, exposed pairs; consistency across halves of the period.
#      Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
dx <- read_csv("raw/diplex/dipex.csv", show_col_types = FALSE) %>% transmute(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2), yr = year, DE = as.integer(DE)) %>% group_by(lo, hi, yr) %>% summarise(DE = max(DE), .groups = "drop")
# carry the latest observation (<= t) forward; before a dyad's first observation DE = 0
obs_years <- sort(unique(dx$yr)); d$yr <- sapply(d$year, function(y) { o <- obs_years[obs_years <= y]; if (length(o)) max(o) else NA })
d <- d %>% left_join(dx, by = c("lo", "hi", "yr")) %>% mutate(DE = replace_na(DE, 0L))
un <- read_csv("raw/diplex/unideal.csv", show_col_types = FALSE) %>% transmute(ccode, year = 1945 + session, ip = IdealPointAll)
d <- d %>% left_join(un %>% transmute(year, lo = ccode, ip_lo = ip), by = c("year", "lo")) %>% left_join(un %>% transmute(year, hi = ccode, ip_hi = ip), by = c("year", "hi")) %>%
  mutate(ipd = abs(ip_lo - ip_hi), onset_m = as.integer(m1 == 1), era = cut(year, c(1869, 1913, 1945, 1989, 2011), labels = c("1870-1913", "1914-45", "1946-89", "1990-2010")))
d$z_ipd <- as.numeric(scale(d$ipd))
nc <- d %>% filter(m == 0, contig == 0); p90 <- quantile(nc$mx, .9); p99 <- quantile(nc$mx, .99); nc$band <- cut(nc$mx, c(-Inf, p90, p99, Inf), labels = c("<p90", "p90-p99", ">=p99"))
fit <- function(x, extra, label, terms) { x$l_dur <- log(x$dur_m); m <- glm(as.formula(paste("onset_m ~ z_mx + z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad)
  for (tm in terms) cat(sprintf("%-60s %-12s %-10s %s | N=%d events=%d\n", label, ifelse(grepl("vol", extra), "+vol+decade", "baseline"), tm, ci(coef(m)[tm], sqrt(V[tm, tm])), nrow(x), sum(x$onset_m))) }
ig <- function(x, extra, label, terms = "zi") { x$zi <- x$z_mx * x$contig; fit(x %>% mutate(contig = contig), paste("+ contig + zi", extra), label, terms) }
sink("results/diplomacy_alignment_output.txt", split = TRUE)
cat("== D1. DIPLOMATIC EXCHANGE ==\n")
cat(sprintf("Non-contiguous dyad-years %d; with exchange %.1f%%\n", nrow(nc), 100 * mean(nc$DE)))
print(as.data.frame(nc %>% group_by(band) %>% summarise(dyad_years = n(), share_DE_pct = round(100 * mean(DE), 1), .groups = "drop")))
print(as.data.frame(nc %>% filter(band != "p90-p99") %>% group_by(era, band) %>% summarise(dyad_years = n(), share_DE_pct = round(100 * mean(DE), 1), .groups = "drop")))
cat(sprintf("Onset rate per 1,000 non-contiguous: DE=1 %.2f, DE=0 %.2f\n", 1000 * mean(nc$onset_m[nc$DE == 1]), 1000 * mean(nc$onset_m[nc$DE == 0])))
for (sp in c("", "+ z_vol + decade")) { fit(nc, sp, "non-contiguous, all (DE not controlled)", "z_mx"); fit(nc %>% filter(DE == 0), sp, "non-contiguous, no diplomatic exchange", "z_mx"); fit(nc %>% filter(DE == 1), sp, "non-contiguous, diplomatic exchange", "z_mx")
  x <- nc %>% mutate(zde = z_mx * DE); fit(x, paste("+ DE + zde", sp), "non-contiguous, interaction share x DE", c("DE", "zde"))
  a <- d %>% filter(m == 0); a$l_dur <- log(a$dur_m); a$zi <- a$z_mx * a$contig; a$zde <- a$z_mx * a$DE
  fit(a, paste("+ contig + zi", sp), "all dyads, geography interaction (no DE)", "zi"); fit(a, paste("+ contig + zi + DE + zde", sp), "all dyads, geography interaction with DE and DE x share", "zi") }
cat("\n== D2. UN VOTING ALIGNMENT (1946-2010) ==\n")
n2 <- nc %>% filter(!is.na(z_ipd)); a2 <- d %>% filter(m == 0, !is.na(z_ipd))
cat(sprintf("Non-contiguous 1946-2010 with ideal points: %d dyad-years, %d onsets (of %d)\n", nrow(n2), sum(n2$onset_m), sum(nc$onset_m[nc$year >= 1946])))
print(as.data.frame(n2 %>% group_by(band) %>% summarise(dyad_years = n(), mean_distance = round(mean(ipd), 3), median_distance = round(median(ipd), 3), .groups = "drop")))
for (sp in c("", "+ z_vol + decade")) { fit(n2, sp, "non-contiguous 1946-2010, no alignment control", "z_mx"); fit(n2, paste("+ z_ipd", sp), "non-contiguous 1946-2010, alignment control", c("z_mx", "z_ipd"))
  x <- n2 %>% mutate(zip = z_mx * z_ipd); fit(x, paste("+ z_ipd + zip", sp), "non-contiguous 1946-2010, interaction share x distance", c("z_mx", "z_ipd", "zip"))
  a2$zi <- a2$z_mx * a2$contig; a2$zip <- a2$z_mx * a2$z_ipd; fit(a2, paste("+ contig + zi", sp), "all dyads 1946-2010, geography interaction (no alignment)", "zi"); fit(a2, paste("+ contig + zi + z_ipd + zip", sp), "all dyads 1946-2010, geography interaction with alignment", "zi") }
cat("\n== D3. HALVES OF THE PERIOD ==\n")
for (h in list(c(1870, 1945), c(1946, 2010))) { a <- d %>% filter(m == 0, year >= h[1], year <= h[2]); a$zi <- a$z_mx * a$contig
  for (sp in c("", "+ z_vol + decade")) { fit(a, paste("+ contig + zi", sp), sprintf("all dyads %d-%d, geography interaction", h[1], h[2]), "zi"); fit(a %>% filter(contig == 0), sp, sprintf("non-contiguous %d-%d, larger share", h[1], h[2]), "z_mx") } }
sink()
