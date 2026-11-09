# 4_exposed_pairs: Section 5.3: the highly exposed non-contiguous pairs (Appendix C, Tables 8 to 12).
# Run from the repository root: Rscript scripts/4_exposed_pairs.R (or run everything with Rscript run_all.R).
# Each block below is self-contained and writes the results/ file named in its sink() call.
source("scripts/setup.R")

# ---- FULL PANEL: who populates the upper tail of exposure among non-contiguous pairs? Description and leave-out checks.
local({
nmc <- read_csv("raw/nmc/NMCv7/abridged/NMC-70-abridged.csv", show_col_types = FALSE) %>% select(ccode, year, cinc) %>% group_by(year) %>% mutate(rk = rank(-cinc, ties.method = "min")) %>% ungroup()
nm <- read_csv("raw/states/States2024/statelist2024.csv", show_col_types = FALSE) %>% distinct(ccode, .keep_all = TRUE) %>% select(ccode, statenme)
d <- d %>% mutate(onset_m = as.integer(m1 == 1)) %>% left_join(nmc %>% transmute(year, lo = ccode, rk_lo = rk), by = c("year", "lo")) %>% left_join(nmc %>% transmute(year, hi = ccode, rk_hi = rk), by = c("year", "hi")) %>%
  mutate(top5 = as.integer(coalesce(rk_lo, 99) <= 5 | coalesce(rk_hi, 99) <= 5), era = cut(year, c(1869, 1913, 1945, 1989, 2011), labels = c("1870-1913", "1914-45", "1946-89", "1990-2010")))
nc <- d %>% filter(m == 0, contig == 0); p99 <- quantile(nc$mx, .99); tl <- nc %>% filter(mx >= p99)
cat_name <- function(x) nm$statenme[match(x, nm$ccode)]
sink("results/upper_tail_output.txt", split = TRUE)
cat(sprintf("Non-contiguous onset risk set: %d dyad-years; p99 of the larger share = %.4f (%.1f%%)\n", nrow(nc), p99, 100 * p99))
cat(sprintf("Upper tail (>= p99): %d dyad-years, %d dyads, %d onsets (rate %.2f per 1,000) vs rest %.2f per 1,000\n", nrow(tl), n_distinct(tl$dyad), sum(tl$onset_m), 1000 * mean(tl$onset_m), 1000 * mean(nc$onset_m[nc$mx < p99])))
cat(sprintf("Share allied (defence pact): %.1f%%; with a top-5 CINC state: %.1f%% (all non-contiguous risk set: %.1f%% allied, %.1f%% top-5)\n", 100 * mean(tl$a), 100 * mean(tl$top5), 100 * mean(nc$a), 100 * mean(nc$top5)))
cat("\nBy era: tail dyad-years, share of tail, share of all non-contiguous dyad-years, tail onsets\n")
print(as.data.frame(tl %>% count(era, name = "tail_n") %>% left_join(nc %>% count(era, name = "all_n"), by = "era") %>% left_join(tl %>% group_by(era) %>% summarise(onsets = sum(onset_m)), by = "era") %>% mutate(tail_share = round(tail_n / sum(tail_n), 3), all_share = round(all_n / sum(all_n), 3))))
st <- bind_rows(tl %>% transmute(ccode = lo, dyad), tl %>% transmute(ccode = hi, dyad)) %>% count(ccode, name = "n") %>% arrange(desc(n)) %>% mutate(name = cat_name(ccode), share = round(n / nrow(tl), 3))
cat("\nMost frequent states in the tail (dyad-years in which the state is a member)\n"); print(as.data.frame(head(st, 12)))
dy <- tl %>% group_by(lo, hi) %>% summarise(n = n(), first = min(year), last = max(year), onsets = sum(onset_m), allied = mean(a), .groups = "drop") %>% arrange(desc(n)) %>% mutate(pair = paste(cat_name(lo), "-", cat_name(hi)), share = round(n / nrow(tl), 3), allied = round(allied, 2))
cat("\nMost frequent dyads in the tail\n"); print(as.data.frame(head(dy %>% select(pair, n, first, last, onsets, allied, share), 15)))
cat(sprintf("\nTop 10 dyads account for %.1f%% of tail dyad-years; top 5 states appear in %.1f%%\n", 100 * sum(head(dy$n, 10)) / nrow(tl), 100 * mean(tl$lo %in% head(st$ccode, 5) | tl$hi %in% head(st$ccode, 5))))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
est <- function(x, label) for (sp in c("", "+ z_vol + decade")) { x$l_dur <- log(x$dur_m)
  m <- glm(as.formula(paste("onset_m ~ z_mx + z_mn + z_cap + l_dur + l_same + l_diff", sp)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad)
  cat(sprintf("%-44s %-16s %s | N=%d events=%d\n", label, ifelse(sp == "", "baseline", "+ vol + decade"), ci(coef(m)["z_mx"], sqrt(V["z_mx", "z_mx"])), nrow(x), sum(x$onset_m))) }
cat("\nNon-contiguous larger-share association, leaving out parts of the sample (dyad-clustered)\n")
est(nc, "all non-contiguous")
for (e in levels(nc$era)) est(nc %>% filter(era != e), paste("leave out era", e))
for (s in head(st$ccode, 5)) est(nc %>% filter(lo != s, hi != s), paste("leave out", cat_name(s)))
est(nc %>% filter(!(paste(lo, hi) %in% paste(head(dy$lo, 10), head(dy$hi, 10)))), "leave out top 10 tail dyads")
est(nc %>% filter(top5 == 0), "leave out any pair with a top-5 state")
sink()
})

# ---- FULL PANEL: what distinguishes the distant, highly exposed pairs? Colonial ties, territorial claims, dispute issue, upper-tail events, coverage (longer-run within-pair exposure: see scripts/5_within_between.R).
local({
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
# ---- new data ----
col <- read_csv("raw/icow/colhist/ICOW Colonial History 1.1/coldata110.csv", show_col_types = FALSE)
colp <- bind_rows(lapply(c("ColRuler", "IndFrom", "SecFrom"), function(v) col %>% transmute(s = State, r = .data[[v]], y0 = ifelse(IndDate > 0, IndDate %/% 100, NA)))) %>%
  filter(r > 0, !is.na(y0)) %>% transmute(lo = pmin(s, r), hi = pmax(s, r), y0) %>% group_by(lo, hi) %>% summarise(y0 = min(y0), .groups = "drop")
terr <- read_csv("raw/icow/terr/ICOWdyadyr.csv", col_select = c(chal, tgt, year), show_col_types = FALSE) %>%
  transmute(year, lo = pmin(chal, tgt), hi = pmax(chal, tgt), claim = 1L) %>% distinct()
iss <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% filter(hihost >= 3) %>%
  transmute(year, lo = pmin(statea, stateb), hi = pmax(statea, stateb), terr_i = as.integer(revtypea == 1 | revtypeb == 1)) %>% group_by(year, lo, hi) %>% summarise(terr_i = max(terr_i), .groups = "drop")
d <- d %>% left_join(colp, by = c("lo", "hi")) %>% mutate(colonial = as.integer(!is.na(y0) & year >= y0)) %>% select(-y0) %>%
  left_join(terr, by = c("year", "lo", "hi")) %>% mutate(claim = replace_na(claim, 0L)) %>%
  left_join(iss %>% transmute(year = year - 1, lo, hi, terr_n = terr_i), by = c("year", "lo", "hi")) %>%
  mutate(onset_m = as.integer(m1 == 1), form_a = as.integer(a1 == 1), onset_terr = as.integer(m1 == 1 & coalesce(terr_n, 0L) == 1), onset_non = as.integer(m1 == 1 & coalesce(terr_n, 0L) == 0),
         era = cut(year, c(1869, 1913, 1945, 1989, 2011), labels = c("1870-1913", "1914-45", "1946-89", "1990-2010")))
nc <- d %>% filter(m == 0, contig == 0); p90 <- quantile(nc$mx, .9); p99 <- quantile(nc$mx, .99)
est1 <- function(x, y, dur, extra, label, term = "z_mx") { x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad)
  cat(sprintf("%-58s %-15s %s | N=%d events=%d\n", label, ifelse(grepl("vol", extra), "+vol+decade", "baseline"), ci(coef(m)[term], sqrt(V[term, term])), nrow(x), sum(x[[y]]))); invisible(m) }
sink("results/what_distinguishes_output.txt", split = TRUE)
cat("== 1. COLONIAL TIES ==\n")
tl <- nc %>% filter(mx >= p99)
cat(sprintf("Non-contiguous risk set %d dyad-years; colonial-tie dyad-years %d (%.1f%%), dyads %d. Top-1%% tail: %d dyad-years, colonial %.1f%%; tail onsets %d, of which in colonial-tie pairs %d. Elsewhere in risk set: colonial %.1f%%\n",
  nrow(nc), sum(nc$colonial), 100 * mean(nc$colonial), n_distinct(nc$dyad[nc$colonial == 1]), nrow(tl), 100 * mean(tl$colonial), sum(tl$onset_m), sum(tl$onset_m[tl$colonial == 1]), 100 * mean(nc$colonial[nc$mx < p99])))
cat(sprintf("Onset rate per 1,000 non-contiguous: colonial %.2f vs other %.2f\n", 1000 * mean(nc$onset_m[nc$colonial == 1]), 1000 * mean(nc$onset_m[nc$colonial == 0])))
for (sp in c("", "+ z_vol + decade")) { est1(nc, "onset_m", "dur_m", sp, "onset, non-contiguous, all")
  est1(nc %>% filter(colonial == 0), "onset_m", "dur_m", sp, "onset, non-contiguous, colonial-tie pairs dropped")
  est1(nc %>% filter(colonial == 1), "onset_m", "dur_m", sp, "onset, non-contiguous, colonial-tie pairs only")
  est1(nc %>% mutate(zc = z_mx * colonial), "onset_m", "dur_m", paste("+ colonial + zc", sp), "onset, non-contiguous, interaction", term = "zc") }
af <- d %>% filter(a == 0, m == 0, contig == 0)
for (sp in c("", "+ z_vol + decade")) { est1(af, "form_a", "dur_a", sp, "alliance formation, non-contiguous, all")
  est1(af %>% filter(colonial == 0), "form_a", "dur_a", sp, "alliance formation, non-contiguous, colonial-tie pairs dropped") }
cat("\n== 2. TERRITORIAL CLAIMS (years <= 2001) ==\n")
n2 <- nc %>% filter(year <= 2001); a2 <- d %>% filter(m == 0, year <= 2001)
cat(sprintf("Non-contiguous <=2001: %d dyad-years, with an active claim %d (%.2f%%); claim onset rate %.2f vs %.2f per 1,000. Top-1%% tail with a claim: %.1f%%\n", nrow(n2), sum(n2$claim), 100 * mean(n2$claim), 1000 * mean(n2$onset_m[n2$claim == 1]), 1000 * mean(n2$onset_m[n2$claim == 0]), 100 * mean((n2 %>% filter(mx >= p99))$claim)))
for (sp in c("", "+ z_vol + decade")) { est1(n2, "onset_m", "dur_m", sp, "onset, non-contiguous <=2001, all"); est1(n2 %>% filter(claim == 0), "onset_m", "dur_m", sp, "onset, non-contiguous <=2001, claim pairs dropped") }
ig <- function(x, extra, label, y = "onset_m") { x$l_dur <- log(x$dur_m); x$zi <- x$z_mx * x$contig; x$zcl <- x$z_mx * x$claim
  m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + zi + l_dur + l_same + l_diff", extra)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad)
  cat(sprintf("%-58s geography interaction %s | N=%d events=%d\n", label, ci(coef(m)["zi"], sqrt(V["zi", "zi"])), nrow(x), sum(x[[y]]))) }
for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+vol+decade")
  ig(a2, sp, paste("all dyads <=2001, no claim control,", s)); ig(a2, paste("+ claim + I(z_mx*claim)", sp), paste("all dyads <=2001, claim and claim x share,", s)) }
cat("\n== 3. DISPUTE ISSUE (territorial revisionism vs other) ==\n")
cat(sprintf("Onsets in the full risk set (no dispute at t): %d; territorial %d, other %d. Non-contiguous onsets: territorial %d, other %d. Contiguous: territorial %d, other %d\n", sum(d$onset_m[d$m == 0]), sum(d$onset_terr[d$m == 0]), sum(d$onset_non[d$m == 0]), sum(nc$onset_terr), sum(nc$onset_non), sum(d$onset_terr[d$m == 0 & d$contig == 1]), sum(d$onset_non[d$m == 0 & d$contig == 1])))
a3 <- d %>% filter(m == 0)
for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+vol+decade")
  ig(a3, sp, paste("territorial-issue onsets,", s), "onset_terr"); ig(a3, sp, paste("other-issue onsets,", s), "onset_non")
  est1(nc, "onset_terr", "dur_m", sp, "territorial-issue onset, non-contiguous"); est1(nc, "onset_non", "dur_m", sp, "other-issue onset, non-contiguous")
  est1(a3 %>% filter(contig == 1), "onset_terr", "dur_m", sp, "territorial-issue onset, contiguous"); est1(a3 %>% filter(contig == 1), "onset_non", "dur_m", sp, "other-issue onset, contiguous") }
cat("\n== 4. UPPER TAIL BY PERIOD: non-contiguous dyad-years and onsets by exposure band ==\n")
nc$band <- cut(nc$mx, c(-Inf, p90, p99, Inf), labels = c("<p90", "p90-p99", ">=p99"))
print(as.data.frame(nc %>% group_by(era, band) %>% summarise(dyad_years = n(), dyads = n_distinct(dyad), onsets = sum(onset_m), per1000 = round(1000 * mean(onset_m), 2), .groups = "drop")))
nm <- read_csv("raw/states/States2024/statelist2024.csv", show_col_types = FALSE) %>% distinct(ccode, .keep_all = TRUE) %>% select(ccode, statenme)
ev <- tl %>% filter(onset_m == 1) %>% transmute(pair = paste(nm$statenme[match(lo, nm$ccode)], "-", nm$statenme[match(hi, nm$ccode)]), year_t = year, colonial, claim, allied = a) %>% arrange(year_t)
cat(sprintf("\nTop-1%% non-contiguous onsets (%d): by dyad\n", nrow(ev))); print(as.data.frame(ev %>% count(pair, name = "onsets") %>% arrange(desc(onsets)) %>% head(20)))
cat("\n== 5. COVERAGE: dyadic MID onsets (hostility >= 3) inside the trade panel ==\n")
allm <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% filter(hihost >= 3, year >= 1870, year <= 2010) %>% transmute(year, lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>% distinct()
inp <- d %>% distinct(year, lo, hi, contig) %>% mutate(inpanel = 1L)
cv <- allm %>% left_join(inp, by = c("year", "lo", "hi")) %>% mutate(inpanel = replace_na(inpanel, 0L), era = cut(year, c(1869, 1913, 1945, 1989, 2011), labels = c("1870-1913", "1914-45", "1946-89", "1990-2010")))
cat(sprintf("Dispute dyad-years (hostility >= 3) 1870-2010: %d; inside the trade panel: %d (%.1f%%)\n", nrow(cv), sum(cv$inpanel), 100 * mean(cv$inpanel)))
print(as.data.frame(cv %>% group_by(era) %>% summarise(dispute_dyad_years = n(), in_panel = sum(inpanel), pct = round(100 * mean(inpanel), 1))))
sink()
})

# ---- FULL PANEL: diplomatic exchange (COW) and UN voting alignment (Bailey-Strezhnev-Voeten) as indicators of political engagement among the distant, exposed pairs; consistency across halves of the period.
local({
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
})
