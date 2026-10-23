# 37: FULL PANEL: what distinguishes the distant, highly exposed pairs? Colonial ties, territorial claims, dispute issue, upper-tail events, coverage (longer-run within-pair exposure: see script 39).
#     Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(fixest); library(zoo) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
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
