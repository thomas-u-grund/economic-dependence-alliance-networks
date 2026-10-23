# 04: Exposure and dispute onset/ending within the triad-embedded population (pairs with >=1 common tie partner). Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })

src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]
eval(parse(text = src[1:(k - 1)]))   # builds d (full panel with l_same, l_diff, z_mx, z_mn, z_cap, z_vol, contig, dur_m, dur_a, a, m, a1, m1, dyad, decade)
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
