# 21: Escalation conditional on onset: among dispute onsets in the embedded population, does the larger share predict use of force/war (vs display) and fatalities? Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
mid <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb), fatlev = ifelse(fatlev < 0, NA, fatlev)) %>%
  group_by(year, lo, hi) %>% summarise(h1 = max(hihost), fat1 = suppressWarnings(max(fatlev, na.rm = TRUE)), .groups = "drop") %>% mutate(fat1 = ifelse(is.infinite(fat1), NA, fat1), year = year - 1)  # year t = onset year - 1
d <- d %>% mutate(tri = (n_same + n_diff) >= 1) %>% filter(tri, m == 0, m1 == 1) %>% left_join(mid, by = c("year", "lo", "hi")) %>% filter(!is.na(h1)) %>% mutate(l_dur = log(dur_m), y_force = as.integer(h1 >= 4), y_fat = as.integer(!is.na(fat1) & fat1 >= 1))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
sink("results/escalation_given_onset_output.txt", split = TRUE)
cat(sprintf("Onsets in embedded pairs with hostility information: %d; use of force or war: %d (%.1f%%); any fatalities: %d (%.1f%%)\n", nrow(d), sum(d$y_force), 100 * mean(d$y_force), sum(d$y_fat), 100 * mean(d$y_fat)))
for (sp in c("", "+ z_vol + decade")) { cat(sprintf("\n== larger share (per SD), spec: %s ==\n", ifelse(sp == "", "baseline", sp)))
  for (y in c("y_force", "y_fat")) { m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff", sp)), d, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = d$dyad))
    cat(sprintf("  %-8s %s p=%.4f events %d of %d\n", y, ci(ct["z_mx", 1], ct["z_mx", 2]), ct["z_mx", 4], sum(d[[y]]), nrow(d))) } }
sink()
