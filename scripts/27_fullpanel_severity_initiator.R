# 27: FULL PANEL: severity of onsets by hostility level (risk set: no MID at t) with stacked difference tests, and initiation. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
mid <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb), fatlev = ifelse(fatlev < 0, NA, fatlev)) %>%
  group_by(year, lo, hi) %>% summarise(h = max(hihost), fat = suppressWarnings(max(fatlev, na.rm = TRUE)), .groups = "drop") %>% mutate(fat = ifelse(is.infinite(fat), NA, fat))
d <- d %>% left_join(mid, by = c("year", "lo", "hi")) %>% mutate(h = replace_na(h, 0), fat = replace_na(fat, 0))
nx <- d %>% transmute(year = year - 1, lo, hi, h1 = h, fat1 = fat); d <- d %>% left_join(nx, by = c("year", "lo", "hi")) %>% filter(!is.na(h1))
d <- d %>% arrange(dyad, year) %>% group_by(dyad) %>% mutate(run = cumsum(c(TRUE, diff(year) != 1 | diff(h > 0) != 0))) %>% group_by(dyad, run) %>% mutate(dur_q = row_number()) %>% ungroup()
R <- d %>% filter(h == 0) %>% mutate(l_dur = log(dur_q), y_thr = as.integer(h1 == 2), y_dis = as.integer(h1 == 3), y_force = as.integer(h1 >= 4), y_any = as.integer(h1 >= 2), y_fat = as.integer(fat1 >= 1))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
sink("results/fullpanel_severity_initiator_output.txt", split = TRUE)
cat(sprintf("Risk set (full panel, no MID at t): %d dyad-years; onsets: threat %d, display %d, force/war %d, any %d, any fatalities %d\n", nrow(R), sum(R$y_thr), sum(R$y_dis), sum(R$y_force), sum(R$y_any), sum(R$y_fat)))
f <- function(y, extra = "") { m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra)), R, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = R$dyad)); list(ct["z_mx", 1], ct["z_mx", 2], ct["z_mx", 4], sum(R[[y]])) }
for (sp in c("", "+ z_vol + decade")) { cat(sprintf("\n== larger share (per SD), spec: %s ==\n", ifelse(sp == "", "baseline", sp)))
  for (y in c("y_thr", "y_dis", "y_force", "y_fat", "y_any")) { r <- f(y, sp); cat(sprintf("  %-8s %s p=%.4f events %d\n", y, ci(r[[1]], r[[2]]), r[[3]], r[[4]])) } }
cat("\n== stacked difference tests ==\n")
for (pair in list(c("y_thr", "y_force"), c("y_thr", "y_dis"), c("y_dis", "y_force"))) {
  S <- bind_rows(R %>% mutate(y = .data[[pair[1]]], g = 0L), R %>% mutate(y = .data[[pair[2]]], g = 1L)); m <- glm(y ~ g * (z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff), S, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = S$dyad))
  cat(sprintf("  %s -> %s : difference in larger-share coefficient %6.3f (se %.3f, p=%.4f)\n", pair[1], pair[2], ct["g:z_mx", 1], ct["g:z_mx", 2], ct["g:z_mx", 4])) }
# ---- initiation
nmc <- read_csv("raw/nmc/NMCv7/abridged/NMC-70-abridged.csv", show_col_types = FALSE) %>% select(ccode, year, cinc)
mid2 <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% filter(hihost >= 3) %>% mutate(init = ifelse(dyad_rolea == 1, statea, ifelse(dyad_roleb == 1, stateb, NA)), lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>%
  filter(!is.na(init)) %>% group_by(year, lo, hi) %>% summarise(n_init = n_distinct(init), init = first(init), .groups = "drop") %>% filter(n_init == 1)
ev <- d %>% filter(m == 0, m1 == 1) %>% mutate(year1 = year + 1) %>% left_join(mid2, by = c("year1" = "year", "lo", "hi")) %>% left_join(nmc %>% rename(lo = ccode, cinc_lo = cinc), by = c("year", "lo")) %>% left_join(nmc %>% rename(hi = ccode, cinc_hi = cinc), by = c("year", "hi"))
e <- ev %>% filter(!is.na(init), !is.na(cinc_lo), !is.na(cinc_hi), cinc_lo > 0, cinc_hi > 0) %>% mutate(lo_init = as.integer(init == lo), dshare = (shl - shh) / 0.0478, dcinc = log(cinc_lo) - log(cinc_hi), more_exp_init = as.integer((shl > shh & init == lo) | (shh > shl & init == hi)), tied = shl == shh)
cat(sprintf("\n== initiation: onsets %d, with a clear initiator and capability data %d (dyads %d) ==\n", nrow(ev), nrow(e), n_distinct(e$dyad)))
x <- e %>% filter(!tied); b <- binom.test(sum(x$more_exp_init), nrow(x)); cat(sprintf("The MORE exposed state initiates in %d of %d (%.3f; 95%% CI %.3f-%.3f); the more exposed state is the weaker one in %.3f; the stronger state initiates in %.3f\n", sum(x$more_exp_init), nrow(x), b$estimate, b$conf.int[1], b$conf.int[2], mean(((x$shl > x$shh) == (x$cinc_lo < x$cinc_hi))), mean((e$cinc_lo > e$cinc_hi) == (e$lo_init == 1))))
for (nm in c("exposure only", "+ capability")) { m <- if (nm == "exposure only") glm(lo_init ~ 0 + dshare, e, family = binomial) else glm(lo_init ~ 0 + dshare + dcinc, e, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = e$dyad)); cat(sprintf("  %-14s", nm)); print(round(ct, 4)) }
sink()
