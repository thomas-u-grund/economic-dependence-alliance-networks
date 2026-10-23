# 15: Hostility gradient of the exposure association (stakes vs interaction). Embedded pairs; risk set = no MID row (hostility >= 2) at t. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })

src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
mid <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb), fatlev = ifelse(fatlev < 0, NA, fatlev)) %>%
  group_by(year, lo, hi) %>% summarise(h = max(hihost), fat = suppressWarnings(max(fatlev, na.rm = TRUE)), .groups = "drop") %>% mutate(fat = ifelse(is.infinite(fat), NA, fat))
d <- d %>% left_join(mid, by = c("year", "lo", "hi")) %>% mutate(h = replace_na(h, 0), fat = replace_na(fat, 0), tri = (n_same + n_diff) >= 1)
nx <- d %>% transmute(year = year - 1, lo, hi, h1 = h, fat1 = fat); d <- d %>% left_join(nx, by = c("year", "lo", "hi")) %>% filter(!is.na(h1))
d <- d %>% arrange(dyad, year) %>% group_by(dyad) %>% mutate(run = cumsum(c(TRUE, diff(year) != 1 | diff(h > 0) != 0))) %>% group_by(dyad, run) %>% mutate(dur_q = row_number()) %>% ungroup()
R <- d %>% filter(tri, h == 0) %>% mutate(l_dur = log(dur_q), y_thr = as.integer(h1 == 2), y_dis = as.integer(h1 == 3), y_force = as.integer(h1 >= 4), y_any = as.integer(h1 >= 2), y_fat = as.integer(fat1 >= 1))
sink("results/hostility_gradient_output.txt", split = TRUE)
cat(sprintf("Risk set (embedded, no MID at t): %d dyad-years; onsets: threat %d, display %d, force/war %d, any %d, any fatalities %d\n", nrow(R), sum(R$y_thr), sum(R$y_dis), sum(R$y_force), sum(R$y_any), sum(R$y_fat)))
f <- function(y, extra = "") { m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra)), R, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = R$dyad)); c(est = ct["z_mx", 1], se = ct["z_mx", 2], lo = ct["z_mx", 1] - 1.96 * ct["z_mx", 2], hi = ct["z_mx", 1] + 1.96 * ct["z_mx", 2], p = ct["z_mx", 4], ev = sum(R[[y]])) }
for (sp in c("", "+ z_vol + decade")) { cat(sprintf("\n== larger share (per SD), spec: %s ==\n", ifelse(sp == "", "baseline", sp)))
  for (y in c("y_thr", "y_dis", "y_force", "y_any", "y_fat")) { r <- f(y, sp); cat(sprintf("  %-8s est %6.3f (se %.3f) CI [%6.3f, %6.3f] p=%.4f events %d\n", y, r["est"], r["se"], r["lo"], r["hi"], r["p"], r["ev"])) } }
cat("\n== stacked difference test: larger share x (force/war vs threat), and (display vs threat) ==\n")
for (pair in list(c("y_thr", "y_force"), c("y_thr", "y_dis"), c("y_dis", "y_force"))) {
  S <- bind_rows(R %>% mutate(y = .data[[pair[1]]], g = 0L), R %>% mutate(y = .data[[pair[2]]], g = 1L)); S <- S %>% filter(y == 1 | row_number() > 0)
  m <- glm(y ~ g * (z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff), S, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = S$dyad))
  cat(sprintf("  %s -> %s : difference in larger-share coefficient %6.3f (se %.3f, p=%.4f)\n", pair[1], pair[2], ct["g:z_mx", 1], ct["g:z_mx", 2], ct["g:z_mx", 4])) }
sink()
