# 24: Full-panel versions of the main results, and lagged definition of embedding. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
cl <- get_cnt(d$year - 1, d$lo, d$hi); d$n_same_l <- cl$s; d$n_diff_l <- cl$d; d$l_same_l <- log1p(d$n_same_l); d$l_diff_l <- log1p(d$n_diff_l)
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), form_a = as.integer(a1 == 1), end_a = as.integer(a1 == 0), tri = (n_same + n_diff) >= 1, tri_l = (n_same_l + n_diff_l) >= 1)
mid <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb), fatlev = ifelse(fatlev < 0, NA, fatlev)) %>%
  group_by(year, lo, hi) %>% summarise(h = max(hihost), fat = suppressWarnings(max(fatlev, na.rm = TRUE)), .groups = "drop") %>% mutate(fat = ifelse(is.infinite(fat), NA, fat))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
est1 <- function(x, y, dur, label, sp = "", same = "l_same", diff = "l_diff", inter = NULL) {
  x$l_dur <- log(x[[dur]]); rhs <- if (is.null(inter)) "z_mx" else paste("z_mx *", inter)
  m <- glm(as.formula(paste(y, "~", rhs, "+ z_mn + z_cap + contig + l_dur +", same, "+", diff, sp)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  out <- sprintf("%-46s larger share %s", label, ci(b["z_mx"], sqrt(V["z_mx", "z_mx"])))
  if (!is.null(inter)) { it <- paste0("z_mx:", inter); out <- paste0(out, sprintf(" | interaction %s", ci(b[it], sqrt(V[it, it])))) }
  cat(sprintf("%s | N=%d events=%d dyads=%d\n", out, nrow(x), sum(x[[y]]), n_distinct(x$dyad[x[[y]] == 1])))
}
pops <- function(dd) list(list("dispute onset", "onset_m", "dur_m", dd %>% filter(m == 0)), list("dispute ending", "end_m", "dur_m", dd %>% filter(m == 1)), list("alliance formation", "form_a", "dur_a", dd %>% filter(a == 0, m == 0)), list("alliance ending", "end_a", "dur_a", dd %>% filter(a == 1, m == 0)))
sink("results/population_selection_output.txt", split = TRUE)
cat("######## (a) FULL PANEL, contemporaneous partner counts as controls ########\n")
for (p in pops(d)) for (sp in c("", "+ z_vol + decade")) est1(p[[4]], p[[2]], p[[3]], paste(p[[1]], ifelse(sp == "", "baseline", "+ vol + decade FE")), sp)
cat("\n-- full panel, larger share x contiguity (onset, formation)\n")
for (p in pops(d)[c(1, 3)]) for (sp in c("", "+ z_vol + decade")) est1(p[[4]], p[[2]], p[[3]], paste(p[[1]], ifelse(sp == "", "baseline", "+ vol + decade FE")), sp, inter = "contig")
cat("\n-- full panel: severity of onsets (risk set: no MID at t)\n")
dm <- d %>% left_join(mid, by = c("year", "lo", "hi")) %>% mutate(h = replace_na(h, 0)) %>% filter(h == 0) %>% left_join(mid %>% transmute(year = year - 1, lo, hi, h1 = h, fat1 = fat), by = c("year", "lo", "hi")) %>%
  mutate(h1 = replace_na(h1, 0), fat1 = replace_na(fat1, 0), y_dis = as.integer(h1 == 3), y_force = as.integer(h1 >= 4), y_fat = as.integer(fat1 >= 1), y_any = as.integer(h1 >= 2))
dm <- dm %>% arrange(dyad, year) %>% group_by(dyad) %>% mutate(run = cumsum(c(TRUE, diff(year) != 1))) %>% group_by(dyad, run) %>% mutate(dur_q = row_number()) %>% ungroup()
for (y in c("y_dis", "y_force", "y_fat", "y_any")) for (sp in c("", "+ z_vol + decade")) est1(dm, y, "dur_q", paste("severity", y, ifelse(sp == "", "baseline", "+ vol + decade FE")), sp)
cat("\n######## (b) LAGGED EMBEDDING: embedded pair and partner counts from the network of year t-1 ########\n")
for (p in pops(d %>% filter(tri_l))) for (sp in c("", "+ z_vol + decade")) est1(p[[4]], p[[2]], p[[3]], paste("lagged-embedded", p[[1]], ifelse(sp == "", "baseline", "+ vol + decade FE")), sp, same = "l_same_l", diff = "l_diff_l")
cat(sprintf("\nShare of contemporaneously embedded dyad-years that are also embedded at t-1: %.3f; of t-1 embedded that are embedded at t: %.3f\n", mean(d$tri_l[d$tri]), mean(d$tri[d$tri_l])))
sink()
