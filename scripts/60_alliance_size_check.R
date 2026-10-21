# 60: Does the allied-dispute pattern (and the exposure association) depend on large regional alliances?
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE)
size <- mem %>% group_by(version4id) %>% summarise(n_mem = n_distinct(ccode), .groups = "drop")
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>% left_join(size, by = "version4id") %>%
  mutate(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>% group_by(year, lo, hi) %>% summarise(min_n = min(n_mem), .groups = "drop") %>%
  mutate(asize = cut(min_n, c(0, 3, 7, Inf), labels = c("small (<=3)", "medium (4-7)", "large (>=8)")))
d <- d0 %>% left_join(al, by = c("year", "lo", "hi"))
cat_ <- function(x) { x$asize <- as.character(x$asize); x }
sink("results/alliance_size_check_output.txt", split = TRUE)
cat("== Counts by size of the smallest alliance covering the dyad (full signed panel) ==\n")
a <- d %>% filter(!is.na(asize)) %>% group_by(asize) %>% summarise(allied_dyad_years = n(), dispute_years = sum(sign == -1), dispute_ends_in_alliance = sum(sign == -1 & flip_i == 1, na.rm = TRUE), alliance_years_no_dispute = sum(sign == 1), alliance_sign_changes = sum(sign == 1 & flip_i == 1, na.rm = TRUE), .groups = "drop")
print(as.data.frame(a))
prep <- function(x) x %>% filter(!is.na(dep), !is.na(mut), !is.na(cap)) %>% mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
fit1 <- function(x, label) { if (sum(x$flip_i) < 30) { cat(sprintf("%-52s too few events (%d)\n", label, sum(x$flip_i))); return(invisible()) }
  m <- glm(flip_i ~ z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal, x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%-52s larger share %6.3f (se %.3f, p=%.4f) | N=%d events=%d dyads=%d\n", label, ct["z_dep",1], ct["z_dep",2], ct["z_dep",4], nobs(m), sum(x$flip_i), n_distinct(x$dyad))) }
cat("\n== Larger-share coefficient by alliance size (trade sample) ==\n")
for (g in c("small (<=3)", "medium (4-7)", "large (>=8)")) {
  fit1(prep(d %>% filter(sign == -1, asize == g)), paste("dispute ends inside alliance |", g))
  fit1(prep(d %>% filter(sign == 1, asize == g)), paste("alliance turns hostile |", g)) }
cat("\n== Pooled interaction larger share x large alliance ==\n")
for (sg in c(-1, 1)) { x <- prep(d %>% filter(sign == sg, !is.na(asize))) %>% mutate(large = as.integer(asize == "large (>=8)"))
  m <- glm(flip_i ~ z_dep * large + z_mut + z_cap + l_dur + l_unbal + l_bal, x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%s: larger share (not large) %.3f (p=%.4f) | x large %.3f (se %.3f, p=%.4f) | N=%d events=%d\n", ifelse(sg == -1, "disputes inside alliances", "alliances"), ct["z_dep",1], ct["z_dep",4], ct["z_dep:large",1], ct["z_dep:large",2], ct["z_dep:large",4], nobs(m), sum(x$flip_i))) }
sink()
