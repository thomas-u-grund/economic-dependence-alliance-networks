# 40: Check 3a redone: drop dispute-years during which an alliance covering the same dyad is in force.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>%
  transmute(year, lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2), has_all = TRUE) %>% distinct()
d <- d0 %>% left_join(al, by = c("year", "lo", "hi")) %>% mutate(overlap = sign == -1 & !is.na(has_all))
cat("dispute-years:", sum(d$sign == -1), " with alliance in force:", sum(d$overlap), sprintf("(%.1f%%)\n", 100 * mean(d$overlap[d$sign == -1])))
zs <- function(x) as.numeric(scale(x))
sink("results/hazard_overlap_exclusion_output.txt", split = TRUE)
cat("dispute-years:", sum(d$sign == -1), " with alliance in force:", sum(d$overlap), "\n")
for (lab in c("all disputes", "disputes WITHOUT concurrent alliance", "disputes WITH concurrent alliance")) {
  s <- d %>% filter(sign == -1, !is.na(dep), !is.na(mut), !is.na(cap))
  s <- if (lab == "all disputes") s else if (grepl("WITHOUT", lab)) s %>% filter(!overlap) else s %>% filter(overlap)
  s <- s %>% mutate(z_dep = zs(dep), z_mut = zs(mut), z_cap = zs(cap), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
  for (y in c("flip_i", "gone_i")) { if (sum(s[[y]]) < 15) { cat(lab, y, ": too few events\n"); next }
    m <- glm(as.formula(paste(y, "~ z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal")), s, family = binomial)
    ct <- coeftest(m, vcov = vcovCL(m, cluster = s$dyad))
    cat(sprintf("%-40s %-7s dep=%6.3f (se %.3f, p=%.4f) | N=%d ev=%d dyads=%d\n", lab, y, ct["z_dep",1], ct["z_dep",2], ct["z_dep",4], nobs(m), sum(s[[y]]), n_distinct(s$dyad))) }
}
sink()
