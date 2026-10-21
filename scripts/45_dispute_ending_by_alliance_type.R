# 45: Dispute ending (sign change or end) within alliances, by alliance type and by bloc vs non-bloc alliance.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>%
  mutate(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>%
  group_by(year, lo, hi) %>% summarise(has_def = max(defense == 1), has_oth = max(neutrality == 1 | nonaggression == 1 | entente == 1), .groups = "drop") %>%
  mutate(atype = ifelse(has_def == 1, "defense pact", "non-defense only (neutrality/non-aggression/entente)"))
mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE)
big <- mem %>% group_by(version4id) %>% summarise(n = n_distinct(ccode)) %>% filter(n >= 8) %>% pull(version4id)
d <- d0 %>% filter(sign == -1) %>% inner_join(al, by = c("year", "lo", "hi")) %>% mutate(end_i = as.integer(flip | gone)) %>%
  filter(!is.na(dep), !is.na(mut), !is.na(cap)) %>%
  mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
# large-bloc membership of both states (any of the >=8-member alliances)
mb <- mem %>% filter(version4id %in% big) %>% transmute(id = version4id, ccode, y0 = mem_st_year, y1 = ifelse(is.na(mem_end_year) | mem_end_year == 0, 9999, mem_end_year))
inbloc <- rep(FALSE, nrow(d))
for (id in big) { m <- mb[mb$id == id, ]
  f <- function(cc, yr) sapply(seq_along(cc), function(i) any(m$ccode == cc[i] & m$y0 <= yr[i] & m$y1 >= yr[i]))
  inbloc <- inbloc | (f(d$lo, d$year) & f(d$hi, d$year)) }
d$bloc_all <- ifelse(inbloc, "both in a large (>=8 member) alliance", "other alliance")
sink("results/dispute_ending_by_alliance_type_output.txt", split = TRUE)
cat("dispute-years with alliance in force:", nrow(d), "\n")
rhs <- "z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal"
for (v in c("atype", "bloc_all")) for (g in sort(unique(d[[v]]))) {
  s <- d[d[[v]] == g, ]
  if (sum(s$end_i) < 20) { cat(sprintf("%-60s events=%d dyads=%d  (too few)\n", g, sum(s$end_i), n_distinct(s$dyad))); next }
  m <- glm(as.formula(paste("end_i ~", rhs)), s, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = s$dyad))
  cat(sprintf("%-60s dep=%6.3f (se %.3f, p=%.4f) | N=%d ends=%d dyads=%d\n", g, ct["z_dep",1], ct["z_dep",2], ct["z_dep",4], nobs(m), sum(s$end_i), n_distinct(s$dyad)))
}
sink()
