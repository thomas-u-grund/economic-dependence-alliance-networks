# 44: Does exposure predict a militarized dispute ENDING at all, irrespective of alliance status,
#     and does that depend on an alliance being in force? (disputes only; outcome = sign change OR ending)
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>%
  transmute(year, lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2), allied = 1L) %>% distinct()
d <- d0 %>% filter(sign == -1) %>% left_join(al, by = c("year", "lo", "hi")) %>% mutate(allied = replace_na(allied, 0L), end_i = as.integer(flip | gone)) %>%
  filter(!is.na(dep), !is.na(mut), !is.na(cap)) %>%
  mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)),
         z_max = as.numeric(scale(mx)), z_min = as.numeric(scale(mut)),
         l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
BASE <- "z_mut + z_cap + l_dur + l_unbal + l_bal"
sink("results/dispute_ending_by_alliance_output.txt", split = TRUE)
cat("dispute-years", nrow(d), " allied", sum(d$allied), " ends", sum(d$end_i), " ends|allied", sum(d$end_i[d$allied == 1]), " ends|not allied", sum(d$end_i[d$allied == 0]), "\n")
show <- function(m, vars, lab) { ct <- coeftest(m, vcov = vcovCL(m, cluster = d$dyad[as.integer(rownames(m$model))]))
  cat(sprintf("%-52s", lab)); for (v in vars) cat(sprintf(" | %s=%6.3f (se %.3f, p=%.4f)", v, ct[v,1], ct[v,2], ct[v,4])); cat(sprintf(" | N=%d\n", nobs(m))) }
show(glm(end_i ~ z_dep + allied + z_mut + z_cap + l_dur + l_unbal + l_bal, d, family = binomial), c("z_dep", "allied"), "Dispute ends (any), + allied indicator")
show(glm(end_i ~ z_dep * allied + z_mut + z_cap + l_dur + l_unbal + l_bal, d, family = binomial), c("z_dep", "allied", "z_dep:allied"), "Dispute ends (any), asymmetry x allied")
for (a in c(1, 0)) { s <- d %>% filter(allied == a)
  m <- glm(as.formula(paste("end_i ~ z_dep +", BASE)), s, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = s$dyad))
  cat(sprintf("Dispute ends (any), %-18s dep=%6.3f (se %.3f, p=%.4f) | N=%d ev=%d dyads=%d\n", ifelse(a == 1, "alliance in force", "no alliance"), ct["z_dep",1], ct["z_dep",2], ct["z_dep",4], nobs(m), sum(s$end_i), n_distinct(s$dyad))) }
cat("\n-- same with the two shares in place of asymmetry --\n")
for (a in c(1, 0)) { s <- d %>% filter(allied == a)
  m <- glm(end_i ~ z_max + z_min + z_cap + l_dur + l_unbal + l_bal, s, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = s$dyad))
  cat(sprintf("%-18s max=%6.3f (p=%.4f) min=%6.3f (p=%.4f)\n", ifelse(a == 1, "alliance in force", "no alliance"), ct["z_max",1], ct["z_max",4], ct["z_min",1], ct["z_min",4])) }
cat("\n-- sign change vs ending by composition: where does the dispute end? (among dispute endings) --\n")
e <- d %>% filter(end_i == 1) %>% mutate(in_alliance = as.integer(flip))
m <- glm(in_alliance ~ z_dep + allied + z_mut + z_cap + l_dur + l_unbal + l_bal, e, family = binomial)
ct <- coeftest(m, vcov = vcovCL(m, cluster = e$dyad)); cat(sprintf("P(ends with alliance retained | ends): dep=%6.3f (se %.3f, p=%.4f) allied=%6.3f N=%d\n", ct["z_dep",1], ct["z_dep",2], ct["z_dep",4], ct["allied",1], nobs(m)))
sink()
