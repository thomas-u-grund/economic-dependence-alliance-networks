# 35: Dyad-year discrete-time hazard. Unit = signed dyad-year (dispute or alliance).
#     Flipping a dyad's sign flips the balance of EVERY triad containing it, so pressure is
#     two-sided: n_unbal = unbalanced triads a flip would fix; n_bal = balanced triads a flip would break.
#     Outcomes at t+1: sign change, dissolution (tie absent), continue. Cluster by dyad.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })

source("scripts/38_hazard_panel_lib.R")
d <- build_panel(base_signed(), "baseline") %>% mutate(end_ = flip | gone, share_unbal = ifelse(any_tri, n_unbal / (n_unbal + n_bal), NA))
cat("dyad-years:", nrow(d), " trade-valid:", sum(!is.na(d$dep)), "\n")
print(d %>% group_by(sign) %>% summarise(n = n(), flip = sum(flip), gone = sum(gone), in_any_triad = mean(any_tri)))
z <- function(x) as.numeric(scale(x))
fit <- function(data, y, rhs, label) {
  data <- data %>% filter(!is.na(dep), !is.na(mut), !is.na(cap)) %>%
    mutate(z_dep = z(dep), z_mut = z(mut), z_cap = z(cap), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal), l_dur = log(dur))
  m <- glm(as.formula(paste(y, "~", rhs)), data = data, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = data$dyad))
  cat("\n---", label, " N=", nobs(m), " events=", sum(data[[y]]), " dyads=", n_distinct(data$dyad), "\n"); print(round(ct, 4)); invisible(m)
}
d <- d %>% mutate(flip_i = as.integer(flip), gone_i = as.integer(gone), end_i = as.integer(end_))
sink("results/dyad_year_hazard_output.txt", split = TRUE)
base <- "z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal"
for (sg in c(-1, 1)) {
  nm <- ifelse(sg == -1, "DISPUTES (negative ties)", "ALLIANCES (positive ties)")
  cat("\n\n########", nm, "########\n")
  s <- d %>% filter(sign == sg)
  fit(s, "end_i",  base, "any end (sign change or dissolution)")
  fit(s, "flip_i", base, "sign change only")
  fit(s, "gone_i", base, "dissolution only")
  fit(s, "flip_i", paste(base, "+ z_dep:l_unbal + z_dep:l_bal"), "sign change, asymmetry x pressure")
  fit(s %>% filter(any_tri), "flip_i", base, "sign change, only dyads in >=1 triad")
}
sink()
