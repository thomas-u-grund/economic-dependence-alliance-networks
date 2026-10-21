# 52: Does sign change associate more with the larger share for disputes than for alliances? Pooled, fully interacted.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
x <- d0 %>% filter(!is.na(dep), !is.na(mut), !is.na(cap)) %>%
  mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal), disp = as.integer(sign == -1))
sink("results/dispute_vs_alliance_bend_output.txt", split = TRUE)
for (y in c("flip_i", "gone_i")) {
  m <- glm(as.formula(paste(y, "~ disp * (z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal)")), x, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%s: larger share alliances %6.3f (se %.3f) | extra for disputes %6.3f (se %.3f, p=%.4f) | N=%d events=%d\n", y,
      ct["z_dep",1], ct["z_dep",2], ct["disp:z_dep",1], ct["disp:z_dep",2], ct["disp:z_dep",4], nobs(m), sum(x[[y]])))
}
sink()
