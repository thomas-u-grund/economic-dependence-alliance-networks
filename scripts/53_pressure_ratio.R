# 53: Structural pressure as the share of unbalanced triads, U/(U+B), versus separate U and B terms. Sign change, larger share held in the model.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
sink("results/pressure_ratio_output.txt", split = TRUE)
for (sg in c(-1, 1)) {
  x <- d0 %>% filter(sign == sg, !is.na(dep), !is.na(mut), !is.na(cap), any_tri) %>%
    mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), l_dur = log(dur),
           l_unbal = log1p(n_unbal), l_bal = log1p(n_bal), share_u = n_unbal / (n_unbal + n_bal), l_tot = log(n_unbal + n_bal),
           z_share = as.numeric(scale(share_u)), net = as.numeric(scale(n_unbal - n_bal)))
  cat(sprintf("\n######## %s: dyad-years in >=1 triad N=%d, sign changes=%d, mean share unbalanced %.2f\n", ifelse(sg == -1, "DISPUTES", "ALLIANCES"), nrow(x), sum(x$flip_i), mean(x$share_u)))
  ctrl <- "z_dep + z_mut + z_cap + l_dur"
  specs <- list("U and B separately" = "l_unbal + l_bal", "share unbalanced U/(U+B)" = "z_share", "share + log total triads" = "z_share + l_tot", "net U - B" = "net")
  for (nm in names(specs)) {
    m <- glm(as.formula(paste("flip_i ~", ctrl, "+", specs[[nm]])), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
    pr <- setdiff(rownames(ct), c("(Intercept)", "z_dep", "z_mut", "z_cap", "l_dur"))
    cat(sprintf("%-28s AIC %8.1f | larger share %6.3f (p=%.4f) | ", nm, AIC(m), ct["z_dep",1], ct["z_dep",4]))
    for (p in pr) cat(sprintf("%s %6.3f (p=%.4f)  ", p, ct[p,1], ct[p,4])); cat("\n") }
}
sink()
