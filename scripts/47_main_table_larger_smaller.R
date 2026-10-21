# 47: Main hazard table with the larger share (maximum exposure) and smaller share entered directly.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
d <- build_panel(base_signed(), "baseline")
zs <- function(x) as.numeric(scale(x))
sink("results/main_table_larger_smaller_output.txt", split = TRUE)
for (sg in c(-1, 1)) for (y in c("flip_i", "gone_i")) {
  x <- d %>% filter(sign == sg, !is.na(dep), !is.na(mut), !is.na(cap)) %>%
    mutate(z_mx = zs(mx), z_mut = zs(mut), z_cap = zs(cap), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
  m <- glm(as.formula(paste(y, "~ z_mx + z_mut + z_cap + l_dur + l_unbal + l_bal")), x, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("\n-- %s / %s  N=%d ev=%d dyads=%d\n", ifelse(sg == -1, "DISPUTES", "ALLIANCES"), y, nobs(m), sum(x[[y]]), n_distinct(x$dyad)))
  print(round(ct, 4))
}
cat("\nSD larger share:", sd(d$mx, na.rm = TRUE), " SD smaller:", sd(d$mut, na.rm = TRUE), " SD asym:", sd(d$dep, na.rm = TRUE), "\n")
sink()
