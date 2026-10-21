# 93: FULL PANEL: dyad fixed-effects logit for dispute onset and ending (uses only dyads with variation in the outcome). Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(fixest) })
src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), decade = as.factor(decade))
sink("results/fullpanel_dyad_fe_output.txt", split = TRUE)
for (spec in list(list("DISPUTE ONSET", "onset_m", d %>% filter(m == 0), "dur_m"), list("DISPUTE ENDING", "end_m", d %>% filter(m == 1), "dur_m"))) {
  x <- spec[[3]]; x$l_dur <- log(x[[spec[[4]]]]); cat(sprintf("\n## %s, full panel\n", spec[[1]]))
  m0 <- feglm(as.formula(paste(spec[[2]], "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff")), x, family = "binomial", cluster = ~dyad); c0 <- coeftable(m0)
  m1 <- feglm(as.formula(paste(spec[[2]], "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff | dyad")), x, family = "binomial", cluster = ~dyad); c1 <- coeftable(m1)
  m2 <- feglm(as.formula(paste(spec[[2]], "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff | dyad + decade")), x, family = "binomial", cluster = ~dyad); c2 <- coeftable(m2)
  cat(sprintf("pooled              larger %6.3f (se %.3f, p=%.4f)\n", c0["z_mx", 1], c0["z_mx", 2], c0["z_mx", 4]))
  cat(sprintf("FE: dyad            larger %6.3f (se %.3f, p=%.4f) | N used=%d\n", c1["z_mx", 1], c1["z_mx", 2], c1["z_mx", 4], m1$nobs))
  cat(sprintf("FE: dyad + decade   larger %6.3f (se %.3f, p=%.4f) | N used=%d\n", c2["z_mx", 1], c2["z_mx", 2], c2["z_mx", 4], m2$nobs)) }
sink()
