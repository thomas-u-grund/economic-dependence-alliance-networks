# 69: Within-dyad (dyad fixed effects) estimates for dispute onset in the embedded population T. ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(fixest) })

src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), tri = (n_same + n_diff) >= 1, l_dur_m = log(dur_m))
sink("results/dyad_fe_in_T_output.txt", split = TRUE)
for (spec in list(list("onset_m", d %>% filter(tri, m == 0), "DISPUTE ONSET in T"), list("end_m", d %>% filter(tri, m == 1), "DISPUTE ENDING in T"))) {
  x <- spec[[2]]; y <- spec[[1]]; cat(sprintf("\n## %s: N=%d events=%d\n", spec[[3]], nrow(x), sum(x[[y]])))
  f0 <- as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur_m + l_same + l_diff")); m0 <- feglm(f0, x, family = "binomial", cluster = ~dyad); co <- coeftable(m0)
  cat(sprintf("pooled              larger %6.3f (se %.3f, p=%.4f)\n", co["z_mx",1], co["z_mx",2], co["z_mx",4]))
  for (fe in c("dyad", "dyad + decade")) { m <- tryCatch(feglm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + l_dur_m + l_same + l_diff |", fe)), x, family = "binomial", cluster = ~dyad), error = function(e) NULL)
    if (!is.null(m)) { co <- coeftable(m); cat(sprintf("FE: %-14s larger %6.3f (se %.3f, p=%.4f) | N used=%d | split %6.3f (p=%.4f) common %6.3f (p=%.4f)\n", fe, co["z_mx",1], co["z_mx",2], co["z_mx",4], m$nobs, co["l_diff",1], co["l_diff",4], co["l_same",1], co["l_same",4])) } }
}
sink()
