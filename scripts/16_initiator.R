# 16: Who initiates? Among dispute onsets in embedded pairs, is the more exposed state the initiator? Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })

src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
nmc <- read_csv("raw/nmc/NMCv7/abridged/NMC-70-abridged.csv", show_col_types = FALSE) %>% select(ccode, year, cinc)
mid <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% filter(hihost >= 3) %>%
  mutate(init = ifelse(dyad_rolea == 1, statea, ifelse(dyad_roleb == 1, stateb, NA)), lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>%
  filter(!is.na(init)) %>% group_by(year, lo, hi) %>% summarise(n_init = n_distinct(init), init = first(init), .groups = "drop") %>% filter(n_init == 1)
d <- d %>% mutate(tri = (n_same + n_diff) >= 1, onset = m == 0 & m1 == 1)
ev <- d %>% filter(tri, onset) %>% mutate(year1 = year + 1) %>% left_join(mid, by = c("year1" = "year", "lo", "hi")) %>%
  left_join(nmc %>% rename(lo = ccode, cinc_lo = cinc), by = c("year", "lo")) %>% left_join(nmc %>% rename(hi = ccode, cinc_hi = cinc), by = c("year", "hi"))
cat(sprintf("Onsets in embedded pairs: %d; with a single clear initiator: %d; with capability data: %d\n", sum(d$tri & d$onset), sum(!is.na(ev$init)), sum(!is.na(ev$init) & !is.na(ev$cinc_lo) & !is.na(ev$cinc_hi))))
e <- ev %>% filter(!is.na(init), !is.na(cinc_lo), !is.na(cinc_hi), cinc_lo > 0, cinc_hi > 0) %>%
  mutate(lo_init = as.integer(init == lo), dshare = (shl - shh) / 0.074, dcinc = log(cinc_lo) - log(cinc_hi), more_exp_init = as.integer((shl > shh & init == lo) | (shh > shl & init == hi)), tied = shl == shh)
sink("results/initiator_output.txt", split = TRUE)
cat(sprintf("Events used: %d (dyads %d); share in which the lower-ccode state initiates: %.3f\n", nrow(e), n_distinct(e$dyad), mean(e$lo_init)))
x <- e %>% filter(!tied); b <- binom.test(sum(x$more_exp_init), nrow(x)); cat(sprintf("Unconditional: the MORE exposed state initiates in %d of %d (%.3f; 95%% CI %.3f-%.3f)\n", sum(x$more_exp_init), nrow(x), b$estimate, b$conf.int[1], b$conf.int[2]))
cat(sprintf("Stronger state (higher CINC) initiates in %.3f of events; the more exposed state is the weaker one in %.3f\n", mean((e$cinc_lo > e$cinc_hi) == (e$lo_init == 1)), mean(((e$shl > e$shh) == (e$cinc_lo < e$cinc_hi))[!e$tied])))
m1 <- glm(lo_init ~ 0 + dshare, e, family = binomial); m2 <- glm(lo_init ~ 0 + dshare + dcinc, e, family = binomial)
for (nm in c("m1", "m2")) { m <- get(nm); ct <- coeftest(m, vcov = vcovCL(m, cluster = e$dyad)); cat(sprintf("\n%s (no intercept; dshare = (share_lo - share_hi)/0.074; dcinc = log CINC lo - hi):\n", nm)); print(round(ct, 4)) }
m3 <- glm(lo_init ~ 0 + dshare + dcinc, e %>% filter(contig == 0), family = binomial); ct3 <- coeftest(m3, vcov = vcovCL(m3, cluster = e %>% filter(contig == 0) %>% pull(dyad))); cat("\nNon-contiguous pairs only:\n"); print(round(ct3, 4)); cat(sprintf("(N = %d)\n", nobs(m3)))
sink()
