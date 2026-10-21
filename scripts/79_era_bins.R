# 79: Onset rates for embedded and other pairs by 20-year period (data for Figure 3). Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr) })

src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(tri = (n_same + n_diff) >= 1, onset = as.integer(m == 0 & m1 == 1),
  per = cut(year, c(1869, 1889, 1909, 1929, 1949, 1969, 1989, 2011), labels = c("1870-89", "1890-1909", "1910-29", "1930-49", "1950-69", "1970-89", "1990-2010")))
r <- d %>% filter(m == 0) %>% group_by(per, tri) %>% summarise(n = n(), onsets = sum(onset), .groups = "drop") %>% group_by(per, tri) %>%
  mutate(rate = 1000 * onsets / n, lo = 1000 * binom.test(onsets, n)$conf.int[1], hi = 1000 * binom.test(onsets, n)$conf.int[2]) %>% ungroup()
sh <- d %>% group_by(per) %>% summarise(pct_dyad_years_embedded = 100 * mean(tri), .groups = "drop") %>%
  left_join(d %>% filter(m == 0, onset == 1) %>% group_by(per) %>% summarise(pct_onsets_embedded = 100 * mean(tri), onsets = n(), .groups = "drop"), by = "per")
write_csv(r, "results/era_bins_rates.csv"); write_csv(sh, "results/era_bins_shares.csv")
print(as.data.frame(r)); print(as.data.frame(sh))
