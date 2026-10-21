# 57: Descriptive diagnostic: dispute-years (signed network) with and without trade shares, by period and alliance status.
source("scripts/38_hazard_panel_lib.R")
d <- build_panel(base_signed(), "")
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>% transmute(year, lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2), allied = 1L) %>% distinct()
x <- d %>% filter(sign == -1) %>% left_join(al, by = c("year", "lo", "hi")) %>%
  mutate(allied = replace_na(allied, 0L), has_trade = !is.na(dep) & !is.na(mut) & !is.na(cap),
         per = cut(year, c(1815, 1913, 1945, 1989, 2012), labels = c("1816-1913", "1914-45", "1946-89", "1990-2011")))
sink("results/missing_trade_diagnostic_output.txt", split = TRUE)
cat("Dispute-years in signed network:", nrow(x), " with complete trade/capability data:", sum(x$has_trade), "\n")
print(x %>% group_by(per) %>% summarise(n = n(), pct_with_trade = round(100 * mean(has_trade), 1), flip_rate_with = round(100 * mean(flip_i[has_trade]), 1), flip_rate_without = round(100 * mean(flip_i[!has_trade]), 1)))
print(x %>% group_by(allied) %>% summarise(n = n(), pct_with_trade = round(100 * mean(has_trade), 1), flip_rate_with = round(100 * mean(flip_i[has_trade]), 1), flip_rate_without = round(100 * mean(flip_i[!has_trade]), 1), flips_without = sum(flip_i[!has_trade])))
cat("Overall flip rate with / without trade data:", round(100 * mean(x$flip_i[x$has_trade]), 1), round(100 * mean(x$flip_i[!x$has_trade]), 1), "\n")
sink()
