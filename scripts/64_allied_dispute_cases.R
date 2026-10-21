# 64: What are the disputes between allies? Alliance type and the dyads that contribute most.
suppressMessages({ library(dplyr); library(readr); library(tidyr) })

al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>% filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>%
  mutate(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>% group_by(year, lo, hi) %>% summarise(def = max(defense == 1), neu = max(neutrality == 1), nag = max(nonaggression == 1), ent = max(entente == 1), .groups = "drop")
mi <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>% filter(hihost >= 3) %>% distinct(year, lo, hi, hihost)
x <- mi %>% inner_join(al, by = c("year", "lo", "hi"))
nm <- read_csv("raw/states/States2024/statelist2024.csv", show_col_types = FALSE); names(nm) <- tolower(names(nm))
cn <- nm %>% distinct(ccode, .keep_all = TRUE) %>% select(ccode, stateabb)
sink("results/allied_dispute_cases_output.txt", split = TRUE)
cat("Dispute-years (hihost>=3) in allied dyads:", nrow(x), " of all dispute-years:", nrow(mi), "\n")
cat("By strongest covering COW commitment: defense pact", sum(x$def == 1), "| neutrality (no defense)", sum(x$def == 0 & x$neu == 1), "| non-aggression only", sum(x$def == 0 & x$neu == 0 & x$nag == 1), "| entente only", sum(x$def == 0 & x$neu == 0 & x$nag == 0 & x$ent == 1), "\n")
cat("By hostility level: ", paste(names(table(x$hihost)), table(x$hihost), collapse = "; "), "\n")
top <- x %>% group_by(lo, hi) %>% summarise(years = n(), first = min(year), last = max(year), .groups = "drop") %>% arrange(desc(years)) %>% left_join(cn, by = c("lo" = "ccode")) %>% rename(a = stateabb) %>% left_join(cn, by = c("hi" = "ccode")) %>% rename(b = stateabb)
cat("Distinct allied dyads with a dispute-year:", nrow(top), "\n"); print(as.data.frame(head(top[, c("a", "b", "years", "first", "last")], 15)))
cat("Share of allied dispute-years contributed by the top 10 dyads:", round(100 * sum(head(top$years, 10)) / nrow(x), 1), "%\n")
sink()
