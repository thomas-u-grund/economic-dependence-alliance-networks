# 74: Assemble the 'what it survives' table for the larger-share coefficient (dispute onset / ending in T) from existing outputs. No new models.
suppressMessages({ library(dplyr); library(readr); library(stringr) })

pr <- readLines("results/proximity_in_T_output.txt"); sec <- NA; out <- list()
for (l in pr) {
  if (str_starts(l, "########")) sec <- str_match(l, "\\) ([A-Z ]+ in T[^:]*):")[, 2] %>% str_squish(); hdr <- str_match(l, "\\(([a-z-]+)\\)")[, 2]
  if (str_starts(l, "########")) { sec <- str_squish(str_match(l, "\\) (.*?): N=")[, 2]); ev <- as.integer(str_match(l, "events=(\\d+)")[, 2]); next }
  m <- str_match(l, "^\\s+([A-D][0-9]? .*?)\\s+larger share\\s+(-?[0-9.]+) \\(se ([0-9.]+)")
  if (!is.na(m[1, 1])) out[[length(out) + 1]] <- tibble(population = sec, events = ev, step = str_squish(m[1, 2]), est = as.numeric(m[1, 3]), se = as.numeric(m[1, 4]))
}
t1 <- bind_bind <- bind_rows(out) %>% group_by(population) %>% mutate(pct_of_baseline = round(100 * est / est[1]), ci_lo = round(est - 1.96 * se, 3), ci_hi = round(est + 1.96 * se, 3)) %>% ungroup()
fe <- readLines("results/dyad_fe_in_T_output.txt"); sec <- NA; o2 <- list()
for (l in fe) { if (str_starts(l, "##")) { sec <- str_squish(str_remove(l, "#+")); next }
  m <- str_match(l, "^(pooled|FE: dyad \\+ decade|FE: dyad)\\s+larger\\s+(-?[0-9.]+) \\(se ([0-9.]+)")
  if (!is.na(m[1, 1])) o2[[length(o2) + 1]] <- tibble(population = sec, step = m[1, 2], est = as.numeric(m[1, 3]), se = as.numeric(m[1, 4])) }
t2 <- bind_rows(o2) %>% group_by(population) %>% mutate(pct_of_baseline = round(100 * est / est[1]), ci_lo = round(est - 1.96 * se, 3), ci_hi = round(est + 1.96 * se, 3)) %>% ungroup()
vd <- read_csv("results/claims_table.csv", show_col_types = FALSE) %>% filter(spec %in% c("baseline", "+ volume + decade FE", "contiguous only"), outcome %in% c("dispute onset", "dispute ending")) %>% select(outcome, spec, est, ci_lo, ci_hi)
sink("results/survives_table_output.txt", split = TRUE)
cat("== Controls path (script 67): larger share, log-odds per SD ==\n"); print(as.data.frame(t1), row.names = FALSE)
cat("\n== Within-dyad (script 69) ==\n"); print(as.data.frame(t2), row.names = FALSE)
cat("\n== Volume + decade path in T (script 73) ==\n"); print(as.data.frame(vd), row.names = FALSE)
sink(); write_csv(t1, "results/survives_table_controls.csv"); write_csv(t2, "results/survives_table_within.csv")
