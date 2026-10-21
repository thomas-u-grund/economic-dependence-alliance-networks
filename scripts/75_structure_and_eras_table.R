# 75: Assemble the structure benchmark (script 68) and era stability (script 70) tables from existing outputs. No new models.
suppressMessages({ library(dplyr); library(readr); library(stringr) })

pr <- readLines("results/predictors_in_T_output.txt"); out <- list(); sec <- NA
for (l in pr) {
  if (str_starts(l, "########")) { sec <- str_squish(str_match(l, "## (.*?):")[, 2]); next }
  m <- str_match(l, "^(z_\\w+|contig)\\s+joint coef.*?\\s(-?[0-9.]+) \\(se ([0-9.]+), p=([0-9.]+)\\)(?: \\| single-predictor held-out AUC ([0-9.]+))?")
  if (!is.na(m[1, 1])) out[[length(out) + 1]] <- tibble(outcome = sec, predictor = m[1, 2], coef_per_SD = as.numeric(m[1, 3]), se = as.numeric(m[1, 4]), p = as.numeric(m[1, 5]), AUC_single = as.numeric(m[1, 6]))
}
t1 <- bind_rows(out) %>% mutate(predictor = recode(predictor, z_dur = "log duration", z_mx = "larger trade share (exposure)", z_mn = "smaller trade share", z_cap = "capability asymmetry", z_same = "common partners (log1p)", z_diff = "split partners (log1p)", contig = "contiguity (0/1)"))
gains <- tibble(line = pr[str_detect(pr, "^Incremental")]) 
er <- readLines("results/embedding_over_time_output.txt")
orr <- tibble(l = er[str_detect(er, "T log-odds")]) %>% mutate(era = str_match(l, "^\\s+(\\S+)")[, 2], logodds = as.numeric(str_match(l, "log-odds\\s+(-?[0-9.]+)")[, 2]), se = as.numeric(str_match(l, "se ([0-9.]+)")[, 2]), OR = as.numeric(str_match(l, "OR ([0-9.]+)")[, 2])) %>% select(-l)
bal <- tibble(l = er[str_detect(er, "^\\s+l_(diff|same)\\s+\\d{4}")][1:8]) %>% mutate(term = str_match(l, "^\\s+(\\S+)")[, 2], era = str_match(l, "^\\s+\\S+\\s+(\\S+)")[, 2], coef = as.numeric(str_match(l, "\\s(-?[0-9.]+) \\(se")[, 2]), se = as.numeric(str_match(l, "se ([0-9.]+)")[, 2])) %>% select(-l)
rate <- er[grepl("^[1-4] ", er)][1:4]
sink("results/structure_and_eras_table_output.txt", split = TRUE)
cat("== Standardised coefficients (all predictors scaled within the same sample), embedded pairs ==\n"); print(as.data.frame(t1), row.names = FALSE)
cat("\n== Held-out information gains (script 68) ==\n"); cat(pr[str_detect(pr, "^(Incremental|Exposure \\+|\\+ third|Structure)")], sep = "\n")
cat("\n== Era stability (script 70): embedded vs not, odds ratio with controls + decade FE ==\n"); print(as.data.frame(orr), row.names = FALSE)
cat("\n== Era-specific balance slopes (onset, within T) ==\n"); print(as.data.frame(bal), row.names = FALSE)
cat("\n== Onset rate by era (rows from script 70 part a) ==\n"); cat(rate, sep = "\n")
sink(); write_csv(t1, "results/structure_benchmark.csv"); write_csv(orr, "results/era_embedded_OR.csv"); write_csv(bal, "results/era_balance_slopes.csv")
