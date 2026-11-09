# run_all.R: reproduces every table, figure and number in the article and the Online Appendix.
# Run from the repository root:   Rscript run_all.R
# Outputs overwrite the files in results/, which hold the versions behind the article, so a re-run can be
# compared with them (for example with `git diff results/`). A timing log is written to results/run_all_log.txt.
# The analysis panel data/tie_trade.csv.gz is included; scripts/build_data.R rebuilds it from raw/ if needed.
Sys.setenv(ALLIANCE_DEF = "defense")
if (!file.exists("data/tie_trade.csv.gz")) source("scripts/build_data.R")
steps <- c("scripts/1_descriptives.R", "scripts/2_main_estimates.R", "scripts/3_geography.R", "scripts/4_exposed_pairs.R",
           "scripts/5_within_between.R", "scripts/6_severity_initiation.R", "scripts/7_figures.R")
log <- "results/run_all_log.txt"
cat(format(Sys.time()), "\n", file = log)
for (f in steps) {
  t0 <- Sys.time()
  source(f)
  cat(sprintf("%-32s done in %4.0f s\n", f, as.numeric(difftime(Sys.time(), t0, units = "secs"))), file = log, append = TRUE)
}
cat("\n", capture.output(sessionInfo()), sep = "\n", file = log, append = TRUE)
cat("All scripts finished. Outputs are in results/.\n")
