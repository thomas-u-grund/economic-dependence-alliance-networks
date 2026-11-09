# Scripts

Run everything with `Rscript run_all.R` from the repository root, or one script at a time with `Rscript scripts/<name>.R`. Each analysis script consists of self-contained blocks; each block writes the `results/` file named in its `sink()` or `ggsave()` call.

| Script | Content | Where used |
|---|---|---|
| `setup.R` | Packages, alliance definition (`ALLIANCE_DEF`, default `defense`), and the analysis panel; sourced by every analysis script | all |
| `panel.R` | Builds the dyad-year panel `d` from `raw/` and `data/tie_trade.csv.gz` (variables: CODEBOOK.md, section 3) | all |
| `1_descriptives.R` | Onset rates by exposure, contiguity and distance; map of the most exposed non-contiguous pairs | Figure 1, Section 5.1; Figure A1, Appendix C |
| `2_main_estimates.R` | Transitions on the full panel; controls path; actor-overlap standard errors | Section 5.2, Table 3; Appendix A Tables 1 and 2; Section 5.6 |
| `3_geography.R` | Contiguity interaction with implied effects and actor-overlap intervals; interactions with distance and embedding; major powers and regions; effect sizes and influence checks; calendar-year effects and spline; minimum-trade rules | Table 1, Sections 4.3 and 5.3; Appendix B Tables 4 to 7 and checks |
| `4_exposed_pairs.R` | Upper tail and leave-outs; colonial ties, territorial claims, issue in dispute, coverage; diplomatic exchange and UN voting | Sections 5.3 and 5.6; Appendix C Tables 8 to 12 |
| `5_within_between.R` | Within-between decomposition by contiguity; distance interaction and seriousness given onset; dyad fixed-effects logit; longer-run exposure; embedded pairs | Table 2, Section 5.4; Appendix A Table 3 and embedded pairs; Appendix B Table 4; Appendix D Table 14 |
| `6_severity_initiation.R` | Severity of onsets and initiation | Section 5.5; Appendix D Tables 13 and 15 |
| `7_figures.R` | Figures 2 and 3 (reads output of scripts 2 and 6) | Figures 2 and 3 |
| `build_data.R` | Rebuilds `data/tie_trade.csv.gz` from `raw/` (not needed to reproduce the results) | data |
