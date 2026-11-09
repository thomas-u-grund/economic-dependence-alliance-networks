# Replication materials

*Trade Dependence and Conflict Beyond Neighbouring States*

Data, commented scripts, codebook, and software output for every analysis in the article and its Online Appendix.

## Contents

```
run_all.R       master script: reproduces every table, figure and number in the article and the Online Appendix
scripts/        the analysis, one script per part of the article (index: scripts/README.md)
raw/            source data as downloaded (Correlates of War alliances, MIDs, trade, contiguity, capabilities,
                state system, diplomatic exchange; UN ideal points; ICOW). Download links: see "Data sources" below
data/           analysis panel tie_trade.csv.gz, built from raw/ by scripts/build_data.R
results/        output of the scripts: printed model output (*_output.txt), tables (CSV), figures
CODEBOOK.md     variables in the data files and in the analysis panel
session_info.txt   R version and package versions
```

## Requirements

R (4.6.0 used; see `session_info.txt`) and the packages `dplyr, readr, tidyr, stringr, sandwich, lmtest, fixest, zoo, ggplot2, patchwork, sf, cshapes, rnaturalearth`. Install with
```
install.packages(c("dplyr", "readr", "tidyr", "stringr", "sandwich", "lmtest", "fixest", "zoo", "ggplot2", "patchwork", "sf", "cshapes", "rnaturalearth"))
```

## How to reproduce

From the repository root:
```
Rscript run_all.R
```
This runs the seven analysis scripts in order and overwrites the files in `results/`. The committed files in `results/` are the ones behind the article, so a re-run can be checked against them (`git diff results/`). A full run takes about 15 minutes on a laptop; timings are written to `results/run_all_log.txt`.

Each analysis script can also be run on its own (`Rscript scripts/3_geography.R`), except `7_figures.R`, which reads output of scripts 2 and 6. Every script sources `scripts/setup.R`, which loads the packages and builds the dyad-year panel from `raw/` and `data/tie_trade.csv.gz` (`scripts/panel.R`). The article's alliance definition is the Correlates of War defence pact (Type I) and is the default; set the environment variable `ALLIANCE_DEF=all` to count all four commitment types instead. `scripts/build_data.R` rebuilds `data/tie_trade.csv.gz` from `raw/`; it is not needed to reproduce the results. No analysis uses random numbers.

## Where each table and figure comes from

The article analyses all dyad-years. The Online Appendix has four parts: A, main estimates and alternative estimators; B, geography and robustness; C, highly exposed relationships; D, secondary outcomes.

| Item in the article | Script | Output in `results/` |
|---|---|---|
| Figure 1 (raw onset rates by exposure, contiguity and distance), Section 5.1; Figure A1 (map of the most exposed non-contiguous pairs), Appendix C | 1_descriptives | `geography_descriptives_output.txt`, `fig_geography_descriptives.*`, `fig_geography_map.*` |
| Section 5.2 (transitions); Appendix A, Table 1 | 2_main_estimates | `population_selection_output.txt` |
| Table 3 (controls) | 2_main_estimates | `fullpanel_controls_output.txt` |
| Appendix A, Table 2 (actor-overlap standard errors), Section 5.6 | 2_main_estimates | `fullpanel_shared_se_output.txt`, `fullpanel_shared_se.csv` |
| Table 1 (contiguity, with actor-overlap intervals), Section 5.3; Appendix B, Table 4 | 3_geography | `fullpanel_contiguity_implied_output.txt`, `geography_actor_overlap_output.txt`, `heterogeneity_output.txt` |
| Appendix B, Table 5 (major powers, regions), Section 5.3 | 3_geography | `distant_pairs_output.txt` |
| Appendix B, Tables 6 and 7 (effect sizes, influence checks) | 3_geography | `effect_sizes_output.txt` |
| Appendix B (calendar-year effects, spline), Section 5.3 | 3_geography | `geography_sensitivity_output.txt` |
| Appendix B (minimum-trade rule), Section 4.3 | 3_geography | `trade_threshold_output.txt` |
| Appendix C, Tables 8 and 9 (upper tail, leave-outs), Section 5.3 | 4_exposed_pairs | `upper_tail_output.txt` |
| Appendix C, Tables 10 and 11 (colonial ties, territorial claims, issue in dispute) and coverage, Sections 5.3 and 5.6 | 4_exposed_pairs | `what_distinguishes_output.txt` |
| Appendix C, Table 12 (diplomatic exchange, UN voting, halves of the period), Section 5.3 | 4_exposed_pairs | `diplomacy_alignment_output.txt` |
| Table 2 and Section 5.4 (within-between, fixed effects, longer-run exposure); Appendix A, Table 3; Appendix B, Table 4 (distance) | 5_within_between | `within_between_by_contiguity_output.txt`, `fullpanel_within_het_output.txt`, `fullpanel_dyad_fe_output.txt`, `longer_run_exposure_output.txt` |
| Appendix A (embedded pairs: composition, transitions, within-between) | 5_within_between | `triad_population_output.txt`, `within_between_output.txt` |
| Appendix D, Table 14 (seriousness of disputes that begin) | 5_within_between | `fullpanel_within_het_output.txt` |
| Section 5.5; Appendix D, Tables 13 and 15 (severity, initiation) | 6_severity_initiation | `fullpanel_severity_initiator_output.txt` |
| Figure 2 (transitions and severity) and Figure 3 (predicted onsets) | 7_figures | `fig_exposure_map_full.*`, `fig_predicted_onset.*`, `predicted_onset_by_percentile.csv` |

## Data sources

All source files are in `raw/`, unmodified, as downloaded. To obtain them again, use the links below. Full references are in the article.

| Data | Version | Folder in `raw/` | Source |
|---|---|---|---|
| Correlates of War (COW) formal alliances | 4.1 (dyadic yearly) | `alliances/` | https://correlatesofwar.org/data-sets/formal-alliances/ |
| COW Dyadic Militarized Interstate Disputes (MIDs) | 4.03 | `mids/` | https://correlatesofwar.org/data-sets/MIDs/ |
| COW bilateral trade | 4.0 | `trade/` | https://correlatesofwar.org/data-sets/bilateral-trade/ |
| COW direct contiguity | 3.20 | `contiguity/` | https://correlatesofwar.org/data-sets/direct-contiguity/ |
| COW National Material Capabilities | 7 | `nmc/` | https://correlatesofwar.org/data-sets/national-material-capabilities/ |
| COW state system membership | 2024 | `states/` | https://correlatesofwar.org/data-sets/state-system-membership/ |
| COW diplomatic exchange | 2006.1 | `diplex/dipex.csv` | https://correlatesofwar.org/data-sets/diplomatic-exchange/ |
| UN General Assembly ideal points (Bailey, Strezhnev and Voeten) | Dataverse release | `diplex/unideal.csv` | https://doi.org/10.7910/DVN/LEJUQZ |
| Issue Correlates of War (ICOW) territorial claims | dyad-year file | `icow/terr/` | https://paulhensel.org/icowterr.html |
| ICOW colonial history | 1.1 | `icow/colhist/` | https://paulhensel.org/icowcol.html |
| CShapes (capital locations and distances) | 2.0, through the R package `cshapes` | not in `raw/` | https://icr.ethz.ch/data/cshapes/ and https://cran.r-project.org/package=cshapes |

