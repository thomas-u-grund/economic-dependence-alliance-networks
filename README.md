# Replication materials

*Trade Dependence and Conflict Beyond Neighbouring States*

Data, commented scripts, codebook, and software output for every analysis in the article and its Online Appendix.

## Contents

```
scripts/        01-40, run in numerical order (index: scripts/README.md)
raw/            source data as downloaded (Correlates of War alliances, MIDs, trade, contiguity, capabilities,
                state system, diplomatic exchange; UN ideal points; ICOW). Download links: see "Data sources" below
data/           analysis panel built from raw/ by script 02 (CSV, gzipped because of its size)
results/        output of the scripts: printed model output (*_output.txt), tables (CSV), figures, model objects
CODEBOOK.md     variables in the data files and in the analysis panel
session_info.txt   R version and package versions
```

## Requirements

R (4.6.0 used; see `session_info.txt`) and the packages `dplyr, readr, tidyr, stringr, purrr, readxl, igraph, sandwich, lmtest, survival, fixest, car, nnet, countrycode, ggplot2, patchwork`, plus `cshapes` (scripts 05, 10, 20, 25, 26, 40). Install with `install.packages(c(...))`.

## How to reproduce

All paths are relative to the repository root, so run every script from there (`Rscript scripts/NN_name.R`).

1. `bash scripts/01_unpack_data.sh` unpacks `data/tie_trade.csv` (kept gzipped because of its size). Alternatively, rebuild it from `raw/` with script 02.
2. The article's alliance definition is the Correlates of War defence pact (Type I). Set `ALLIANCE_DEF=defense` for all scripts (`export ALLIANCE_DEF=defense`); `scripts/run_all.sh` does this and runs scripts 03 to 42 in order. Leave it unset to use all four commitment types.
3. Scripts 03 to 42 build the dyad-year analysis panel in memory from `raw/` and `data/tie_trade.csv` (the shared preamble is in `scripts/03_layered_hazards.R`) and write their output to `results/`. Scripts 12, 13, and 18 only assemble earlier output and need scripts 04 to 08, 10, 11 and 15 to 17 to have run.

Example:
```
export ALLIANCE_DEF=defense
Rscript scripts/04_triad_population.R
Rscript scripts/11_claims_table.R
Rscript scripts/18_figures.R
```
Seeds are set where randomness is used (cross-validation, bootstrap). Running a script overwrites its output in `results/`. The outputs in `results/` are those behind the article, so a re-run can be checked against them.

## Where each table and figure comes from

The article analyses all dyad-years. The Online Appendix has four parts: A, main estimates and alternative estimators; B, geography and robustness; C, highly exposed relationships; D, secondary outcomes.

| Item in the article | Script | Output in `results/` |
|---|---|---|
| Figure 1 (raw onset rates by exposure, contiguity and distance), Section 5.1; Figure A1 (map of the most exposed non-contiguous pairs), Appendix C | 40 | `geography_descriptives_output.txt`, `fig_geography_descriptives.*`, `fig_geography_map.*` |
| Figure 2A (transitions), Section 5.2; Appendix A, Table 1 | 24, 30 | `population_selection_output.txt`, `fig_exposure_map_full.*` |
| Appendix A, Table 2 (actor-overlap standard errors), Section 5.6 | 29 | `fullpanel_shared_se_output.txt`, `fullpanel_shared_se.csv` |
| Table 1 (contiguity, with actor-overlap intervals), Section 5.3 | 32, 33 | `fullpanel_contiguity_implied_output.txt`, `geography_actor_overlap_output.txt` |
| Table 2 and Section 5.4 (within-between, fixed effects, longer-run exposure); Appendix A, Table 3 | 26, 31, 36, 39 | `fullpanel_within_het_output.txt`, `fullpanel_dyad_fe_output.txt`, `within_between_by_contiguity_output.txt`, `longer_run_exposure_output.txt` |
| Table 3 (controls) | 25 | `fullpanel_controls_output.txt` |
| Appendix A (embedded pairs: composition, transitions, embedding interaction, within-between) | 04, 19, 20 | `triad_population_output.txt`, `within_between_output.txt`, `heterogeneity_output.txt` |
| Appendix B, Table 4 (interactions with contiguity, distance and embedding) | 20, 26, 32 | `heterogeneity_output.txt`, `fullpanel_within_het_output.txt`, `fullpanel_contiguity_implied_output.txt` |
| Appendix B, Table 5 (major powers, regions), Section 5.3 | 22 | `distant_pairs_output.txt` |
| Figure 3 (predicted onsets); Appendix B, Tables 6 and 7 (effect sizes, influence checks) | 23, 30 | `effect_sizes_output.txt`, `predicted_onset_by_percentile.csv`, `fig_predicted_onset.*` |
| Appendix B (calendar-year effects, spline), Section 5.3 | 41 | `geography_sensitivity_output.txt` |
| Appendix B (minimum-trade rule), Section 4.3 | 42 | `trade_threshold_output.txt` |
| Appendix C, Tables 8 and 9 (upper tail, leave-outs), Section 5.3 | 34 | `upper_tail_output.txt` |
| Appendix C, Tables 10 and 11 (colonial ties, territorial claims, issue in dispute) and coverage, Sections 5.3 and 5.6 | 37 | `what_distinguishes_output.txt` |
| Appendix C, Table 12 (diplomatic exchange, UN voting, halves of the period), Section 5.3 | 38 | `diplomacy_alignment_output.txt` |
| Figure 2B (severity) and Section 5.5; Appendix D, Tables 13 and 15 (severity, initiation) | 27, 30 | `fullpanel_severity_initiator_output.txt` |
| Appendix D, Table 14 (seriousness of disputes that begin) | 26 | `fullpanel_within_het_output.txt` |

Further analyses that informed the article but are not reported in the appendix (embedded pairs in detail, embedding over time, predictors of onset and ending, held-out comparisons) are produced by scripts 05-08, 10-18, 21 and 28; their output is in `results/` as well.

## Files not included

`results/hazard_models.rds` (over 100 MB; not needed for the article).

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

