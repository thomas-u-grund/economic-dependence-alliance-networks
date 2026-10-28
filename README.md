# Replication materials

*Economic Dependence and Interstate Conflict: Geography and Differences Between Relationships*

Data, commented scripts, codebook, and software output for every analysis in the article and its Online Appendix.

## Contents

```
scripts/        01-40, run in numerical order (index: scripts/README.md)
raw/            source data as downloaded (Correlates of War alliances, MIDs, trade, contiguity, capabilities,
                state system, diplomatic exchange; UN ideal points; ICOW; ATOP 5.1, Maddison GDP and Polity 5, not
                used by the scripts). Download links: see "Data sources" below
data/           intermediate files built from raw/ by script 02 (CSV); two large ones are gzipped
results/        output of the scripts: printed model output (*_output.txt), tables (CSV), figures, model objects
CODEBOOK.md     variables in the data files and in the analysis panel
session_info.txt   R version and package versions
```

## Requirements

R (4.6.0 used; see `session_info.txt`) and the packages `dplyr, readr, tidyr, stringr, purrr, readxl, igraph, sandwich, lmtest, survival, fixest, car, nnet, countrycode, ggplot2, patchwork`, plus `cshapes` (scripts 05, 10, 20, 25, 26, 40). Install with `install.packages(c(...))`.

## How to reproduce

All paths are relative to the repository root, so run every script from there (`Rscript scripts/NN_name.R`).

1. `bash scripts/01_unpack_data.sh` unpacks `data/tie_trade.csv` and `data/triads_all_years.csv` (kept gzipped because of their size). Alternatively, rebuild them from `raw/` with script 02.
2. The article's alliance definition is the Correlates of War defence pact (Type I). Set `ALLIANCE_DEF=defense` for all scripts (`export ALLIANCE_DEF=defense`); `scripts/run_all.sh` does this and runs scripts 03 to 42 in order. Leave it unset to use all four commitment types.
3. Scripts 03 to 42 build the dyad-year analysis panel in memory from `raw/` and `data/tie_trade.csv` (the shared preamble is in `scripts/03_layered_hazards.R`) and write their output to `results/`. Scripts 12, 13, and 18 only assemble earlier output and need scripts 04 to 11 and 15 to 17 to have run.

Example:
```
export ALLIANCE_DEF=defense
Rscript scripts/04_triad_population.R
Rscript scripts/11_claims_table.R
Rscript scripts/18_figures.R
```
Seeds are set where randomness is used (cross-validation, bootstrap). Running a script overwrites its output in `results/`. The outputs in `results/` are those behind the article, so a re-run can be checked against them.

## Where each table and figure comes from

The principal population of the article is all dyad-years; the results for embedded pairs alone are in Appendix A.

| Item in the article | Script | Output in `results/` |
|---|---|---|
| Appendix A, Table 2 (composition by embedding and alliance) | 04 | `triad_population_output.txt` |
| Figure 1 (raw onset rates by exposure, contiguity and distance) and Figure A2 (map of the most exposed non-contiguous pairs), Section 5.1 and Appendix H | 40 | `geography_descriptives_output.txt`, `fig_geography_descriptives.*`, `fig_geography_map.*` |
| Calendar-year effects, trimmed and spline checks on the geographical interaction, Section 5.3 and Appendix G | 41 | `geography_sensitivity_output.txt` |
| Minimum-trade rule (none, 0.01% and 0.05% of recorded world trade), Section 4.3 and Appendix H | 42 | `trade_threshold_output.txt` |
| Figure 2A (transitions, all dyad-years), Section 5.2, Appendix D | 24, 29, 30 | `population_selection_output.txt`, `fullpanel_shared_se_output.txt`, `fig_exposure_map_full.*` |
| Figure 2B (severity), Sections 5.2 and 5.5, Appendices C and D | 27, 30 | `fullpanel_severity_initiator_output.txt` |
| Table 1 and Section 5.3 (contiguity, distance), Appendix E | 26, 32, 33 | `fullpanel_within_het_output.txt`, `fullpanel_contiguity_implied_output.txt` |
| Section 5.3, Appendix F (major powers, regions) | 22 | `distant_pairs_output.txt` |
| Figure 3 (predicted onsets), Appendix G | 23, 30 | `effect_sizes_output.txt`, `predicted_onset_by_percentile.csv`, `fig_predicted_onset.*` |
| Section 5.3 (colonial ties, territorial claims, issue in dispute, upper tail by era, coverage), Appendix H | 37 | `what_distinguishes_output.txt` |
| Section 5.4 (longer-run within-pair exposure), Appendix E | 39 | `longer_run_exposure_output.txt` |
| Section 5.3 (diplomatic exchange, UN voting alignment, halves of the period), Appendix H | 38 | `diplomacy_alignment_output.txt` |
| Table 2 and Section 5.4 (within-between, fixed effects), Appendix E | 26, 31, 36 | `fullpanel_within_het_output.txt`, `fullpanel_dyad_fe_output.txt`, `within_between_by_contiguity_output.txt` |
| Table 3 (controls) | 25 | `fullpanel_controls_output.txt` |
| Appendix D, Table 12 (predictors, held-out comparisons) | 28 | `fullpanel_predictors_output.txt` |
| Figure A1 (periods), Appendix D | 08, 17, 18 | `embedding_over_time_output.txt`, `era_bins_rates.csv`, `era_bins_shares.csv`, `fig_eras.*` |
| Appendix A, Section 5.2 (embedded pairs only) | 04, 06, 11, 15, 19 | `triad_population_output.txt`, `predictors_in_T_output.txt`, `claims_table_output.txt`, `hostility_gradient_output.txt`, `within_between_output.txt` |
| not in the current appendix (removed 2026-10-27) (disputes inside unbalanced triads) | 09, 10 | `negative_tie_resolution_output.txt` |
| not in the current appendix (removed 2026-10-27), Sections 5.2 and 5.6 (lagged embedding and partner counts) | 24, 35 | `population_selection_output.txt`, `lagged_network_output.txt` |
| Appendix B, Section 5.6 (actor-overlap intervals) | 33 | `geography_actor_overlap_output.txt` |
| Appendix H (upper tail of exposure) | 34 | `upper_tail_output.txt` |

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

`raw/` also holds files that the scripts of the article do not read: ATOP 5.1 (https://www.atopdata.org/data.html), the Maddison Project Database 2023 (https://www.rug.nl/ggdc/historicaldevelopment/maddison/), and Polity 5 (https://www.systemicpeace.org/inscrdata.html). Where a source page offers several versions, use the version in the table.
