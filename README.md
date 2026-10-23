# Replication materials

*Economic Dependence and Interstate Conflict: Geography and Differences Between Relationships*

Data, commented scripts, codebook, and software output for every analysis in the article and its Online Appendix.

## Contents

```
scripts/        01-102, one script per analysis (index: scripts/README.md)
raw/            source data as downloaded (Correlates of War alliances, MIDs, trade, contiguity, capabilities;
                ATOP 5.1; Penn World Table / Maddison GDP; Polity 5; state system)
data/           intermediate files built from raw/ by scripts 01-05 (CSV); two large ones are gzipped
results/        output of the scripts: printed model output (*_output.txt), tables (CSV), figures, model objects
CODEBOOK.md     variables in the data files and in the analysis panel
session_info.txt   R version and package versions
```

## Requirements

R (4.6.0 used; see `session_info.txt`) and the packages `dplyr, readr, tidyr, stringr, purrr, readxl, igraph, sandwich, lmtest, survival, fixest, car, nnet, countrycode, ggplot2, patchwork`, plus `cshapes` (scripts 63, 67, 72, 73, 82, 87, 88, 92, 102) and `pwt10` (script 55) from CRAN. Install with `install.packages(c(...))`.

## How to reproduce

All paths are relative to the repository root, so run every script from there (`Rscript scripts/NN_name.R`).

1. `bash scripts/00_unpack_data.sh` unpacks `data/tie_trade.csv` and `data/triads_all_years.csv` (kept gzipped because of their size). Alternatively, rebuild them from `raw/` with scripts 01 and 02.
2. The article's alliance definition is the Correlates of War defence pact (Type I). Set `ALLIANCE_DEF=defense` for all scripts from 35 onward (`export ALLIANCE_DEF=defense`). Leave it unset to use all four commitment types. `scripts/run_all_defense.sh` re-runs scripts 35 to 64 with the switch on. Scripts 66 to 94 expect it on.
3. Scripts 66 to 102 build the dyad-year analysis panel in memory from `raw/` and `data/tie_trade.csv` (the shared preamble is in `scripts/49_layered_hazards.R`) and write their output to `results/`. Scripts 74, 75, and 80 only assemble earlier output and need 66 to 73 and 77 to 79 to have run.

Example:
```
export ALLIANCE_DEF=defense
Rscript scripts/66_triad_population.R
Rscript scripts/73_claims_table.R
Rscript scripts/80_figures.R
```
Seeds are set where randomness is used (cross-validation, bootstrap). Running a script overwrites its output in `results/`. The outputs in `results/` are those behind the article, so a re-run can be checked against them.

## Where each table and figure comes from

The principal population of the article is all dyad-years; the results for embedded pairs alone are in Appendix A.

| Item in the article | Script | Output in `results/` |
|---|---|---|
| Appendix A, Table 2 (composition by embedding and alliance) | 66 | `triad_population_output.txt` |
| Figure 1 (raw onset rates by exposure, contiguity and distance) and Figure A2 (map of the most exposed non-contiguous pairs), Section 5.1 and Appendix J | 102 | `geography_descriptives_output.txt`, `fig_geography_descriptives.*`, `fig_geography_map.*` |
| Figure 2A (transitions, all dyad-years), Section 5.2, Appendix E | 86, 91, 92 | `population_selection_output.txt`, `fullpanel_shared_se_output.txt`, `fig_exposure_map_full.*` |
| Figure 2B (severity), Sections 5.2 and 5.5, Appendices D and E | 89, 92 | `fullpanel_severity_initiator_output.txt` |
| Section 5.3 (contiguity, distance), Appendix F | 88, 94 | `fullpanel_within_het_output.txt`, `fullpanel_contiguity_implied_output.txt` |
| Section 5.3, Appendix G (major powers, regions) | 84 | `distant_pairs_output.txt` |
| Figure 3 (predicted onsets), Appendix H | 85, 92 | `effect_sizes_output.txt`, `predicted_onset_by_percentile.csv`, `fig_predicted_onset.*` |
| Section 5.3 (colonial ties, territorial claims, issue in dispute, upper tail by era, coverage), Appendix J | 99 | `what_distinguishes_output.txt` |
| Section 5.4 (longer-run within-pair exposure), Appendix F | 101 | `longer_run_exposure_output.txt` |
| Section 5.3 (diplomatic exchange, UN voting alignment, halves of the period), Appendix J | 100 | `diplomacy_alignment_output.txt` |
| Section 5.4 (within-between, fixed effects), Appendix F | 88, 93, 98 | `fullpanel_within_het_output.txt`, `fullpanel_dyad_fe_output.txt`, `within_between_by_contiguity_output.txt` |
| Table 1 (controls) | 87 | `fullpanel_controls_output.txt` |
| Appendix E, Table 14 (predictors, held-out comparisons) | 90 | `fullpanel_predictors_output.txt` |
| Figure A1 (periods), Appendix E | 70, 79, 80 | `embedding_over_time_output.txt`, `era_bins_rates.csv`, `era_bins_shares.csv`, `fig_eras.*` |
| Appendix A, Section 5.2 (embedded pairs only) | 66, 68, 73, 77, 81 | `triad_population_output.txt`, `predictors_in_T_output.txt`, `claims_table_output.txt`, `hostility_gradient_output.txt`, `within_between_output.txt` |
| Appendix B (disputes inside unbalanced triads) | 71, 72 | `negative_tie_resolution_output.txt` |
| Appendix I, Sections 5.2 and 5.6 (lagged embedding and partner counts) | 86, 97 | `population_selection_output.txt`, `lagged_network_output.txt` |
| Appendix C, Section 5.6 (actor-overlap intervals) | 95 | `geography_actor_overlap_output.txt` |
| Appendix J (upper tail of exposure) | 96 | `upper_tail_output.txt` |

Scripts 01 to 65 build the data and are not needed for the article's tables and figures beyond the inputs listed in the steps above; their outputs are not stored.

## Files not included

`results/hazard_models.rds` (over 100 MB; regenerated by script 04). The four-treaty-type outputs of scripts 35 to 64 are not stored; unset `ALLIANCE_DEF` and re-run to obtain them.

## Data sources

Correlates of War: Alliance data v4.1 (dyadic yearly), Dyadic Militarized Interstate Disputes v4.03, Trade v4.0, Direct Contiguity v3.20, National Material Capabilities v7. Alliance Treaty Obligations and Provisions (ATOP) v5.1. Penn World Table 10.01 and the Maddison Project Database 2023 (GDP). Polity 5. CShapes 2.0 (capital distances, through the R package). Full references are in the article. Source files are unmodified in `raw/`.
