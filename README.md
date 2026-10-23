# Replication materials

*Economic Dependence and Interstate Conflict: Geography and Differences Between Relationships*

Data, commented scripts, codebook, and software output for every analysis in the article and its Online Appendix.

## Contents

```
scripts/        01-40, run in numerical order (index: scripts/README.md)
raw/            source data as downloaded (Correlates of War alliances, MIDs, trade, contiguity, capabilities;
                ATOP 5.1; Penn World Table / Maddison GDP; Polity 5; state system)
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
2. The article's alliance definition is the Correlates of War defence pact (Type I). Set `ALLIANCE_DEF=defense` for all scripts (`export ALLIANCE_DEF=defense`); `scripts/run_all.sh` does this and runs scripts 03 to 40 in order. Leave it unset to use all four commitment types.
3. Scripts 03 to 40 build the dyad-year analysis panel in memory from `raw/` and `data/tie_trade.csv` (the shared preamble is in `scripts/03_layered_hazards.R`) and write their output to `results/`. Scripts 12, 13, and 18 only assemble earlier output and need scripts 04 to 11 and 15 to 17 to have run.

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
| Figure 1 (raw onset rates by exposure, contiguity and distance) and Figure A2 (map of the most exposed non-contiguous pairs), Section 5.1 and Appendix J | 40 | `geography_descriptives_output.txt`, `fig_geography_descriptives.*`, `fig_geography_map.*` |
| Figure 2A (transitions, all dyad-years), Section 5.2, Appendix E | 24, 29, 30 | `population_selection_output.txt`, `fullpanel_shared_se_output.txt`, `fig_exposure_map_full.*` |
| Figure 2B (severity), Sections 5.2 and 5.5, Appendices D and E | 27, 30 | `fullpanel_severity_initiator_output.txt` |
| Section 5.3 (contiguity, distance), Appendix F | 26, 32 | `fullpanel_within_het_output.txt`, `fullpanel_contiguity_implied_output.txt` |
| Section 5.3, Appendix G (major powers, regions) | 22 | `distant_pairs_output.txt` |
| Figure 3 (predicted onsets), Appendix H | 23, 30 | `effect_sizes_output.txt`, `predicted_onset_by_percentile.csv`, `fig_predicted_onset.*` |
| Section 5.3 (colonial ties, territorial claims, issue in dispute, upper tail by era, coverage), Appendix J | 37 | `what_distinguishes_output.txt` |
| Section 5.4 (longer-run within-pair exposure), Appendix F | 39 | `longer_run_exposure_output.txt` |
| Section 5.3 (diplomatic exchange, UN voting alignment, halves of the period), Appendix J | 38 | `diplomacy_alignment_output.txt` |
| Section 5.4 (within-between, fixed effects), Appendix F | 26, 31, 36 | `fullpanel_within_het_output.txt`, `fullpanel_dyad_fe_output.txt`, `within_between_by_contiguity_output.txt` |
| Table 1 (controls) | 25 | `fullpanel_controls_output.txt` |
| Appendix E, Table 14 (predictors, held-out comparisons) | 28 | `fullpanel_predictors_output.txt` |
| Figure A1 (periods), Appendix E | 08, 17, 18 | `embedding_over_time_output.txt`, `era_bins_rates.csv`, `era_bins_shares.csv`, `fig_eras.*` |
| Appendix A, Section 5.2 (embedded pairs only) | 04, 06, 11, 15, 19 | `triad_population_output.txt`, `predictors_in_T_output.txt`, `claims_table_output.txt`, `hostility_gradient_output.txt`, `within_between_output.txt` |
| Appendix B (disputes inside unbalanced triads) | 09, 10 | `negative_tie_resolution_output.txt` |
| Appendix I, Sections 5.2 and 5.6 (lagged embedding and partner counts) | 24, 35 | `population_selection_output.txt`, `lagged_network_output.txt` |
| Appendix C, Section 5.6 (actor-overlap intervals) | 33 | `geography_actor_overlap_output.txt` |
| Appendix J (upper tail of exposure) | 34 | `upper_tail_output.txt` |

## Files not included

`results/hazard_models.rds` (over 100 MB; not needed for the article).

## Data sources

Correlates of War: Alliance data v4.1 (dyadic yearly), Dyadic Militarized Interstate Disputes v4.03, Trade v4.0, Direct Contiguity v3.20, National Material Capabilities v7. Alliance Treaty Obligations and Provisions (ATOP) v5.1. Penn World Table 10.01 and the Maddison Project Database 2023 (GDP). Polity 5. CShapes 2.0 (capital distances, through the R package). Full references are in the article. Source files are unmodified in `raw/`.
