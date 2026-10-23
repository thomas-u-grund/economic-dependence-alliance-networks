# Codebook

## 1. Source files (`raw/`, unmodified)
| Folder | Content |
|---|---|
| `alliances/version4.1_csv/` | COW formal alliances v4.1, dyad-year file used here (`alliance_v4.1_by_dyad_yearly.csv`) and member file (`..._by_member.csv`). `defense == 1` marks a defence pact (Type I). |
| `mids/dyadic_mid_4.03_update/` | COW dyadic Militarized Interstate Disputes v4.03 (`dyadic_mid_4.03.csv`) and codebook. Used: `hihost` (highest hostility level, 2 threat, 3 display, 4 use of force, 5 war), `fatlev`, `dyad_rolea/b` (1 = initiator of the dyadic dispute), `sideaa/sideab`. |
| `trade/COW_Trade_4.0/` | COW bilateral trade v4.0 (current US$ millions). |
| `contiguity/DirectContiguity320/` | COW direct contiguity v3.20 (`contdird.csv`; `conttype` 1 = land border, 2-5 = water). |
| `nmc/NMCv7/` | National Material Capabilities v7 (`cinc`). |
| `atop_download/`, `polity/`, `gdp_download/`, `states/` | ATOP 5.1; Polity 5; GDP downloads (Penn World Table, Maddison); COW state system. |

## 2. Built data files (`data/`)
| File | Content |
|---|---|
| `tie_trade.csv(.gz)` | One row per dyad-year (1870-2011) with trade and capability measures: `ccode_low`, `ccode_high` (COW codes, low < high), `trade_share_low`, `trade_share_high` (share of each state's own total trade that is with the partner; missing if total trade under US$50 million or if the dyadic flow exceeds either state's total), `total_dyadic_trade` (US$ millions), `log_trade_volume`, `cinc_low`, `cinc_high`, `cap_asymmetry` (absolute difference in log CINC), `dep_exposed`, `dep_sheltered`, `dep_asymmetry` (larger share, smaller share, their difference), and further transformations. |
| `signed_ties.csv` | Signed network of the original triad design: `year`, `ccode_low`, `ccode_high`, `sign` (+1 alliance, -1 dispute; dispute overrides alliance). Built with all four commitment types. The `multiplex` column is not used (it is always FALSE). The defence-pact network is built in memory by `scripts/03_layered_hazards.R`. |
| `state_years.csv`, `gdp_maddison.csv` | State-years (COW code, year); Maddison GDP by COW code. |
| `triads_all_years.csv(.gz)`, `.rds`, `triads_summary_by_year.csv` | Triads of the signed network by year (original design). |
| `person_year_risk_table.csv/.rds`, `spell_bounds.rds` | Triad-spell risk table of the triad-level design (Appendices A-D, G). |

## 3. Analysis panel (built in memory by scripts 03 to 40)
One row per dyad-year with trade data, 1870-2010 (611,687 rows). Defence-pact definition (`ALLIANCE_DEF=defense`).
| Variable | Definition |
|---|---|
| `year` | Calendar year t |
| `lo`, `hi`, `dyad` | COW codes of the two states (lo < hi); dyad identifier "lo hi" (clustering unit) |
| `a` | 1 if one defence pact covers both states in year t (all four commitment types if `ALLIANCE_DEF` is unset) |
| `m` | 1 if the two states are in a militarised dispute of hostility 3 or higher in year t |
| `a1`, `m1` | Same, in year t+1 (the panel keeps dyad-years observed in t+1) |
| `contig` | 1 if directly contiguous by land or water in year t |
| `shl`, `shh` | Trade share of the lower-/higher-code state (see `tie_trade`) |
| `mx`, `mn` | Larger and smaller of `shl`, `shh` (the exposure of the more and the less exposed state) |
| `cap` | Capability asymmetry |
| `vol`, `lvol` | Total dyadic trade (US$ millions); log(1 + vol) |
| `dur_a`, `dur_m` | Consecutive observed years the dyad has been in its current alliance / dispute state |
| `n_same`, `n_diff` | Number of third states whose ties (alliance or dispute) to the two states have the same sign (common partners) or opposite signs (split partners) in the signed network of year t; `l_same`, `l_diff` = log(1 + n) |
| `tri` / `T` | 1 if `n_same + n_diff >= 1` (embedded pair) |
| `z_mx`, `z_mn`, `z_cap`, `z_vol` | `mx`, `mn`, `cap`, `lvol` standardised over the whole panel (mean 0, SD 1) |
| `onset_m`, `end_m`, `form_a`, `end_a` | Transition indicators: dispute begins (`m == 0` and `m1 == 1`), dispute ends (`m == 1`, `m1 == 0`), alliance forms (`a == 0`, `a1 == 1`), alliance ends (`a == 1`, `a1 == 0`) |
| `decade` | Decade factor |
Scripts 10 and 15-16 add: `n_pp` (common allies), `n_nn` (common enemies), `capdist` (capital distance, CShapes, km), `pact` (largest covering defence pact of five or more members, otherwise "other"/"none"), `h`/`fat` (highest hostility and fatality level in the dyad-year), `init` (initiating state).

## 4. Output files (`results/`)
`*_output.txt` are the printed output of the script that produces them (see the README) (captured with `sink`); `*.csv` are tables behind figures and appendix tables; `*.png/.pdf` are figures; `*.rds` are saved model objects. See README for the map from article items to files.
