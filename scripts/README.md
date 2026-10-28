# Scripts
Run from the repository root (`Rscript scripts/NN_name.R`). Set `ALLIANCE_DEF=defense` (defence-pact definition of alliance); see the top-level README. `03_layered_hazards.R` is sourced by most later scripts.

| Script | Purpose | Where used |
|---|---|---|
| 01 | `01_unpack_data.sh`: unpack the gzipped intermediate files | setup |
| 02 | Build the dyad-year trade/capability file and the signed network from `raw/` (optional after unpacking) | data builder |
| 03 | Full dyad-year panel (preamble sourced by 04-40) and layered hazards | panel |
| 04 | Embedded-pair composition and transitions | Appendix A, Table 2 |
| 05-06 | Proximity diagnostics; predictors in embedded pairs | Appendix A |
| 07 | Dyad fixed effects, embedded pairs | Appendix A |
| 08 | Embedding across eras | Appendix D, Figure A1 |
| 09-10 | Disputes inside unbalanced triads: counts and models | not in the current appendix (removed 2026-10-27) |
| 11-13 | Claims table and tables assembled from earlier output (embedded pairs) | Appendix A |
| 14 | Dyad-clustered versus actor-overlap standard errors | Appendix B |
| 15-16 | Severity of onsets and initiation, embedded pairs | Appendix A |
| 17-18 | Period data and figure | Appendix D, Figure A1 |
| 19, 21 | Within-between decomposition and escalation, embedded pairs | Appendix A |
| 20 | Interactions with contiguity, embedding, and distance | Section 5.3, Appendix E |
| 22 | Major powers and regions among non-contiguous pairs | Section 5.3, Appendix F |
| 23 | Effect sizes at percentiles of the larger share; influence checks | Figure 3, Appendix G |
| 24 | Full-panel transitions, severity, and lagged embedding | Figure 2, Appendices E and I |
| 25 | Controls path, full panel | Table 1 |
| 26 | Within-between, distance interaction, escalation, full panel | Sections 5.3-5.5, Appendix E |
| 27 | Severity and initiation, full panel | Figure 2B, Section 5.5, Appendix C |
| 28 | Predictors and held-out comparisons, full panel | Appendix C, Table 13 |
| 29 | Actor-overlap standard errors, full panel | Appendices C and E |
| 30 | Figures 2 and 3 | Figures 2, 3 |
| 31 | Dyad fixed-effects logit, full panel | Section 5.4, Appendix E |
| 32 | Contiguity interaction with implied effects, full panel | Section 5.3, Appendix E |
| 33 | Contiguity interaction with actor-overlap intervals | Appendix B |
| 34 | Upper tail of exposure among non-contiguous pairs; leave-out checks | Section 5.3, Appendix H |
| 35 | Lagged partner counts | not in the current appendix (removed 2026-10-27) |
| 36 | Within-between by contiguity | Section 5.4, Appendix E |
| 37 | Colonial ties, territorial claims, issue in dispute, coverage | Section 5.3, Appendix H |
| 38 | Diplomatic exchange, UN voting alignment, halves of the period | Section 5.3, Appendix H |
| 39 | Longer-run exposure windows | Section 5.4, Appendix E |
| 40 | Geographic descriptives and map | Figure 1, Section 5.1; Figure A2, Appendix H |
| 41 | Geography interaction: calendar-year effects, top 1% and 5% dropped, spline | Section 5.3, Appendix G (influence) |
| 42 | Geography interaction under alternative minimum-trade rules (none, 0.01% and 0.05% of world trade) | Section 4.3, Appendix H |
`run_all.sh` runs scripts 03 to 42 in order with `ALLIANCE_DEF=defense`.
