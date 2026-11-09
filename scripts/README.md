# Scripts
Run from the repository root (`Rscript scripts/NN_name.R`). Set `ALLIANCE_DEF=defense` (defence-pact definition of alliance); see the top-level README. `03_layered_hazards.R` is sourced by most later scripts.

| Script | Purpose | Where used |
|---|---|---|
| 01 | `01_unpack_data.sh`: unpack the gzipped analysis panel | setup |
| 02 | Build the dyad-year trade/capability file and the signed network from `raw/` (optional after unpacking) | data builder |
| 03 | Full dyad-year panel (preamble sourced by 04-42) and layered hazards | panel |
| 04 | Embedded-pair composition and transitions; pooled and fully interacted embedding comparisons | Appendix A (embedded pairs) |
| 05-07 | Proximity diagnostics, predictors, and dyad fixed effects among embedded pairs | further analysis, not in the appendix |
| 08 | Embedding across eras | further analysis, not in the appendix |
| 10-13 | Claims table and tables assembled from earlier output (embedded pairs) | further analysis, not in the appendix |
| 14 | Dyad-clustered versus actor-overlap standard errors, embedded pairs | further analysis, not in the appendix |
| 15-16 | Severity of onsets and initiation, embedded pairs | further analysis, not in the appendix |
| 17-18 | Period data and figure (embedding over time) | further analysis, not in the appendix |
| 19 | Within-between decomposition, embedded pairs | Section 5.4, Appendix A (embedded pairs) |
| 20 | Interactions with contiguity, embedding, and distance | Section 5.3, Appendix A (embedded pairs), Appendix B Table 4 |
| 21 | Escalation given onset, embedded pairs | further analysis, not in the appendix |
| 22 | Major powers and regions among non-contiguous pairs | Section 5.3, Appendix B Table 5 |
| 23 | Effect sizes at percentiles of the larger share; influence checks | Figure 3, Appendix B Tables 6 and 7 |
| 24 | Full-panel transitions | Figure 2A, Appendix A Table 1 |
| 25 | Controls path, full panel | Table 3 |
| 26 | Within-between, distance interaction, escalation given onset, full panel | Sections 5.3-5.5, Appendix A Table 3, Appendix B Table 4, Appendix D Table 14 |
| 27 | Severity and initiation, full panel | Figure 2B, Section 5.5, Appendix D Tables 13 and 15 |
| 28 | Predictors and held-out comparisons, full panel | further analysis, not in the appendix |
| 29 | Actor-overlap standard errors, full panel | Section 5.6, Appendix A Table 2 |
| 30 | Figures 2 and 3 | Figures 2, 3 |
| 31 | Dyad fixed-effects logit, full panel | Section 5.4, Appendix A |
| 32 | Contiguity interaction with implied effects, full panel | Table 1, Appendix B Table 4 |
| 33 | Contiguity interaction with actor-overlap intervals | Table 1, Appendix A |
| 34 | Upper tail of exposure among non-contiguous pairs; leave-out checks | Section 5.3, Appendix C Tables 8 and 9 |
| 36 | Within-between by contiguity | Table 2, Section 5.4, Appendix A |
| 37 | Colonial ties, territorial claims, issue in dispute, coverage | Sections 5.3 and 5.6, Appendix C Tables 10 and 11 |
| 38 | Diplomatic exchange, UN voting alignment, halves of the period | Section 5.3, Appendix C Table 12 |
| 39 | Longer-run exposure windows | Section 5.4, Appendix A |
| 40 | Geographic descriptives and map | Figure 1, Section 5.1; Figure A1, Appendix C |
| 41 | Geography interaction: calendar-year effects, top 1% and 5% dropped, spline | Section 5.3, Appendix B |
| 42 | Geography interaction under alternative minimum-trade rules (none, 0.01% and 0.05% of world trade) | Section 4.3, Appendix B |
`run_all.sh` runs scripts 03 to 42 (numbers 09 and 35 are unused) in order with `ALLIANCE_DEF=defense`.
