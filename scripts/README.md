# Scripts
Run from the repository root (`Rscript scripts/NN_name.R`). Set `ALLIANCE_DEF=defense` for scripts 35 onward (defence-pact definition of alliance); see the top-level README. `38_hazard_panel_lib.R` and `49_layered_hazards.R` are sourced by many later scripts.

| Range | Purpose | Where used |
|---|---|---|
| 00 | `00_unpack_data.sh`: unpack the gzipped intermediate files | setup |
| 01-05 | Build the dyad-year trade/capability file, signed network, triads, and first hazard models | data builders |
| 06-65 | Earlier triad-level and signed-network analyses; panel library in 38 | not used for the article's tables and figures |
| 48-49 | Four-state transitions; full dyad-year panel (preamble sourced by 66-102) | panel |
| 66 | Embedded-pair composition and transitions | Appendix A, Table 2 |
| 67-68 | Proximity diagnostics; predictors in embedded pairs | Appendix A |
| 69 | Dyad fixed effects, embedded pairs | Appendix A |
| 70 | Embedding across eras | Appendix E, Figure A1 |
| 71-72 | Disputes inside unbalanced triads: counts and models | Appendix B |
| 73-75 | Claims table and tables assembled from earlier output (embedded pairs) | Appendix A |
| 76 | Dyad-clustered versus actor-overlap standard errors | Appendix C |
| 77-78 | Severity of onsets and initiation, embedded pairs | Appendix A |
| 79-80 | Period data and figure | Appendix E, Figure A1 |
| 81, 83 | Within-between decomposition and escalation, embedded pairs | Appendix A |
| 82 | Interactions with contiguity, embedding, and distance | Section 5.3, Appendix F |
| 84 | Major powers and regions among non-contiguous pairs | Section 5.3, Appendix G |
| 85 | Effect sizes at percentiles of the larger share; influence checks | Figure 3, Appendix H |
| 86 | Full-panel transitions, severity, and lagged embedding | Figure 2, Appendices E and I |
| 87 | Controls path, full panel | Table 1 |
| 88 | Within-between, distance interaction, escalation, full panel | Sections 5.3-5.5, Appendix F |
| 89 | Severity and initiation, full panel | Figure 2B, Section 5.5, Appendix D |
| 90 | Predictors and held-out comparisons, full panel | Appendix E, Table 14 |
| 91 | Actor-overlap standard errors, full panel | Appendices C and E |
| 92 | Figures 2 and 3 | Figures 2, 3 |
| 93 | Dyad fixed-effects logit, full panel | Section 5.4, Appendix F |
| 94 | Contiguity interaction with implied effects, full panel | Section 5.3, Appendix F |
| 95 | Contiguity interaction with actor-overlap intervals | Appendix C |
| 96 | Upper tail of exposure among non-contiguous pairs; leave-out checks | Section 5.3, Appendix J |
| 97 | Lagged partner counts | Appendix I |
| 98 | Within-between by contiguity | Section 5.4, Appendix F |
| 99 | Colonial ties, territorial claims, issue in dispute, coverage | Section 5.3, Appendix J |
| 100 | Diplomatic exchange, UN voting alignment, halves of the period | Section 5.3, Appendix J |
| 101 | Longer-run exposure windows | Section 5.4, Appendix F |
| 102 | Geographic descriptives and map | Figure 1, Section 5.1; Figure A2, Appendix J |
`run_all_defense.sh` re-runs scripts 35 to 64 with `ALLIANCE_DEF=defense`.
