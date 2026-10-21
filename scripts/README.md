# Scripts
Run from the repository root (`Rscript scripts/NN_name.R`). Set `ALLIANCE_DEF=defense` for scripts 35 onward (defence-pact definition of alliance); see the top-level README. `38_hazard_panel_lib.R` and `49_layered_hazards.R` are sourced by many later scripts.

| Range | Purpose | Where used |
|---|---|---|
| 00 | `00_unpack_data.sh`: unpack the gzipped intermediate files | setup |
| 01-05 | Build the dyad-year trade/capability file, signed network, triads, person-year panel, first hazard models | data builders; Appendices A-D, G |
| 06-34 | Robustness checks, clustering, figures, triad-level design | Appendices A-D, G |
| 35-47, 50-62 | Dyad-year signed-network hazard design; panel library in 38 | Appendices E, F, H, I, J, K |
| 48-49 | Four-state transitions; full dyad-year panel (preamble sourced by 66-80) | panel |
| 63-65 | Proximity and interaction among allied dyads, allied dispute cases, descriptives | Appendix K |
| 66 | Embedded-pair composition and transitions | Table 1, Appendix L |
| 67 | Proximity and interaction diagnostics in embedded pairs | Table 2 |
| 68 | Predictors; exposure versus structure | Table 3 |
| 69 | Dyad fixed effects | Sections 5.3, 5.9 |
| 70 | Embedding and balance across eras | Section 5.5, Appendix P |
| 71-72 | Disputes inside unbalanced triads: counts and models | Section 5.6, Appendix M |
| 73 | Claims table (estimates, intervals, minimum detectable effects) | Figure 2A, Appendix P |
| 74-75 | Tables assembled from earlier output | Tables 2, 3, Appendix P |
| 76 | Dyad-clustered versus actor-overlap standard errors | Appendix N |
| 77-78 | Severity of onsets and initiation | Figure 2B, Section 5.7, Appendix O |
| 79-80 | Period data and Figures 2-3 | Figures 2, 3 |
`run_all_defense.sh` re-runs scripts 35 to 64 with `ALLIANCE_DEF=defense`.
