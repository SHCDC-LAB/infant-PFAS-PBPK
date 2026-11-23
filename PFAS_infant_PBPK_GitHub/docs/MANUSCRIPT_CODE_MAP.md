# Mapping between the current manuscript and R scripts

This mapping is based on the supplied current main manuscript and Supplementary Material.

| Manuscript component | Code | Notes |
|---|---|---|
| S1 six-compound infant PBPK simulation | `01_run_main_model.R` | Main model; requires authorized cohort workbook |
| Model qualification / external evidence | `02_four_layer_evaluation.R`, `07_Varsi_study_specific_revalidation.R` | Supporting evaluation modules; confirm which output tables feed current Supplementary Tables S5–S6 |
| Milk trajectory proxy sensitivity | `06_trajectory_proxy_sensitivity.R` | Optional robustness analysis for PFNA, PFDA and PFUdA proxy trajectory |
| Integrated sensitivity and numerical/renal diagnostics | `08_integrated_sensitivity_analysis.R` | Optional supporting analysis; confirm exact output-to-table mapping |
| Global sensitivity (FAST175) | `10_global_sensitivity_uncertainty_FAST175.R` | Supporting sensitivity analysis; manuscript Supplementary Table S7 reports LHS/PRCC results |
| FAST175 convergence / rank robustness | `10b_FAST175_posthoc_convergence_noODE.R` | Post-hoc CSV-only diagnostic; not part of the primary model run |
| Milk-derived EDI and four-PFAS internal benchmark | `15_infant_PFAS_risk_characterization_NO_WATER.R` | Reconstructs EDI from model-recorded ingestion rate and computes manuscript-specific serum-equivalent ratios |
| Main Figures 1–2 | `16_export_main_figures_1_2_NO_WATER.R` | Figure 1 EDI distributions; Figure 2 plasma/liver/kidney trajectories |
| Main Figures 3–4 | `17_export_IBQ_organ_contribution_figures_NO_WATER.R` | Figure 3 within-organ six-PFAS CTWA composition; Figure 4 four-PFAS internal benchmark quotient |
| Supplementary Figure S1 | `18_export_Supplementary_Figure_S1_NO_WATER.R` | Derives rankings from generated EDI and PBPK outputs rather than hard-coded ranks |

## Deliberately excluded from this manuscript repository

- `03_source_attribution_sensitivity.R`: source-attribution analysis is not reported in the supplied current manuscript.
- `04_export_EP_main_figures.R`: superseded figure export script.
- `05_human_organ_benchmark.R`: no corresponding independent infant organ validation result is presented in the current manuscript; adult tissue literature is discussed as context only.
- `09_water_PFAS_environmental_linkage.R` and `09b_water_ND_sensitivity.R`: drinking-water analyses excluded from the current manuscript.
- `11_water_milk_21PFAS_regional_analysis.R` and `12_water_milk_PBPK_integrated_analysis.R`: regional water–milk–PBPK analyses excluded.
- `13_cross_matrix_PFAS_prioritization_noODE.R`: cross-matrix prioritization depends on the excluded water analyses.
- `14_PBPK_incremental_information_noODE.R`: incremental-information diagnostic is not reported in the current manuscript.
- `15_infant_PFAS_risk_characterization.R` and `15_infant_PFAS_risk_characterization_v2.R`: earlier water-containing risk-characterization versions; replaced in this repository by the manuscript-aligned NO-WATER script.
- `16_export_infant_PFAS_risk_figures_journal.R`: older exporter includes a drinking-water figure and an obsolete internal benchmark scenario; replaced by the NO-WATER Figure 1–2 exporter.

## Manuscript-to-output verification still required

The scripts generate candidate outputs. Before release, compare the output workbook values with every corresponding main-text and supplementary table, and compare the exported figure panels/labels with the submitted manuscript. A successful script run alone is not proof that the output reproduces the manuscript.
