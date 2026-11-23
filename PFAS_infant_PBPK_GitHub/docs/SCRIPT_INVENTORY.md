# Script inventory and provenance

| Repository file | Source upload | Repository treatment |
|---|---|---|
| `R/00_run_pipeline.R` | `00_run_pipeline.R` | Rewritten to run the manuscript-aligned no-water pipeline |
| `R/01_run_main_model.R` | `01_run_main_model(1).R` | Copied; filename normalized only |
| `R/02_four_layer_evaluation.R` | `02_four_layer_evaluation.R` | Copied |
| `R/06_trajectory_proxy_sensitivity.R` | `06_trajectory_proxy_sensitivity.R` | Copied |
| `R/07_Varsi_study_specific_revalidation.R` | `07_Varsi_study_specific_revalidation.R` | Copied |
| `R/08_integrated_sensitivity_analysis.R` | `08_integrated_sensitivity_analysis.R` | Copied |
| `R/10_global_sensitivity_uncertainty_FAST175.R` | `10_global_sensitivity_uncertainty.R` | Copied; filename normalized |
| `R/10b_FAST175_posthoc_convergence_noODE.R` | `10_FAST175_posthoc_convergence_noODE.R` | Copied; filename normalized |
| `R/15_infant_PFAS_risk_characterization_NO_WATER.R` | New, derived from current manuscript methods and uploaded PBPK output schema | Reconstructed; no water input |
| `R/16_export_main_figures_1_2_NO_WATER.R` | New, based on current Figure 1–2 definitions and PBPK output schema | Reconstructed; no water input |
| `R/17_export_IBQ_organ_contribution_figures_NO_WATER.R` | `17_export_IBQ_organ_contribution_figures_NO_WATER.R` | Copied; input/output paths aligned to repository |
| `R/18_export_Supplementary_Figure_S1_NO_WATER.R` | New, based on current Supplementary Figure S1 description and available output fields | Reconstructed to calculate ranks from outputs |

The reconstructed scripts are candidate replacements and require execution plus table/figure reconciliation before publication. They have not been run against the private cohort workbook in this environment.
