#!/usr/bin/env Rscript
# Run from repository root. Private input data must be supplied locally.
options(stringsAsFactors = FALSE)
ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", getwd()), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(ROOT, "R", "01_run_main_model.R")))
  stop("Set working directory/PFAS_PBPK_ROOT to repository root.")
input <- Sys.getenv("PFAS_PBPK_INPUT", unset = file.path(ROOT, "data", "input-PBPK.xlsx"))
if (!file.exists(input)) stop("Missing authorized cohort input: ", input,
  "\nDo not commit participant-level data to the public repository.")
exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run <- function(script) {
  f <- file.path(ROOT, "R", script)
  if (!file.exists(f)) stop("Missing script: ", f)
  message("\n========== RUN ", script, " ==========")
  status <- system2(exe, shQuote(f), stdout = "", stderr = "")
  if (status != 0) stop("Script failed: ", script, " (status ", status, ")")
}
# Primary manuscript pipeline
run("01_run_main_model.R")
run("15_infant_PFAS_risk_characterization_NO_WATER.R")
run("16_export_main_figures_1_2_NO_WATER.R")
run("17_export_IBQ_organ_contribution_figures_NO_WATER.R")
run("18_export_Supplementary_Figure_S1_NO_WATER.R")
# Optional manuscript-supporting analyses; activate explicitly.
if (Sys.getenv("PFAS_RUN_FOUR_LAYER", "0") == "1") run("02_four_layer_evaluation.R")
if (Sys.getenv("PFAS_RUN_VARSI", "0") == "1") run("07_Varsi_study_specific_revalidation.R")
if (Sys.getenv("PFAS_RUN_TRAJECTORY_SENSITIVITY", "0") == "1") run("06_trajectory_proxy_sensitivity.R")
if (Sys.getenv("PFAS_RUN_INTEGRATED_SENSITIVITY", "0") == "1") run("08_integrated_sensitivity_analysis.R")
if (Sys.getenv("PFAS_RUN_FAST175", "0") == "1") run("10_global_sensitivity_uncertainty_FAST175.R")
if (Sys.getenv("PFAS_RUN_FAST175_POSTHOC", "0") == "1") run("10b_FAST175_posthoc_convergence_noODE.R")
message("\nPrimary NO-WATER manuscript pipeline completed. Review every QC sheet and version note.")
