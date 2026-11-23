# ============================================================
# 06_trajectory_proxy_sensitivity.R -- COMPLETE S1 beta ensemble
# PBPK infant PFAS / total-plasma kidney exchange / v2
# ============================================================
# PURPOSE
#   Compare four breast-milk trajectory assumptions for PFNA,
#   PFDA and PFUdA, using the ACTUAL S1 PBPK simulation outputs.
#   PFOA, PFHxS and PFOS trajectories remain unchanged.
#   Computes PLASMA, LIVER, KIDNEY and BODY BURDEN endpoints,
#   selected-age summaries, mode-vs-proxy changes and paired ID QC.
#
# REQUIRED FILES
#   R/01_run_main_model.R
#   data/input-PBPK.xlsx (only needed if additional modes run)
#   Results are written by 01 to:
#   results/PFAS_infant_PBPK_S1_renal_totalplasma_v2[/_beta_MODE]/tables/
#     PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx
#     PFAS_infant_PBPK_final_freeze_v1_summary.xlsx
#
# DEFAULT BEHAVIOR: auto-run ONLY missing modes, one at a time.
#   Sys.setenv(PFAS_RUN_BETA_ENSEMBLE="1")  # run missing modes (default)
#   Sys.setenv(PFAS_RUN_BETA_ENSEMBLE="0")  # audit only; FAIL if missing
#   Sys.setenv(PFAS_BETA_FORCE_RERUN="1")    # rerun all 4, even if present
#   Sys.setenv(PFAS_PBPK_ROOT="...")        # project root; else getwd()
#
# CAUTION: 3 or 4 full PBPK simulations can be computationally costly.
# This code deliberately never treats absent alternatives as completed.
# ============================================================

options(stringsAsFactors = FALSE)

ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", unset = getwd()),
                      winslash = "/", mustWork = TRUE)
MAIN <- file.path(ROOT, "R", "01_run_main_model.R")
if (!file.exists(MAIN)) stop("Main model script not found: ", MAIN,
                             "\nSet PFAS_PBPK_ROOT to the GitHub repository root.")

if (!requireNamespace("readxl", quietly = TRUE)) stop("Please install.packages('readxl')")
if (!requireNamespace("writexl", quietly = TRUE)) stop("Please install.packages('writexl')")

BASE_VERSION <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
MODES <- c("proxy", "flat", "half", "strong")
TARGET_DAYS <- c(0, 30, 90, 180, 270, 365.25)
COMPOUNDS <- c("PFOA", "PFNA", "PFDA", "PFUdA", "PFHxS", "PFOS")
PROXY_COMPOUNDS <- c("PFNA", "PFDA", "PFUdA")
EXPECTED_N <- 175L
RUN_MISSING <- Sys.getenv("PFAS_RUN_BETA_ENSEMBLE", unset = "1") == "1"
FORCE_RERUN <- Sys.getenv("PFAS_BETA_FORCE_RERUN", unset = "0") == "1"

out_dir <- file.path(ROOT, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_dir <- file.path(out_dir, "beta_ensemble_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

mode_version <- function(mode) {
  if (identical(mode, "proxy")) BASE_VERSION else paste0(BASE_VERSION, "_beta_", mode)
}
mode_tables <- function(mode) file.path(out_dir, mode_version(mode), "tables")
mode_paths <- function(mode) list(
  full = file.path(mode_tables(mode), "PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx"),
  summary = file.path(mode_tables(mode), "PFAS_infant_PBPK_final_freeze_v1_summary.xlsx")
)
mode_complete <- function(mode) {
  p <- mode_paths(mode)
  all(file.exists(unlist(p))) && all(file.info(unlist(p))$size > 1000)
}

status <- data.frame(
  Trajectory_Mode = MODES,
  Model_version = vapply(MODES, mode_version, character(1)),
  Full_Results_Present_Before = vapply(MODES, function(m) file.exists(mode_paths(m)$full), logical(1)),
  Summary_Present_Before = vapply(MODES, function(m) file.exists(mode_paths(m)$summary), logical(1)),
  Action = "reused",
  Completed = FALSE,
  Full_Results_Path = vapply(MODES, function(m) mode_paths(m)$full, character(1)),
  stringsAsFactors = FALSE
)

message("\n===== S1 milk beta ensemble =====")
message("Project root: ", ROOT)
message("Run missing scenarios: ", RUN_MISSING)
message("Force rerun: ", FORCE_RERUN)
for (m in MODES) message(sprintf("  %-7s : %s", m, if (mode_complete(m)) "both files found" else "MISSING file(s)"))

# Execute each mode in an INDEPENDENT R process because 01 is a full script
# and its constants, parallel plan and outputs must not leak across runs.
run_one_mode <- function(mode) {
  exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (!file.exists(exe)) stop("Rscript executable not found: ", exe)
  stdout_log <- file.path(log_dir, paste0("beta_", mode, ".log"))
  stderr_log <- file.path(log_dir, paste0("beta_", mode, "_errors.log"))
  old_root <- Sys.getenv("PFAS_PBPK_ROOT", unset = NA_character_)
  old_mode <- Sys.getenv("PFAS_TRAJECTORY_MODE", unset = NA_character_)
  on.exit({
    if (is.na(old_root)) Sys.unsetenv("PFAS_PBPK_ROOT") else Sys.setenv(PFAS_PBPK_ROOT = old_root)
    if (is.na(old_mode)) Sys.unsetenv("PFAS_TRAJECTORY_MODE") else Sys.setenv(PFAS_TRAJECTORY_MODE = old_mode)
  }, add = TRUE)
  Sys.setenv(PFAS_PBPK_ROOT = ROOT, PFAS_TRAJECTORY_MODE = mode)
  message("\nStarting full PBPK run: ", mode)
  message("Log: ", stdout_log)
  exit_code <- tryCatch(
    system2(exe, args = shQuote(MAIN), stdout = stdout_log, stderr = stderr_log, wait = TRUE),
    error = function(e) stop("Failed starting mode '", mode, "': ", conditionMessage(e))
  )
  if (length(exit_code) != 1L || is.na(exit_code) || exit_code != 0) {
    last_lines <- if (file.exists(stderr_log)) tail(readLines(stderr_log, warn = FALSE), 18) else character()
    stop("PBPK simulation failed for mode '", mode, "' (exit code ", exit_code, ").\n",
         "See logs: ", stdout_log, " and ", stderr_log, "\n",
         paste(last_lines, collapse = "\n"))
  }
  if (!mode_complete(mode)) {
    stop("R process exited successfully but one or both expected output files are missing for '",
         mode, "'. Check: ", mode_tables(mode), "\nConsult: ", stdout_log)
  }
  message("Completed: ", mode)
}

for (i in seq_along(MODES)) {
  mode <- MODES[i]
  needs_run <- FORCE_RERUN || !mode_complete(mode)
  if (needs_run && RUN_MISSING) {
    run_one_mode(mode)
    status$Action[i] <- if (FORCE_RERUN) "recalculated_forced" else "calculated_missing"
  } else if (needs_run) {
    status$Action[i] <- "missing_not_run"
  }
  status$Completed[i] <- mode_complete(mode) && status$Action[i] != "missing_not_run"
}

status_path <- file.path(out_dir, "beta_proxy_run_status.csv")
utils::write.csv(status, status_path, row.names = FALSE, fileEncoding = "UTF-8")
if (!all(status$Completed)) {
  missing <- status$Trajectory_Mode[!status$Completed]
  stop("INCOMPLETE BETA ENSEMBLE: ", paste(missing, collapse = ", "),
       ". No misleading summary has been written.\n",
       "Status CSV: ", status_path, "\n",
       "To calculate missing modes: Sys.setenv(PFAS_RUN_BETA_ENSEMBLE='1'); ",
       "source('R/06_trajectory_proxy_sensitivity.R')")
}

# Read ONLY needed columns from the large (>50K row) PBPK workbook.
# This dramatically reduces memory and enables all four mode analyses.
needed <- c("ID", "Scenario", "Compound", "Time_days", "C_Plasma",
            "C_Liver", "C_Kidney", "A_Body_Burden")
read_full_select <- function(path) {
  headers <- names(readxl::read_excel(path, sheet = "PBPK_Results", n_max = 0,
                                      .name_repair = "minimal"))
  missing <- setdiff(needed, headers)
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  if (anyDuplicated(headers)) stop("Duplicate headers in ", path)
  col_types <- ifelse(headers %in% needed, "guess", "skip")
  d <- as.data.frame(readxl::read_excel(
    path, sheet = "PBPK_Results", col_types = col_types,
    .name_repair = "minimal", progress = FALSE
  ), stringsAsFactors = FALSE)
  d <- d[, needed, drop = FALSE]
  d$Compound <- as.character(d$Compound)
  d$Scenario <- as.character(d$Scenario)
  d$ID <- as.character(d$ID)
  for (nm in c("Time_days", "C_Plasma", "C_Liver", "C_Kidney", "A_Body_Burden")) {
    d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
  }
  d
}

# ID- and compound-specific exact last time is selected; no arbitrary 364-day
# substitution for the model's explicit 365.25-day endpoint.
endpoint_rows <- function(d, mode) {
  d <- d[!is.na(d$Time_days) & d$Scenario == "S1", , drop = FALSE]
  if (nrow(d) == 0) stop("No S1 records found: ", mode)
  if (!setequal(unique(d$Compound), COMPOUNDS)) {
    stop("Unexpected compounds in mode ", mode, ": ", paste(unique(d$Compound), collapse = ", "))
  }
  endpoint <- d[abs(d$Time_days - 365.25) < 0.001, , drop = FALSE]
  if (!nrow(endpoint)) stop("365.25-day endpoint not present for mode ", mode)
  k <- paste(endpoint$ID, endpoint$Compound, sep = "||")
  if (anyDuplicated(k)) stop("Duplicated ID/Compound endpoint in mode ", mode)
  counts <- table(factor(endpoint$Compound, levels = COMPOUNDS))
  if (any(counts != EXPECTED_N)) {
    stop("Unexpected endpoint subject counts for mode ", mode, ": ",
         paste(names(counts), counts, collapse = "; "))
  }
  endpoint$Trajectory_Mode <- mode
  endpoint$Model_version <- mode_version(mode)
  endpoint
}

summarize_dist <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(Median = NA_real_, P025 = NA_real_, P975 = NA_real_))
  c(Median = stats::median(x),
    P025 = as.numeric(stats::quantile(x, .025, names = FALSE)),
    P975 = as.numeric(stats::quantile(x, .975, names = FALSE)))
}

measurements <- c(
  Plasma = "C_Plasma", Liver = "C_Liver",
  Kidney = "C_Kidney", Body_Burden = "A_Body_Burden"
)
units <- c(Plasma = "ng/L", Liver = "ng/L", Kidney = "ng/L", Body_Burden = "ng")

summarize_rows <- function(d, mode, group_label = "Endpoint_12m") {
  result <- vector("list", length(COMPOUNDS) * length(measurements))
  z <- 0L
  for (cmp in COMPOUNDS) {
    sub <- d[d$Compound == cmp, , drop = FALSE]
    for (meas in names(measurements)) {
      z <- z + 1L
      q <- summarize_dist(sub[[measurements[[meas]]]])
      result[[z]] <- data.frame(
        Trajectory_Mode = mode, Model_version = mode_version(mode),
        Compound = cmp, Variable = meas, Units = units[[meas]],
        N_ID = length(unique(sub$ID)),
        Median = unname(q["Median"]), P025 = unname(q["P025"]),
        P975 = unname(q["P975"]), stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, result)
}

all_endpoint <- vector("list", length(MODES))
all_summary <- vector("list", length(MODES))
selected_summary <- list()
paired_endpoints <- list()
qc <- vector("list", length(MODES))

for (j in seq_along(MODES)) {
  m <- MODES[j]
  message("Reading simulated plasma/liver/kidney: ", m)
  d <- read_full_select(mode_paths(m)$full)
  e <- endpoint_rows(d, m)
  qc[[j]] <- data.frame(
    Trajectory_Mode = m, N_Rows = nrow(d),
    N_Endpoint = nrow(e), N_ID = length(unique(e$ID)),
    N_Compounds = length(unique(e$Compound)),
    Max_Reported_Day = max(d$Time_days, na.rm = TRUE),
    Negative_Plasma = sum(d$C_Plasma < -1e-9, na.rm = TRUE),
    Negative_Liver = sum(d$C_Liver < -1e-9, na.rm = TRUE),
    Negative_Kidney = sum(d$C_Kidney < -1e-9, na.rm = TRUE),
    Nonfinite_Endpoint = sum(!is.finite(as.matrix(e[,unname(measurements),drop=FALSE]))),
    stringsAsFactors = FALSE
  )
  all_endpoint[[j]] <- e
  all_summary[[j]] <- summarize_rows(e, m)
  # Each requested age is matched to the closest ACTUAL model output day.
  observed_days <- sort(unique(d$Time_days[is.finite(d$Time_days)]))
  for (target in TARGET_DAYS) {
    actual <- observed_days[which.min(abs(observed_days - target))]
    sub <- d[abs(d$Time_days - actual) < 1e-5 & d$Scenario == "S1", , drop = FALSE]
    ss <- summarize_rows(sub, m)
    ss$Requested_Day <- target
    ss$Actual_Model_Day <- actual
    selected_summary[[length(selected_summary) + 1L]] <- ss
  }
  rm(d, e)
  gc(verbose = FALSE)
}

end <- do.call(rbind, all_summary)
selected <- do.call(rbind, selected_summary)
qc <- do.call(rbind, qc)

# Compare alternative population medians to the proxy population median.
base <- end[end$Trajectory_Mode == "proxy", c("Compound", "Variable", "Median", "N_ID")]
names(base)[3:4] <- c("Proxy_Median", "Proxy_N_ID")
comp <- merge(end, base, by = c("Compound", "Variable"), sort = FALSE)
comp$Delta_vs_Proxy_pct <- ifelse(
  is.finite(comp$Proxy_Median) & comp$Proxy_Median != 0,
  100 * (comp$Median / comp$Proxy_Median - 1), NA_real_
)
comp$Ratio_to_Proxy <- ifelse(
  is.finite(comp$Proxy_Median) & comp$Proxy_Median != 0,
  comp$Median / comp$Proxy_Median, NA_real_
)
comp <- comp[order(match(comp$Trajectory_Mode, MODES),
                   match(comp$Compound, COMPOUNDS),
                   match(comp$Variable, names(measurements))), , drop = FALSE]

# Subject-paired endpoint comparison avoids confusing a ratio of medians
# with the median of within-infant ratios (they are different statistics).
proxy_end <- all_endpoint[[match("proxy", MODES)]]
paired_results <- list()
for (m in setdiff(MODES, "proxy")) {
  alt <- all_endpoint[[match(m, MODES)]]
  a <- merge(proxy_end[, c("ID", "Compound", unname(measurements))],
             alt[, c("ID", "Compound", unname(measurements))],
             by = c("ID", "Compound"), suffixes = c("_proxy", "_alt"))
  if (nrow(a) != length(COMPOUNDS) * EXPECTED_N) {
    stop("Paired ID mismatch for mode ", m, ": matched ", nrow(a), " expected ",
         length(COMPOUNDS) * EXPECTED_N)
  }
  for (cmp in COMPOUNDS) {
    one <- a[a$Compound == cmp, , drop = FALSE]
    for (meas in names(measurements)) {
      fld <- measurements[[meas]]
      x <- one[[paste0(fld,"_proxy")]]
      y <- one[[paste0(fld,"_alt")]]
      valid <- is.finite(x) & is.finite(y) & x > 0
      dif <- if (any(valid)) 100 * (y[valid] / x[valid] - 1) else numeric()
      qs <- summarize_dist(dif)
      paired_results[[length(paired_results) + 1L]] <- data.frame(
        Trajectory_Mode = m, Compound = cmp, Variable = meas,
        N_Paired = sum(valid),
        Median_Within_ID_Change_pct = unname(qs["Median"]),
        P025_Within_ID_Change_pct = unname(qs["P025"]),
        P975_Within_ID_Change_pct = unname(qs["P975"]),
        stringsAsFactors = FALSE
      )
    }
  }
}
paired <- do.call(rbind, paired_results)

# Central scientific assumption: only these three compounds change their
# source milk beta in response to mode, not the PBPK's renal mechanisms.
assumptions <- data.frame(
  Trajectory_Mode = MODES,
  Relative_Beta_Multiplier_For_PFNA_PFDA_PFUdA = c(1, 0, 0.5, 1.5),
  Affected_Compounds = "PFNA; PFDA; PFUdA",
  Fixed_Compounds = "PFOA; PFHxS; PFOS",
  Note = c(
    "Primary PFOA-derived proxy monthly decay",
    "No monthly breast-milk decline (beta=0)",
    "Half of PFOA-derived beta",
    "1.5 times PFOA-derived beta"
  ), stringsAsFactors = FALSE
)
# Detect accidental unintended changes in compounds whose milk beta should
# be identical in all modes, allowing trivial solver noise.
control <- comp[comp$Compound %in% setdiff(COMPOUNDS, PROXY_COMPOUNDS) &
                  comp$Trajectory_Mode != "proxy", , drop = FALSE]
control$Invariant_Within_0_1pct <- is.finite(control$Delta_vs_Proxy_pct) &
  abs(control$Delta_vs_Proxy_pct) <= 0.1

interpretation <- data.frame(
  Topic = c("Purpose", "Reference mode", "Primary scenario",
            "Statistical estimand", "Paired ratio", "Fixed-chemical QC",
            "External validation", "Important limitations"),
  Description = c(
    "Full independent S1 PBPK runs under four beta modes; actual simulated organ concentrations are used.",
    "Proxy uses PFOA-derived decay for PFNA/PFDA/PFUdA; PFOA/PFHxS/PFOS retain their original S1 beta.",
    "This is assumption sensitivity, not four externally calibrated model variants.",
    "Delta_vs_Proxy_pct = 100*(median across 175 infants in alternative / median across 175 infants in proxy - 1).",
    "Within-ID paired change quantifies median and 2.5-97.5% of infant-specific percentage changes.",
    "Unchanged compounds should not materially vary across beta modes; examine Control_Compounds_QC.",
    "This analysis alone does not externally validate longitudinal milk beta parameters.",
    "Predicted 2.5-97.5% are cohort percentiles, not parameter uncertainty or model 95% confidence intervals."
  ), stringsAsFactors = FALSE
)

out_file <- file.path(out_dir, "beta_proxy_sensitivity_summary.xlsx")
sheets <- list(
  Run_Status = status,
  Endpoint_By_Mode = end,
  Relative_To_Proxy = comp,
  Paired_Infant_Change = paired,
  Selected_Days = selected,
  Control_Compounds_QC = control,
  Numerical_QC = qc,
  Beta_Assumptions = assumptions,
  Interpretation = interpretation
)
stopifnot(all(nchar(names(sheets)) <= 31L))
writexl::write_xlsx(sheets, out_file)
utils::write.csv(comp, file.path(out_dir, "beta_proxy_relative_to_proxy.csv"), row.names = FALSE)
utils::write.csv(paired, file.path(out_dir, "beta_proxy_paired_infant_change.csv"), row.names = FALSE)

message("\n===== COMPLETE FOUR-MODE ANALYSIS =====")
message("Workbook: ", out_file)
message("Modes: ", paste(unique(end$Trajectory_Mode), collapse = ", "))
message("Output includes plasma, liver, kidney and body burden.")
message("Run status: ", status_path)
print(comp[comp$Compound %in% PROXY_COMPOUNDS &
             comp$Variable %in% c("Plasma", "Liver", "Kidney"),
           c("Trajectory_Mode", "Compound", "Variable", "Median", "Delta_vs_Proxy_pct")],
      row.names = FALSE)
