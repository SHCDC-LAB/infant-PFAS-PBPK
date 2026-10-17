# FAST3: checkpointed, subject-level, equivalent mathematical model.
# Data/settings checksum guards checkpoint reuse; stop/resume supported.
# NOTE: FAST3 requires independent R execution and numerical equivalence QC.
# Run with PFAS_ATTR_MAX_SUBJECTS=2 for quick end-to-end smoke test.

# PERFORMANCE-OPTIMIZED VERSION (2026-10-17)
# Exact same PBPK equations, LSODA tolerances, daily output times, sensitivity
# grids, Shapley source-attribution definitions and output file structure.
# PRENATAL_ONLY and POSTNATAL_ONLY counterfactual ODE results are reused when
# their actual chemical inputs are identical within one subject/compound job.
# TOTAL runs are not cached. No sensitivities are dropped or approximated.
# A byte-for-byte comparison to prior frozen output has NOT been executed here.
# ===========================================================================

# ============================================================
# PFAS infant PBPK model — FINAL FREEZE v1.1
# Frozen after adult renal qualification, infant numerical QC,
# study-specific Varsi evaluation, prenatal/postnatal counterfactual,
# exposure×Tm sensitivity, and targeted PFNA/PFHxS R_BM audit.
#
# KEY STRUCTURAL REVISION
# -----------------------
# The renal module now follows the PFAS PBPK-template architecture:
#
#   Plasma --Qfil--> renal filtrate --kust--> urinary storage --ku--> urine
#                  |                      
#                  +-- saturable Tm/Kt reabsorption --> Plasma
#
# Thus the chemical-specific urinary elimination parameter (KurineC)
# is NOT used directly on Afiltrate.  A separate filtrate-to-urinary-
# storage transfer rate (kust) is used.  KurineC controls elimination
# from the urinary-storage compartment.
#
# Main features
#   • S0: constant breast-milk concentration
#   • S1: PFAS-specific longitudinal milk trajectory
#   • PFOA, PFNA, PFDA, PFUdA(PFUnDA), PFHxS, PFOS
#   • One common renal Tm/Kt framework for all six PFAS
#   • Age-dependent infant physiology and GFR
#   • BW^0.75 scaling for Tm
#   • BW^-0.25 scaling for Ku and Kust coefficients
#   • Separate urinary-storage compartment
#   • Strict mass balance including urinary storage
#   • Adult-derived/uncertain parameters remain explicitly labeled
#   • PFDA/PFUdA do NOT use GFR/Vfil as chemical-specific Ku
#
# IMPORTANT SCIENTIFIC NOTE
# --------------------------
# The Sweeney/Bernstein PFAS PBPK template describes filtrate -> urinary
# storage using a first-order kust process and urinary-storage -> urine
# using ku. The template also notes that its kidney representation is a
# useful common mathematical framework rather than a fully physiological
# representation of renal microanatomy. We therefore retain the unified
# template architecture while avoiding a compound-specific OAT/PTC
# compartment in the base model.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ============================================================
# v3 CHANGELOG — READ THIS FIRST
# ============================================================
# 1) GFR corrected to the published pediatric BSA equation:
#      GFR = -17.74 + 99.054*BSA - 6.1604*BSA^2   [mL/min]
#    The extra ×0.06 factor in v2 was removed.
# 2) Placental transfer R_CM is explicitly sourced to Liu et al. 2024
#    J Hazard Mater 466:133501, Table 3 median paired whole-blood ratios.
# 3) FINAL FREEZE R_BM: PFNA and PFHxS updated to independent pooled-literature
#    milk:serum ratios (0.0300 and 0.0133); PFOA/PFDA/PFUdA/PFOS retain the
#    prior v3 chemical-specific R_BM values. No infant serum outcome was used
#    to fit these final transfer ratios.
# 4) PFDA/PFUdA effective Tmc values are updated from adult benchmark v5:
#      PFDA  = 14.700113 ug/h/kg^0.75
#      PFUdA = 260.217197 ug/h/kg^0.75
# 5) PFDA/PFUdA KurineC remains NA because adult v5 showed it is not
#    identifiable from plasma kinetics in the current urinary-storage structure.
# 6) Mass-balance bookkeeping is corrected: cumulative breast-milk input is now
#    an auxiliary ODE state (Ainput), integrated by LSODA. Because cumulative
#    urine/feces are already included in A_Total_Model, they are not subtracted
#    a second time from Mass_Balance_Expected.
# 7) Safe summary helpers return NA (not -Inf) when urinary-clearance diagnostics
#    are entirely missing, as expected for PFDA/PFUdA while KurineC is unresolved.
# 8) All other conflicts with the old manuscript draft are resolved in favor
#    of the current code structure/parameters. No OAT ontogeny or old Vmax_app
#    infant calibration has been reintroduced.
# ============================================================

# 7) Infant body-weight individualization revised:
#      • Measured birth weight is matched exactly at age 0.
#      • The proportional birth-weight deviation from the sex-specific reference
#        growth curve decays exponentially during infancy rather than remaining
#        fixed for the whole first year.
#      • Primary assumption: 10% of the initial proportional deviation remains
#        at 12 months (explicit model assumption; sensitivity values exported).
#      • This choice is motivated by evidence that LGA infants commonly show
#        catch-down/regression toward the mean, predominantly during the first
#        6 months; the exact 10% value is not claimed as a directly measured
#        biological constant.
# 8) Birth-weight flags added for sensitivity analysis:
#      • >=4.0 kg: macrosomia flag (commonly used absolute threshold).
#      • >=4.5 kg: extreme/high macrosomia sensitivity threshold. ACOG PB 216
#        stratifies macrosomia as 4,000–4,499 g, 4,500–4,999 g and >=5,000 g;
#        4.5 kg is therefore used here as the prespecified extreme-BW threshold.
#      Extreme-BW subjects are retained in the primary analysis and excluded
#      only in a secondary influence/sensitivity summary.
# ============================================================

# ============================================================
# 0. PACKAGES
# ============================================================
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "deSolve",
  "writexl", "future", "future.apply"
)

new_pkgs <- required_pkgs[
  !(required_pkgs %in% installed.packages()[, "Package"])
]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(deSolve)
  library(writexl)
  library(future)
  library(future.apply)
})

# ============================================================
# 1. USER SETTINGS
# ============================================================
PROJECT_ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", unset = getwd()),
                              winslash = "/", mustWork = FALSE)
file_input <- Sys.getenv("PFAS_PBPK_INPUT", unset = file.path(PROJECT_ROOT,
                                                        "data", "input-PBPK.xlsx"))
file_abraham_reference <- "E:/桌面2/PFAS-母乳PBPK/文献/1-s2.0-S0160412024006330-mmc2.xlsx"

BASE_RESULTS_DIR <- file.path(PROJECT_ROOT, "results", "source_attribution_setup")
MODEL_FREEZE_VERSION <- "PFAS_RCM_RBM_source_attribution_sensitivity_v1_setup"
MODEL_FREEZE_DATE <- "2026-09-09"
OUTPUT_DIR <- file.path(BASE_RESULTS_DIR, MODEL_FREEZE_VERSION)

INPUT_MILK_UNIT <- "ng/L"
INPUT_BW_UNIT <- "g"

SIM_YEARS <- 1
TIME_STEP_HOURS <- 168       # 7-day reporting grid; LSODA remains adaptive

# ---- fast engine ----
USE_PARALLEL <- TRUE
DETECT_PHYSICAL_CORES <- TRUE
PARALLEL_WORKERS_MANUAL <- 4L

PARALLEL_WORKERS <- if (DETECT_PHYSICAL_CORES) {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
} else {
  PARALLEL_WORKERS_MANUAL
}
if (!is.finite(PARALLEL_WORKERS) || PARALLEL_WORKERS < 1L) {
  PARALLEL_WORKERS <- 1L
}

# Daily age grid for cached physiology and renal parameters.
RENAL_CACHE_STEP_DAYS <- 1

# Only endpoint is currently needed for contribution analysis.
CONTRIBUTION_MONTHS <- 12

# ODE tolerances.
ODE_RTOL <- 1e-6
ODE_ATOL <- 1e-8

# ---- renal model ----
TM_SCALING_EXPONENT <- 0.75
KU_SCALING_EXPONENT <- -0.25
KUST_SCALING_EXPONENT <- -0.25

# Effective filtrate volume used to compute the transit rate into urinary storage.
# IMPORTANT: this is NOT the chemical-specific urinary elimination constant.
VFIL_L_PER_KG_BW <- 4.0e-4

# Qfil mode. The current infant physiology provides developmental GFR.
RENAL_FILTRATION_FLOW_MODE <- "GFR"

# kust describes filtrate -> urinary storage. A physiological residence-rate
# approximation is used because a chemical-specific kust is not available here:
#     kust = Qfil / Vfil
# This is deliberately kept separate from chemical-specific KurineC.
KUST_MODE <- "QFIL_OVER_VFIL"

# Literature/derived chemical-specific urinary-storage elimination rate
# KurineC is used only when available. PFDA/PFUdA remain NA here and are
# therefore flagged as uncertain/calibration candidates rather than receiving
# an artificial GFR/Vfil urinary elimination rate.
ALLOW_MISSING_KU <- TRUE

# ---- QC thresholds ----
BW_BIRTH_MEDIAN_MAX_KG <- 10
MIN_GFR_ML_MIN <- 1e-8
MIN_VOLUME_L <- 1e-10
MASS_BALANCE_REL_TOL <- 1e-5
REABSORPTION_RATIO_WARNING <- 1.01
TIME_GRID_REL_TOL <- 1e-3

# ---- infant growth individualization ----
# Primary model: measured birth weight is matched exactly at age 0, then the
# proportional deviation from the sex-specific reference curve decays toward 1.
# Evidence base for direction/timing: LGA infants frequently show catch-down
# growth and regression toward the mean, with most change occurring before
# 6 months (e.g., Zhou & Guo 2024; Dunn et al. 2021).
# The exact remaining fraction is a transparent model assumption, not a directly
# measured PFAS/PBPK parameter. It is therefore accompanied by sensitivity output.
BW_INDIVIDUALIZATION_MODE <- "DECAY_TO_REFERENCE"
BW_BIRTH_EFFECT_REMAINING_AT_12MO <- 0.10
BW_DECAY_K_PER_YEAR <- -log(BW_BIRTH_EFFECT_REMAINING_AT_12MO)
BW_REMAINING_FRACTION_SENSITIVITY <- c(0.05, 0.10, 0.25, 1.00)

# Birth-weight classification / influence thresholds.
# 4.0 kg is a common absolute macrosomia threshold. ACOG Practice Bulletin 216
# presents risk categories 4,000–4,499 g, 4,500–4,999 g, and >=5,000 g.
# We therefore prespecify >=4.5 kg as the 'extreme high birth weight' threshold
# for sensitivity/influence analysis. Subjects are NOT removed from the primary run.
MACROSOMIA_THRESHOLD_KG <- 4.0
EXTREME_HIGH_BW_THRESHOLD_KG <- 4.5
VERY_EXTREME_BW_THRESHOLD_KG <- 5.0

# ---- plotting ----
PLOT_DPI <- 600
PLOT_WIDTH <- 9
PLOT_HEIGHT <- 6.5

SCENARIOS <- c("S0", "S1")
pfas_cols <- c("PFOA", "PFNA", "PFDA", "PFUdA", "PFHxS", "PFOS")

# ---- output directories ----
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "milk_trajectories"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "qc"), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1B. SAFE SUMMARY UTILITIES
# ============================================================
# These helpers return NA rather than +/-Inf when a whole group contains
# only NA/non-finite values. This is expected for PFDA/PFUdA urinary
# clearance diagnostics while chemical-specific KurineC remains unresolved.
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  
  as.numeric(
    stats::quantile(
      x,
      probs = p,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )
}

safe_write_xlsx <- function(x, path) {
  
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  excel_limit <- 1048576L
  
  nr <- vapply(
    x,
    function(z) {
      if (is.data.frame(z)) nrow(z) else 0L
    },
    integer(1)
  )
  
  oversized <- names(nr)[nr + 1L > excel_limit]
  
  if (length(oversized)) {
    stop(
      "Excel row limit exceeded in sheet(s): ",
      paste(oversized, collapse = ", "),
      ". Split the data before export."
    )
  }
  
  ok <- tryCatch(
    {
      writexl::write_xlsx(
        x = x,
        path = path
      )
      TRUE
    },
    error = function(e) {
      message(
        "write_xlsx failed for: ",
        path,
        "\nReason: ",
        conditionMessage(e)
      )
      FALSE
    }
  )
  
  if (ok) return(path)
  
  stamp <- format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )
  
  alt_path <- file.path(
    dirname(path),
    paste0(
      tools::file_path_sans_ext(
        basename(path)
      ),
      "_",
      stamp,
      ".xlsx"
    )
  )
  
  writexl::write_xlsx(
    x = x,
    path = alt_path
  )
  
  warning(
    "Could not write requested workbook. Saved instead to:\n",
    alt_path
  )
  
  alt_path
}

# ============================================================
# 2. LITERATURE-INFORMED MILK TRAJECTORIES
# ============================================================
trajectory_params <- tibble::tribble(
  ~Compound, ~supported_by_blomberg, ~early_change_mean, ~early_ci_low, ~early_ci_high,
  ~monthly_change_mean, ~monthly_ci_low, ~monthly_ci_high, ~source_note,
  "PFOA",  TRUE,  -0.167, -0.280, -0.035, -0.120, -0.220, -0.015,
  "Blomberg 2023 Tables 3–4",
  "PFNA",  FALSE, 0.000, NA, NA, 0.000, NA, NA,
  "No supported longitudinal trajectory estimate in Blomberg 2023; flat scenario fallback",
  "PFDA",  FALSE, 0.000, NA, NA, 0.000, NA, NA,
  "No supported longitudinal trajectory estimate in Blomberg 2023; flat scenario fallback",
  "PFUdA", FALSE, 0.000, NA, NA, 0.000, NA, NA,
  "No supported longitudinal trajectory estimate in Blomberg 2023; flat scenario fallback",
  "PFHxS", TRUE, -0.116, -0.244, 0.033, -0.025, -0.097, 0.054,
  "Blomberg 2023 Tables 3–4",
  "PFOS",  TRUE,  0.211,  0.089, 0.347, -0.032, -0.084, 0.023,
  "Blomberg 2023 Tables 3–4"
)

DAYS_PER_MONTH <- 30.4375
monthly_to_beta <- function(monthly_change) {
  ifelse(
    is.finite(monthly_change) & monthly_change > -1,
    log1p(monthly_change) / DAYS_PER_MONTH,
    NA_real_
  )
}

trajectory_params <- trajectory_params %>%
  mutate(
    beta_mean_day = monthly_to_beta(monthly_change_mean),
    beta_low_day = monthly_to_beta(monthly_ci_low),
    beta_high_day = monthly_to_beta(monthly_ci_high)
  )

write_xlsx(
  list(Milk_Trajectory_Parameters = trajectory_params),
  path = file.path(OUTPUT_DIR, "tables", "milk_trajectory_parameters.xlsx")
)

# ============================================================
# 3. UNIFIED RENAL Tm/Kt PARAMETER SET
# ============================================================
# Tm(BW) = Tmc * BW^0.75
#
# PFDA/PFUdA starting Tm values were Fàbrega-derived, but v3 carries forward
# the effective Tmc values minimally calibrated in the adult v5 benchmark
# against Abraham 2024 terminal plasma kinetics (primary terminal window >30 d):
#   PFDA  Tmc = 14.700113 ug/h/kg^0.75
#   PFUdA Tmc = 260.217197 ug/h/kg^0.75
# Kt and Free remain fixed at the Fàbrega-derived values.
#
# IMPORTANT: PFDA/PFUdA KurineC is NOT identifiable from plasma kinetics in
# the current Afil -> Aust -> urine structure. It therefore remains NA and is
# not interpreted quantitatively for urinary excretion in the base infant model.
# ============================================================
TMKT_PARAMS <- tibble::tribble(
  ~Compound, ~Tmc_ug_h_bw075, ~Tm_basis, ~Kt_ug_L, ~Free_fraction,
  ~KurineC_h_bwminus025, ~Kurine_source, ~Tm_source, ~Parameter_class,
  
  "PFOA", 6000.0000, "literature Tmc", 55.0, 0.020,
  0.0003, "EFSA 2018 / PFAS PBK implementation", "EFSA 2018", "human-calibrated",
  
  "PFOS", 3500.0000, "literature Tmc", 23.0, 0.025,
  0.0010, "EFSA 2018 / PFAS PBK implementation", "EFSA 2018", "human-calibrated",
  
  "PFNA", 7900.0000, "literature/calibrated Tmc", 55.0, 0.020,
  0.0003, "Fragki 2023 / Iulini 2024; based on PFOA", "Fragki 2023 / Iulini 2024", "human half-life calibrated",
  
  "PFHxS", 7000.0000, "literature/calibrated Tmc", 23.0, 0.025,
  0.0010, "Fragki 2023 / Iulini 2024; based on PFOS", "Fragki 2023 / Iulini 2024", "human half-life calibrated",
  
  "PFDA", 14.700113, "Fàbrega-informed; Abraham 2024 adult-v5 calibrated effective Tmc", 0.6, 0.010,
  NA_real_, "Not identifiable from plasma; urine/fecal data required", "Fàbrega 2015 + Abraham 2024 adult v5", "current-study calibrated Tmc",
  
  "PFUdA", 260.217197, "Fàbrega-informed; Abraham 2024 adult-v5 calibrated effective Tmc", 9.0, 0.030,
  NA_real_, "Not identifiable from plasma; urine/fecal data required", "Fàbrega 2015 + Abraham 2024 adult v5", "current-study calibrated Tmc"
)

write_xlsx(
  list(Unified_TmKt_Parameters = TMKT_PARAMS),
  path = file.path(OUTPUT_DIR, "tables", "unified_TmKt_parameters.xlsx")
)

# ============================================================
# 4. INFANT PHYSIOLOGICAL PARAMETERS
# ============================================================
# Reference growth curves are the same sex-specific equations used previously.
# Individualization no longer multiplies the entire first-year trajectory by a
# constant birth-weight ratio. Instead:
#   birth_ratio = measured_BW_birth / reference_BW_at_birth
#   ratio(age)   = 1 + (birth_ratio - 1) * exp(-k * age_years)
#   BW(age)      = reference_BW(age) * ratio(age)
# This guarantees BW(0) = measured birth weight and allows catch-down/catch-up
# toward the reference trajectory during infancy.

get_reference_height_weight <- function(age_yr, gender_code) {
  gender <- ifelse(as.numeric(gender_code) == 1, "Male", "Female")
  
  if (gender == "Male") {
    HT <- -0.000347983 * age_yr^6 +
      0.0168412 * age_yr^5 -
      0.315285 * age_yr^4 +
      2.8874 * age_yr^3 -
      13.5636 * age_yr^2 +
      37.5714 * age_yr + 50.4
    
    WT_base <- 8.0384 * (1.0 - exp(age_yr * -0.7657)) +
      exp((HT * 0.02184) + (0.0135 * age_yr))
  } else {
    HT <- -0.000436143 * age_yr^6 +
      0.0194832 * age_yr^5 -
      0.340837 * age_yr^4 +
      2.95769 * age_yr^3 -
      13.3793 * age_yr^2 +
      36.6627 * age_yr + 49.7
    
    WT_base <- 4.1974 * (1.0 - exp(age_yr * -2.7662)) +
      exp((HT * 0.02318) + (0.0093 * age_yr))
  }
  
  list(gender = gender, HT = HT, WT_base = WT_base)
}

get_growth_adjustment <- function(age_yr, gender_code, BW_birth_kg,
                                  remaining_at_12mo = BW_BIRTH_EFFECT_REMAINING_AT_12MO,
                                  mode = BW_INDIVIDUALIZATION_MODE) {
  ref_now <- get_reference_height_weight(age_yr, gender_code)
  ref_birth <- get_reference_height_weight(0, gender_code)$WT_base
  
  if (!is.finite(BW_birth_kg) || BW_birth_kg <= 0) {
    stop("Invalid measured birth weight supplied to physiology function.")
  }
  
  birth_ratio <- BW_birth_kg / ref_birth
  
  if (mode == "DECAY_TO_REFERENCE") {
    if (!is.finite(remaining_at_12mo) || remaining_at_12mo <= 0 || remaining_at_12mo > 1) {
      stop("remaining_at_12mo must be in (0, 1].")
    }
    k <- -log(remaining_at_12mo)
    ratio_age <- 1 + (birth_ratio - 1) * exp(-k * age_yr)
  } else if (mode == "CONSTANT_RATIO") {
    ratio_age <- birth_ratio
  } else if (mode == "REFERENCE_ONLY") {
    ratio_age <- 1
  } else {
    stop("Unknown BW_INDIVIDUALIZATION_MODE: ", mode)
  }
  
  list(
    reference_BW_kg = ref_now$WT_base,
    reference_birth_BW_kg = ref_birth,
    birth_ratio_to_reference = birth_ratio,
    adjustment_ratio = ratio_age,
    individualized_BW_kg = ref_now$WT_base * ratio_age
  )
}

get_infant_phys_meta <- function(age_yr, gender_code, BW_birth_kg,
                                 remaining_at_12mo = BW_BIRTH_EFFECT_REMAINING_AT_12MO,
                                 mode = BW_INDIVIDUALIZATION_MODE) {
  PNA_days <- max(1, age_yr * 365)
  ref <- get_reference_height_weight(age_yr, gender_code)
  gender <- ref$gender
  HT <- ref$HT
  WT_base <- ref$WT_base
  
  if (gender == "Male") {
    Hct <- 53 - (
      (43.0 * age_yr^1.12 /
         (0.05^1.12 + age_yr^1.12)) *
        (1 + (-0.93 * age_yr^0.25 /
                (0.10^0.25 + age_yr^0.25)))
    )
  } else {
    Hct <- 53 - (
      (37.4 * age_yr^1.12 /
         (0.05^1.12 + age_yr^1.12)) *
        (1 + (-0.80 * age_yr^0.25 /
                (0.10^0.25 + age_yr^0.25)))
    )
  }
  
  growth <- get_growth_adjustment(
    age_yr = age_yr,
    gender_code = gender_code,
    BW_birth_kg = BW_birth_kg,
    remaining_at_12mo = remaining_at_12mo,
    mode = mode
  )
  WT <- max(growth$individualized_BW_kg, 0.5)
  
  milk_ml_day <- 51 - 1.4 * PNA_days + 180 * log(PNA_days)
  KMilk_L_h <- max(0, milk_ml_day) / 1000 / 24
  
  CO_mL_min <- (214.81 * WT + 76.57) * 0.06
  BSA <- 0.09395 * WT^0.7032
  
  V_li <- ifelse(
    gender == "Male",
    0.6974 * BSA^1.0388,
    0.722 * BSA^1.176
  )
  
  Kidney_WT_g <- (4.214 * WT^0.823) /
    (4.456 + WT^0.795) * 1.05
  
  V_kd <- Kidney_WT_g / 1000 / 1.05
  
  BlV <- ifelse(
    gender == "Male",
    (1.362 * BSA)^1.787 + 0.375,
    (1.714 * BSA)^1.313 + 0.155
  )
  
  V_pl <- BlV * (1 - Hct / 100)
  
  V_re <- max(
    (WT / 1.05) - (V_li + V_kd + V_pl),
    0.10 * WT
  )
  
  if (gender == "Male") {
    RBF_perc <- 7.94 +
      (11 * age_yr^3.34) /
      (0.97^3.34 + age_yr^3.34)
  } else {
    RBF_perc <- 6.35 +
      (10.89 * age_yr^1.26) /
      (1.51^1.26 + age_yr^1.26)
  }
  
  Q_li_mL_min <- 0.26 * CO_mL_min
  Q_kd_mL_min <- (RBF_perc / 100) * CO_mL_min
  Q_re_mL_min <- max(
    CO_mL_min - Q_li_mL_min - Q_kd_mL_min,
    0
  )
  
  # Literature pediatric GFR equation (mL/min):
  #   GFR = -17.74 + 99.054*BSA - 6.1604*BSA^2
  GFR_mL_min <- max(
    -17.74 + 99.054 * BSA - 6.1604 * BSA^2,
    0
  )
  
  data.frame(
    BW = WT,
    BW_reference = WT_base,
    BW_adjustment_ratio = growth$adjustment_ratio,
    BW_birth_ratio_to_reference = growth$birth_ratio_to_reference,
    Hct = Hct,
    BSA = BSA,
    V_li = V_li,
    V_kd = V_kd,
    V_pl = V_pl,
    V_re = V_re,
    CO_mL_min = CO_mL_min,
    Q_li_mL_min = Q_li_mL_min,
    Q_kd_mL_min = Q_kd_mL_min,
    Q_re_mL_min = Q_re_mL_min,
    GFR_mL_min = GFR_mL_min,
    KMilk_L_h = KMilk_L_h
  )
}

# ============================================================
# 5. PLACENTAL TRANSFER / INITIAL BURDEN
# ============================================================
# Transfer-efficiency definitions used for birth-burden reconstruction:
#   R_CM = cord blood / maternal blood (placental transfer efficiency)
#   R_BM = breast milk / maternal serum (breastfeeding transfer efficiency)
#
# R_CM (v3): Liu et al. 2024, Journal of Hazardous Materials 466:133501,
# Table 3 median paired cord-to-maternal whole-blood ratios.
#   PFOA 0.290; PFNA 0.230; PFDA 0.170; PFUnDA 0.230; PFHxS 0.220; PFOS 0.090.
#
# R_BM FINAL-FREEZE DECISION
# ------------------------------------------------------------
# PFOA, PFDA, PFUdA and PFOS retain the prior chemical-specific values used in
# infant v3 because study-specific Varsi evaluation did not identify a strong
# transfer-related failure for those compounds.
#
# PFNA and PFHxS are updated to independent literature-pooled milk:serum ratios:
#   PFNA  0.0300
#   PFHxS 0.0133
#
# These values were selected BEFORE the final freeze because the targeted
# PFNA/PFHxS transfer audit showed:
#   - current v3 R_BM values (PFNA 0.103; PFHxS 0.047) produced large/high
#     Varsi-specific postnatal exposure,
#   - independent paired/pooled literature supported substantially lower ratios,
#   - lowering R_BM improved PFNA/PFHxS predictions without altering renal
#     kinetics, consistent with the earlier exposure×Tm sensitivity analysis.
#
# IMPORTANT:
# These are literature-driven final defaults, NOT values fitted to infant serum.
# The final model does not back-calculate "optimal" R_BM from Varsi outcomes.
#
# In the primary cohort, measured milk PFAS is used both as postnatal exposure
# and, through R_BM, to infer maternal plasma for birth-burden initialization:
#
#   C_mat_plasma = C_milk_measured / R_BM
#   C_fetal      = C_mat_plasma × R_CM
#
# Therefore updating PFNA/PFHxS R_BM also changes inferred prenatal burden in
# the primary cohort. This is intentional and must be reported transparently.
# The Varsi study-specific validation instead initializes birth burden directly
# from maternal serum × R_CM, thereby separating prenatal and postnatal transfer.
#
transfer_data <- data.frame(
  Compound = c("PFOA", "PFNA", "PFDA", "PFUdA", "PFHxS", "PFOS"),
  
  R_CM = c(
    0.290,
    0.230,
    0.170,
    0.230,
    0.220,
    0.090
  ),
  
  R_BM = c(
    0.053,
    0.0300,
    0.023,
    0.022,
    0.0133,
    0.042
  ),
  
  R_CM_Source = rep(
    "Liu et al. 2024 J Hazard Mater 466:133501, Table 3 median RCM",
    6
  ),
  
  R_BM_Source = c(
    "Prior v3 chemical-specific value retained; Zheng et al. 2022-attributed matched maternal serum/breast milk parameterization",
    "FINAL FREEZE: literature-pooled breast milk/serum ratio; pooled central tendency across international studies",
    "Prior v3 chemical-specific value retained; Zheng et al. 2022-attributed matched maternal serum/breast milk parameterization",
    "Prior v3 chemical-specific value retained; Zheng et al. 2022-attributed matched maternal serum/breast milk parameterization",
    "FINAL FREEZE: literature-pooled breast milk/serum ratio; pooled central tendency across international studies",
    "Prior v3 chemical-specific value retained; Zheng et al. 2022-attributed matched maternal serum/breast milk parameterization"
  ),
  
  R_BM_Evidence_Class = c(
    "retained prior parameter",
    "independent pooled literature",
    "retained prior parameter",
    "retained prior parameter",
    "independent pooled literature",
    "retained prior parameter"
  ),
  
  R_CM_Matrix = rep(
    "umbilical cord whole blood / maternal whole blood",
    6
  ),
  
  R_BM_Matrix = rep(
    "breast milk / maternal serum or plasma",
    6
  ),
  
  Matrix_Compatibility_Note = c(
    "Birth-burden reconstruction treats R_CM and R_BM as effective transfer factors from different studies/matrices.",
    "PFNA final R_BM is pooled literature-based; primary-cohort birth burden is inferred from measured milk and therefore changes when R_BM changes.",
    "Birth-burden reconstruction treats R_CM and R_BM as effective transfer factors from different studies/matrices.",
    "Birth-burden reconstruction treats R_CM and R_BM as effective transfer factors from different studies/matrices.",
    "PFHxS final R_BM is pooled literature-based; primary-cohort birth burden is inferred from measured milk and therefore changes when R_BM changes.",
    "Birth-burden reconstruction treats R_CM and R_BM as effective transfer factors from different studies/matrices."
  ),
  
  Freeze_Status = c(
    "FROZEN",
    "FROZEN_UPDATED",
    "FROZEN",
    "FROZEN",
    "FROZEN_UPDATED",
    "FROZEN"
  ),
  
  stringsAsFactors = FALSE
)

write_xlsx(
  list(Transfer_Parameters = transfer_data),
  path = file.path(OUTPUT_DIR, "tables", "maternal_cord_milk_transfer_parameters_v3.xlsx")
)


# ============================================================
# FINAL MODEL FREEZE MANIFEST
# ============================================================

model_freeze_manifest <- tibble::tribble(
  ~Domain, ~Item, ~Final_Value_or_Status, ~Freeze_Rationale,
  
  "Version",
  "Model version",
  MODEL_FREEZE_VERSION,
  "Final primary infant model after staged qualification and targeted transfer audit.",
  
  "Simulation",
  "Primary horizon",
  "0-1 year",
  "Final infant analysis horizon.",
  
  "Simulation",
  "Reporting interval",
  paste0(TIME_STEP_HOURS / 24, " days"),
  "Weekly reporting; LSODA remains adaptive internally.",
  
  "Renal",
  "Renal structure",
  "Filtration -> saturable Tm/Kt reabsorption -> urinary storage -> urine",
  "Unified renal architecture retained after adult benchmark and infant QC.",
  
  "Renal",
  "PFOA/PFOS/PFNA/PFHxS Tm/Kt",
  "Frozen from final infant v3",
  "Not re-fitted to infant validation outcomes.",
  
  "Renal",
  "PFDA effective Tmc",
  "14.700113 ug/h/kg^0.75",
  "Adult plasma-calibrated effective Tmc; window sensitivity previously checked.",
  
  "Renal",
  "PFUdA effective Tmc",
  "260.217197 ug/h/kg^0.75",
  "Adult plasma-calibrated effective Tmc; window sensitivity previously checked.",
  
  "Renal",
  "PFDA/PFUdA Ku",
  "Unresolved / NA",
  "Plasma kinetics remain usable; urinary excretion is not quantitatively interpreted.",
  
  "Physiology",
  "GFR",
  "Published pediatric BSA equation",
  "Corrected equation retained; no extra x0.06 factor.",
  
  "Physiology",
  "Birth-weight individualization",
  "Exponential convergence to sex-reference growth curve",
  "Primary assumption: 10% of proportional birth-weight deviation remains at 12 months.",
  
  "Transfer",
  "R_CM",
  "Liu et al. 2024 values",
  "Frozen direct placental transfer ratios.",
  
  "Transfer",
  "PFNA R_BM",
  "0.0300",
  "Independent pooled-literature milk:serum ratio; targeted audit supported transfer uncertainty as main PFNA mismatch driver.",
  
  "Transfer",
  "PFHxS R_BM",
  "0.0133",
  "Independent pooled-literature milk:serum ratio; targeted audit supported transfer uncertainty as important PFHxS mismatch driver.",
  
  "Transfer",
  "PFOA/PFDA/PFUdA/PFOS R_BM",
  "0.053 / 0.023 / 0.022 / 0.042",
  "Prior v3 values retained because no targeted update was justified.",
  
  "Exposure",
  "S0",
  "Constant measured milk concentration",
  "Frozen primary scenario.",
  
  "Exposure",
  "S1",
  "Literature-informed longitudinal milk trajectory",
  "PFOA/PFHxS/PFOS use literature trend; PFNA/PFDA/PFUdA remain flat where compound-specific longitudinal evidence is unavailable.",
  
  "Numerics",
  "Mass balance",
  "Ainput integrated as ODE state",
  "Exact cumulative input bookkeeping retained.",
  
  "Qualification",
  "Adult renal benchmark",
  "Completed",
  "Adult terminal kinetics used to qualify/finalize renal parameterization.",
  
  "Qualification",
  "Infant numerical QC",
  "Completed",
  "Mass balance, growth, negative-state and transient reabsorption diagnostics completed.",
  
  "Qualification",
  "Prenatal/postnatal counterfactual",
  "Completed",
  "Postnatal exposure becomes dominant during early infancy for most PFAS.",
  
  "Qualification",
  "Exposure x Tm sensitivity",
  "Completed",
  "Plasma predictions are substantially more sensitive to postnatal exposure than to Tm scaling over the explored range.",
  
  "Qualification",
  "Varsi study-specific evaluation",
  "Completed before final transfer update",
  "Four PFAS reproduced reasonably; PFNA/PFHxS mismatches triggered targeted transfer audit.",
  
  "Qualification",
  "PFNA/PFHxS transfer audit",
  "Completed",
  "Independent literature-based lower R_BM values substantially improved PFNA/PFHxS exposure consistency without changing renal kinetics.",
  
  "Freeze rule",
  "Further parameter tuning",
  "PROHIBITED unless new independent evidence justifies revision",
  "Prevents validation-driven overfitting after final freeze."
)

write_xlsx(
  list(
    Model_Freeze_Manifest = model_freeze_manifest,
    Final_Transfer_Parameters = transfer_data
  ),
  path = file.path(
    OUTPUT_DIR,
    "tables",
    "model_freeze_manifest_v1.xlsx"
  )
)

pfas_params <- list(
  PFOA  = c(P_li = 2.20, P_kd = 1.05, P_re = 0.12),
  PFNA  = c(P_li = 1.46, P_kd = 0.60, P_re = 0.12),
  PFDA  = c(P_li = 0.412, P_kd = 0.550, P_re = 0.334),
  PFUdA = c(P_li = 0.512, P_kd = 0.600, P_re = 0.442),
  PFHxS = c(P_li = 0.85, P_kd = 0.30, P_re = 0.20),
  PFOS  = c(P_li = 3.72, P_kd = 0.80, P_re = 0.20)
)

get_infant_initial_burden <- function(
    compound_name,
    C_milk_measured,
    V_pl_0, V_li_0, V_kd_0, V_re_0,
    P_li, P_kd, P_re
) {
  params <- transfer_data[
    transfer_data$Compound == compound_name,
    , drop = FALSE
  ]
  
  if (nrow(params) != 1) {
    stop("Compound not found in transfer_data: ", compound_name)
  }
  
  C_mat_plasma <- C_milk_measured / params$R_BM
  C_fetal_plasma <- C_mat_plasma * params$R_CM
  
  c(
    Al = C_fetal_plasma * P_li * V_li_0,
    Ak = C_fetal_plasma * P_kd * V_kd_0,
    Ar = C_fetal_plasma * P_re * V_re_0,
    Ap = C_fetal_plasma * V_pl_0
  )
}

# ============================================================
# 6. UNIT CONVERSION
# ============================================================
convert_milk_to_ng_L <- function(x, unit) {
  unit <- tolower(trimws(unit))
  if (unit == "ng/l") return(as.numeric(x))
  if (unit %in% c("ug/l", "µg/l", "μg/l")) return(as.numeric(x) * 1000)
  if (unit == "pg/ml") return(as.numeric(x))
  if (unit == "pg/g") return(as.numeric(x))
  stop("Unsupported INPUT_MILK_UNIT = ", unit)
}

convert_bw_to_kg <- function(x, unit) {
  unit <- tolower(trimws(unit))
  if (unit %in% c("kg", "kgs", "kilogram", "kilograms")) return(as.numeric(x))
  if (unit %in% c("g", "gram", "grams")) return(as.numeric(x) / 1000)
  stop("Unsupported INPUT_BW_UNIT = ", unit, ". Use 'g' or 'kg'.")
}

parse_numeric_censored <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x_chr <- trimws(as.character(x))
  out <- rep(NA_real_, length(x_chr))
  
  for (j in seq_along(x_chr)) {
    val <- x_chr[j]
    if (is.na(val) || val == "") next
    
    if (grepl("^<", val)) {
      num <- suppressWarnings(as.numeric(sub("^<\\s*", "", val)))
      out[j] <- if (is.finite(num)) num / sqrt(2) else NA_real_
    } else {
      val2 <- gsub("^>\\s*", "", val)
      out[j] <- suppressWarnings(as.numeric(val2))
    }
  }
  
  out
}

# ============================================================
# 7. READ / CLEAN INPUT
# ============================================================
if (!file.exists(file_input)) {
  stop("PBPK input file not found: ", file_input)
}

raw_data <- read_excel(file_input, sheet = "Sheet1")
required_input_cols <- c("BW", "infsex", pfas_cols)
missing_input_cols <- setdiff(required_input_cols, names(raw_data))

if (length(missing_input_cols) > 0) {
  stop(
    "Missing required input columns: ",
    paste(missing_input_cols, collapse = ", ")
  )
}

for (col in pfas_cols) {
  raw_data[[col]] <- parse_numeric_censored(raw_data[[col]])
}

raw_data$BW <- suppressWarnings(as.numeric(raw_data$BW))
raw_data$infsex <- suppressWarnings(as.numeric(raw_data$infsex))
raw_data$BW_kg <- convert_bw_to_kg(raw_data$BW, INPUT_BW_UNIT)

mean_BW_by_sex <- raw_data %>%
  group_by(infsex) %>%
  summarise(mean_BW_kg = mean(BW_kg, na.rm = TRUE), .groups = "drop")

raw_data <- raw_data %>%
  left_join(mean_BW_by_sex, by = "infsex") %>%
  mutate(BW_kg = ifelse(is.na(BW_kg), mean_BW_kg, BW_kg)) %>%
  select(-mean_BW_kg)

cleaned_data <- raw_data %>%
  mutate(
    BW_birth_kg = BW_kg,
    # Legacy ratio retained for audit only; it is NOT used for longitudinal growth.
    wt_ratio_legacy = BW_birth_kg / ifelse(infsex == 1, 3.35, 3.22),
    Macrosomia_4kg_Flag = BW_birth_kg >= MACROSOMIA_THRESHOLD_KG,
    Extreme_High_BW_4p5kg_Flag = BW_birth_kg >= EXTREME_HIGH_BW_THRESHOLD_KG,
    Very_Extreme_BW_5kg_Flag = BW_birth_kg >= VERY_EXTREME_BW_THRESHOLD_KG,
    ID = row_number(),
    mix_ratio = 1
  ) %>%
  mutate(
    across(all_of(pfas_cols), ~ convert_milk_to_ng_L(.x, INPUT_MILK_UNIT))
  )

if ("BFpc" %in% names(cleaned_data)) {
  cleaned_data$mix_ratio <- suppressWarnings(as.numeric(cleaned_data$BFpc))
  cleaned_data$mix_ratio[
    is.na(cleaned_data$mix_ratio) |
      cleaned_data$mix_ratio < 0 |
      cleaned_data$mix_ratio > 1
  ] <- 1
}

cleaned_data <- cleaned_data %>%
  filter(
    is.finite(BW_birth_kg),
    BW_birth_kg > 0,
    is.finite(infsex),
    is.finite(BW_birth_kg)
  ) %>%
  filter(if_any(all_of(pfas_cols), ~ !is.na(.x)))

if (nrow(cleaned_data) == 0) {
  stop("No usable rows remain after input cleaning.")
}

# ============================================================
# 8. PRE-RUN QC
# ============================================================
cat("\n============================================================\n")
cat("PRE-RUN BODY-WEIGHT / PHYSIOLOGY QC\n")
cat("============================================================\n")
cat("INPUT_BW_UNIT = ", INPUT_BW_UNIT, "\n", sep = "")
cat(
  "Birth BW: median = ",
  round(median(cleaned_data$BW_birth_kg, na.rm = TRUE), 3),
  " kg; range = ",
  round(min(cleaned_data$BW_birth_kg, na.rm = TRUE), 3),
  "–",
  round(max(cleaned_data$BW_birth_kg, na.rm = TRUE), 3),
  " kg\n",
  sep = ""
)

if (median(cleaned_data$BW_birth_kg, na.rm = TRUE) > BW_BIRTH_MEDIAN_MAX_KG) {
  stop(
    "FATAL QC: median birth BW = ",
    median(cleaned_data$BW_birth_kg, na.rm = TRUE),
    " kg. Check INPUT_BW_UNIT."
  )
}

if (any(cleaned_data$BW_birth_kg <= 0 | !is.finite(cleaned_data$BW_birth_kg))) {
  stop("FATAL QC: non-positive/non-finite birth BW detected.")
}

qc_sub <- cleaned_data[1, , drop = FALSE]
qc_phys <- get_infant_phys_meta(0, qc_sub$infsex, qc_sub$BW_birth_kg)

cat("ID 1 / first subject physiology at birth:\n")
cat("  BW  = ", round(qc_phys$BW, 4), " kg\n", sep = "")
cat("  BSA = ", round(qc_phys$BSA, 5), " m2\n", sep = "")
cat("  Vpl = ", round(qc_phys$V_pl, 6), " L\n", sep = "")
cat("  Vkd = ", round(qc_phys$V_kd, 6), " L\n", sep = "")
cat("  GFR = ", round(qc_phys$GFR_mL_min, 6), " mL/min\n", sep = "")

if (
  qc_phys$BW <= 0 ||
  qc_phys$BW > 10 ||
  qc_phys$BSA <= 0 ||
  qc_phys$V_pl <= MIN_VOLUME_L ||
  qc_phys$V_kd <= MIN_VOLUME_L ||
  qc_phys$GFR_mL_min <= MIN_GFR_ML_MIN
) {
  stop("FATAL QC: implausible infant physiology at birth.")
}

cat("Pre-run physiology QC passed.\n")

# Birth-weight classification summary and 12-month growth QC.
bw_class_summary <- tibble::tibble(
  Classification = c(
    "All",
    ">=4.0 kg macrosomia",
    ">=4.5 kg extreme/high macrosomia",
    ">=5.0 kg very extreme"
  ),
  Threshold_kg = c(NA_real_, MACROSOMIA_THRESHOLD_KG, EXTREME_HIGH_BW_THRESHOLD_KG, VERY_EXTREME_BW_THRESHOLD_KG),
  N = c(
    nrow(cleaned_data),
    sum(cleaned_data$Macrosomia_4kg_Flag, na.rm = TRUE),
    sum(cleaned_data$Extreme_High_BW_4p5kg_Flag, na.rm = TRUE),
    sum(cleaned_data$Very_Extreme_BW_5kg_Flag, na.rm = TRUE)
  ),
  Percent = 100 * c(
    1,
    mean(cleaned_data$Macrosomia_4kg_Flag, na.rm = TRUE),
    mean(cleaned_data$Extreme_High_BW_4p5kg_Flag, na.rm = TRUE),
    mean(cleaned_data$Very_Extreme_BW_5kg_Flag, na.rm = TRUE)
  ),
  Basis = c(
    "Full analytic cohort",
    "Common absolute macrosomia threshold",
    "Prespecified sensitivity threshold; ACOG PB216 risk category begins at 4.5 kg",
    "ACOG PB216 highest displayed category begins at 5.0 kg"
  )
)

bw_12mo_primary <- cleaned_data %>%
  rowwise() %>%
  mutate(
    BW_12mo_kg = get_infant_phys_meta(
      age_yr = 1,
      gender_code = infsex,
      BW_birth_kg = BW_birth_kg
    )$BW
  ) %>%
  ungroup() %>%
  select(
    ID, infsex, BW_birth_kg, wt_ratio_legacy,
    Macrosomia_4kg_Flag, Extreme_High_BW_4p5kg_Flag, Very_Extreme_BW_5kg_Flag,
    BW_12mo_kg
  )

bw_growth_sensitivity <- dplyr::bind_rows(lapply(
  BW_REMAINING_FRACTION_SENSITIVITY,
  function(rem) {
    cleaned_data %>%
      rowwise() %>%
      mutate(
        Remaining_Fraction_12mo = rem,
        BW_12mo_kg = get_infant_phys_meta(
          age_yr = 1,
          gender_code = infsex,
          BW_birth_kg = BW_birth_kg,
          remaining_at_12mo = rem,
          mode = ifelse(rem == 1, "CONSTANT_RATIO", "DECAY_TO_REFERENCE")
        )$BW
      ) %>%
      ungroup() %>%
      select(ID, infsex, BW_birth_kg, Extreme_High_BW_4p5kg_Flag, Remaining_Fraction_12mo, BW_12mo_kg)
  }
))

bw_growth_sensitivity_summary <- bw_growth_sensitivity %>%
  group_by(Remaining_Fraction_12mo) %>%
  summarise(
    N = n(),
    Median_BW_12mo_kg = median(BW_12mo_kg, na.rm = TRUE),
    P95_BW_12mo_kg = quantile(BW_12mo_kg, 0.95, na.rm = TRUE, names = FALSE),
    Max_BW_12mo_kg = max(BW_12mo_kg, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf(
  "Birth-weight flags: >=4.0 kg: %d/%d; >=4.5 kg: %d/%d; >=5.0 kg: %d/%d\n",
  sum(cleaned_data$Macrosomia_4kg_Flag), nrow(cleaned_data),
  sum(cleaned_data$Extreme_High_BW_4p5kg_Flag), nrow(cleaned_data),
  sum(cleaned_data$Very_Extreme_BW_5kg_Flag), nrow(cleaned_data)
))
cat(sprintf(
  "Primary growth model 12-mo BW: median %.3f kg; 95th percentile %.3f kg; max %.3f kg\n",
  median(bw_12mo_primary$BW_12mo_kg),
  quantile(bw_12mo_primary$BW_12mo_kg, 0.95, names = FALSE),
  max(bw_12mo_primary$BW_12mo_kg)
))

write_xlsx(
  list(
    Cleaned_Input = cleaned_data,
    Birth_Weight_Classification = bw_class_summary,
    Primary_12mo_Growth_QC = bw_12mo_primary,
    Growth_Assumption_Sensitivity = bw_growth_sensitivity,
    Growth_Sensitivity_Summary = bw_growth_sensitivity_summary
  ),
  path = file.path(OUTPUT_DIR, "tables", "cleaned_input_and_growth_QC_v3.xlsx")
)

# ============================================================
# 9. MILK FUNCTIONS
# ============================================================
make_milk_function <- function(p_name, c_milk0, scenario) {
  if (scenario == "S0") {
    return(function(t_days) rep(c_milk0, length(t_days)))
  }
  
  if (scenario == "S1") {
    p <- trajectory_params[
      match(p_name, trajectory_params$Compound),
      , drop = FALSE
    ]
    
    if (nrow(p) != 1) {
      stop("No trajectory parameter row for ", p_name)
    }
    
    if (!isTRUE(p$supported_by_blomberg)) {
      return(function(t_days) rep(c_milk0, length(t_days)))
    }
    
    beta <- as.numeric(p$beta_mean_day)
    return(function(t_days) c_milk0 * exp(beta * t_days))
  }
  
  stop("Unknown scenario: ", scenario)
}

# ============================================================
# 10. RENAL PARAMETER SCALING
# ============================================================
get_tmkt_row <- function(compound) {
  p <- TMKT_PARAMS[
    TMKT_PARAMS$Compound == compound,
    , drop = FALSE
  ]
  if (nrow(p) != 1) {
    stop("No unique Tm/Kt parameter row for ", compound)
  }
  p
}

scale_tmc_to_tm <- function(tmc_ug_h_bw075, BW_kg) {
  BW_safe <- pmax(as.numeric(BW_kg), 0.1)
  as.numeric(tmc_ug_h_bw075) * BW_safe^TM_SCALING_EXPONENT * 1000
}

# EFSA/PFAS template allometric form:
#   ku = KuC * BW^(-0.25)
scale_ku <- function(ku_ref, BW_kg) {
  if (!is.finite(ku_ref)) return(rep(NA_real_, length(BW_kg)))
  BW_safe <- pmax(as.numeric(BW_kg), 0.1)
  as.numeric(ku_ref) * BW_safe^KU_SCALING_EXPONENT
}

# Filtrate -> urinary-storage transfer rate.
# This is NOT the chemical-specific Ku parameter.
scale_kust <- function(kust_ref, BW_kg) {
  BW_safe <- pmax(as.numeric(BW_kg), 0.1)
  as.numeric(kust_ref) * BW_safe^KUST_SCALING_EXPONENT
}

# ============================================================
# 11. TIME GRIDS
# ============================================================
times <- seq(
  0,
  SIM_YEARS * 365.25 * 24,
  by = TIME_STEP_HOURS
)

age_seq_yr <- times / 8760

renal_age_seq_yr <- seq(
  0,
  SIM_YEARS,
  by = RENAL_CACHE_STEP_DAYS / 365.25
)

if (tail(renal_age_seq_yr, 1) < SIM_YEARS) {
  renal_age_seq_yr <- c(renal_age_seq_yr, SIM_YEARS)
}

# ============================================================
# 12. INTERPOLATION HELPER
# ============================================================
safe_approxfun <- function(x, y, default = 0) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 2L) {
    val <- if (sum(ok) == 1L) y[which(ok)[1]] else default
    return(function(z) rep(val, length(z)))
  }
  
  stats::approxfun(x[ok], y[ok], rule = 2)
}

# ============================================================
# 13. SUBJECT PHYSIOLOGY CACHE
# ============================================================
make_subject_physiology_cache <- function(sub) {
  phys_df <- dplyr::bind_rows(lapply(renal_age_seq_yr, function(a) {
    get_infant_phys_meta(
      age_yr = a,
      gender_code = sub$infsex,
      BW_birth_kg = sub$BW_birth_kg
    )
  }))
  
  list(
    phys_df = phys_df,
    fns = list(
      f_BW = safe_approxfun(renal_age_seq_yr, phys_df$BW),
      f_V_pl = safe_approxfun(renal_age_seq_yr, phys_df$V_pl),
      f_V_li = safe_approxfun(renal_age_seq_yr, phys_df$V_li),
      f_V_kd = safe_approxfun(renal_age_seq_yr, phys_df$V_kd),
      f_V_re = safe_approxfun(renal_age_seq_yr, phys_df$V_re),
      f_Q_li = safe_approxfun(renal_age_seq_yr, phys_df$Q_li_mL_min),
      f_Q_kd = safe_approxfun(renal_age_seq_yr, phys_df$Q_kd_mL_min),
      f_Q_re = safe_approxfun(renal_age_seq_yr, phys_df$Q_re_mL_min),
      f_GFR = safe_approxfun(renal_age_seq_yr, phys_df$GFR_mL_min),
      f_KMilk = safe_approxfun(renal_age_seq_yr, phys_df$KMilk_L_h)
    )
  )
}

# ============================================================
# 14. SUBJECT × COMPOUND RENAL CACHE
# ============================================================
make_renal_cache_unified <- function(compound, phys_df) {
  p <- get_tmkt_row(compound)
  age <- renal_age_seq_yr
  
  BW <- pmax(as.numeric(phys_df$BW), 0.1)
  
  GFR_L_h <- pmax(
    as.numeric(phys_df$GFR_mL_min),
    0
  ) * 60 / 1000
  
  # Qfil is the glomerular filtration flow used by the mechanistic renal module.
  Qfil_L_h <- GFR_L_h
  
  Vfil_L <- pmax(
    VFIL_L_PER_KG_BW * BW,
    MIN_VOLUME_L
  )
  
  Tm_ng_h <- scale_tmc_to_tm(
    p$Tmc_ug_h_bw075,
    BW
  )
  
  Kt_ng_L <- rep(
    as.numeric(p$Kt_ug_L) * 1000,
    length(BW)
  )
  
  Free <- rep(
    as.numeric(p$Free_fraction),
    length(BW)
  )
  
  # k_ust: filtrate -> urinary storage.
  # We use Qfil/Vfil as a residence-rate approximation because this is a
  # flow/storage parameter, not a chemical-specific renal clearance.
  if (KUST_MODE == "QFIL_OVER_VFIL") {
    kust_h <- ifelse(Vfil_L > 0, Qfil_L_h / Vfil_L, 0)
    kust_source <- "Physiological residence-rate approximation: Qfil/Vfil"
  } else {
    stop("Unsupported KUST_MODE = ", KUST_MODE)
  }
  
  # ku: urinary storage -> cumulative urine.
  ku_h <- scale_ku(
    p$KurineC_h_bwminus025,
    BW
  )
  
  list(
    model_type = "Unified_Tm_Kt_with_Urinary_Storage",
    source_row = p,
    f_Free = safe_approxfun(age, Free),
    f_Tm = safe_approxfun(age, Tm_ng_h),
    f_Kt = safe_approxfun(age, Kt_ng_L),
    f_Qfil = safe_approxfun(age, Qfil_L_h),
    f_kust = safe_approxfun(age, kust_h),
    f_ku = safe_approxfun(age, ku_h, default = NA_real_),
    f_Vfil = safe_approxfun(age, Vfil_L),
    f_GFR = safe_approxfun(age, GFR_L_h),
    metadata = list(
      Tmc_ug_h_bw075 = p$Tmc_ug_h_bw075,
      Kt_ug_L = p$Kt_ug_L,
      Free = p$Free_fraction,
      KurineC_h_bwminus025 = p$KurineC_h_bwminus025,
      Kurine_source = p$Kurine_source,
      Kust_source = kust_source,
      Tm_basis = p$Tm_basis,
      Tm_source = p$Tm_source,
      Parameter_Class = p$Parameter_class,
      Ku_missing = !is.finite(p$KurineC_h_bwminus025)
    )
  )
}

# ============================================================
# 15. UNIFIED RENAL PBPK ODE
# ============================================================
# States:
#   Al       = liver amount (ng)
#   Ak       = kidney-tissue amount (ng)
#   Ar       = rest/deep amount (ng)
#   Ap       = plasma amount (ng)
#   Afil     = renal filtrate amount (ng)
#   Aust     = urinary-storage amount (ng)
#   Aurine   = cumulative urinary excretion (ng)
#   Afeces   = cumulative fecal excretion (ng; zero in current base model)
#   Ainput    = cumulative postnatal breast-milk input (ng), integrated by LSODA
#
# Renal equations:
#   R_filtration = Qfil * Free * Cp
#   R_reabsorption = Tm * Cfil / (Kt + Cfil)
#   R_ust = kust * Afil
#   R_urine = ku * Aust
#
#   dAfil = R_filtration - R_reabsorption - R_ust
#   dAust = R_ust - R_urine
#   dAurine = R_urine
#
# Reabsorbed material is returned directly to plasma, representing the
# net post-tubular re-entry into systemic circulation.
# ============================================================
pbpk_ode_unified_tmkt <- function(t, state, parameters) {
  with(as.list(parameters), {
    age_yr <- t / 8760
    
    V_pl <- max(f_V_pl(age_yr), MIN_VOLUME_L)
    V_li <- max(f_V_li(age_yr), MIN_VOLUME_L)
    V_kd <- max(f_V_kd(age_yr), MIN_VOLUME_L)
    V_re <- max(f_V_re(age_yr), MIN_VOLUME_L)
    V_fil <- max(f_V_fil(age_yr), MIN_VOLUME_L)
    
    Q_li <- max(f_Q_li(age_yr), 0) * 60 / 1000
    Q_kd <- max(f_Q_kd(age_yr), 0) * 60 / 1000
    Q_re <- max(f_Q_re(age_yr), 0) * 60 / 1000
    
    Qfil <- max(f_Qfil(age_yr), 0)
    
    KMilk <- max(f_KMilk(age_yr), 0)
    
    Free <- max(f_Free(age_yr), 0)
    Tm <- max(f_Tm(age_yr), 0)
    Kt <- max(f_Kt(age_yr), 1e-18)
    kust <- max(f_kust(age_yr), 0)
    ku <- f_ku(age_yr)
    
    if (!is.finite(ku) || ku < 0) {
      if (ALLOW_MISSING_KU) {
        ku <- 0
      } else {
        stop("Missing/invalid chemical-specific Ku for this compound.")
      }
    }
    
    Al <- state[1]
    Ak <- state[2]
    Ar <- state[3]
    Ap <- state[4]
    Afil <- state[5]
    Aust <- state[6]
    
    Cp <- max(Ap / V_pl, 0)
    Cl <- max(Al / V_li, 0)
    Ck <- max(Ak / V_kd, 0)
    Cr <- max(Ar / V_re, 0)
    Cfil <- max(Afil / V_fil, 0)
    
    C_milk <- f_C_milk(t / 24)
    Ingestion <- KMilk * C_milk * mix_ratio
    
    # ----------------------------------------------------------
    # Renal processes
    # ----------------------------------------------------------
    R_filtration <- Qfil * Free * Cp
    
    R_reabsorption <- ifelse(
      Cfil > 0,
      Tm * Cfil / (Kt + Cfil),
      0
    )
    
    R_ust <- kust * Afil
    R_urine <- ku * Aust
    
    # ----------------------------------------------------------
    # Tissue exchange
    # ----------------------------------------------------------
    dAl <- Q_li * (Cp - Cl / P_li)
    
    dAk <- Q_kd * (Free * Cp - Ck / P_kd)
    
    dAr <- Q_re * (Cp - Cr / P_re)
    
    # ----------------------------------------------------------
    # Plasma mass balance
    # ----------------------------------------------------------
    dAp <-
      Q_li * Cl / P_li +
      Q_kd * Ck / P_kd +
      Q_re * Cr / P_re -
      (Q_li + Q_re) * Cp -
      Q_kd * Free * Cp -
      R_filtration +
      R_reabsorption +
      Ingestion
    
    # ----------------------------------------------------------
    # Renal storage compartments
    # ----------------------------------------------------------
    dAfil <- R_filtration - R_reabsorption - R_ust
    dAust <- R_ust - R_urine
    
    # Cumulative excretion
    dAurine <- R_urine
    dAfeces <- 0
    
    # Exact cumulative postnatal input.  Integrating this as an auxiliary
    # ODE state ensures that mass-balance bookkeeping uses the same adaptive
    # LSODA integration as the PBPK equations rather than a 7-day trapezoid.
    dAinput <- Ingestion
    
    list(c(
      dAl,
      dAk,
      dAr,
      dAp,
      dAfil,
      dAust,
      dAurine,
      dAfeces,
      dAinput
    ))
  })
}

# ============================================================
# 16. SINGLE SUBJECT / COMPOUND / SCENARIO
# ============================================================
run_single_person_unified <- function(
    sub, p_name, scenario, subject_phys, renal_cache
) {
  
  c_milk0 <- as.numeric(sub[[p_name]])
  if (!is.finite(c_milk0) || c_milk0 <= 0) return(NULL)
  
  mix_ratio <- as.numeric(sub$mix_ratio)
  if (!is.finite(mix_ratio) || mix_ratio < 0 || mix_ratio > 1) {
    mix_ratio <- 1
  }
  
  f_BW <- subject_phys$fns$f_BW
  f_V_pl <- subject_phys$fns$f_V_pl
  f_V_li <- subject_phys$fns$f_V_li
  f_V_kd <- subject_phys$fns$f_V_kd
  f_V_re <- subject_phys$fns$f_V_re
  f_Q_li <- subject_phys$fns$f_Q_li
  f_Q_kd <- subject_phys$fns$f_Q_kd
  f_Q_re <- subject_phys$fns$f_Q_re
  f_KMilk <- subject_phys$fns$f_KMilk
  
  f_C_milk <- make_milk_function(p_name, c_milk0, scenario)
  
  pars <- list(
    P_li = as.numeric(pfas_params[[p_name]]["P_li"]),
    P_kd = as.numeric(pfas_params[[p_name]]["P_kd"]),
    P_re = as.numeric(pfas_params[[p_name]]["P_re"]),
    mix_ratio = mix_ratio,
    
    f_V_pl = f_V_pl,
    f_V_li = f_V_li,
    f_V_kd = f_V_kd,
    f_V_re = f_V_re,
    f_Q_li = f_Q_li,
    f_Q_kd = f_Q_kd,
    f_Q_re = f_Q_re,
    f_KMilk = f_KMilk,
    
    f_Free = renal_cache$f_Free,
    f_Tm = renal_cache$f_Tm,
    f_Kt = renal_cache$f_Kt,
    f_Qfil = renal_cache$f_Qfil,
    f_kust = renal_cache$f_kust,
    f_ku = renal_cache$f_ku,
    f_V_fil = renal_cache$f_Vfil,
    
    f_C_milk = function(t_days) f_C_milk(t_days)
  )
  
  init4 <- get_infant_initial_burden(
    compound_name = p_name,
    C_milk_measured = c_milk0,
    V_pl_0 = f_V_pl(0),
    V_li_0 = f_V_li(0),
    V_kd_0 = f_V_kd(0),
    V_re_0 = f_V_re(0),
    P_li = pars$P_li,
    P_kd = pars$P_kd,
    P_re = pars$P_re
  )
  
  init_full <- c(
    Al = init4[["Al"]],
    Ak = init4[["Ak"]],
    Ar = init4[["Ar"]],
    Ap = init4[["Ap"]],
    Afil = 0,
    Aust = 0,
    Aurine = 0,
    Afeces = 0,
    Ainput = 0
  )
  
  run_ode <- function(init, milk_fun, times_run = times) {
    pars_run <- pars
    pars_run$f_C_milk <- function(t_days) milk_fun(t_days)
    
    deSolve::ode(
      y = init,
      times = times_run,
      func = pbpk_ode_unified_tmkt,
      parms = pars_run,
      method = "lsoda",
      rtol = ODE_RTOL,
      atol = ODE_ATOL
    )
  }
  
  out_full <- tryCatch(
    run_ode(init_full, f_C_milk, times),
    error = function(e) {
      warning(
        "ODE failed: ID=", sub$ID,
        ", compound=", p_name,
        ", scenario=", scenario,
        ": ", e$message
      )
      NULL
    }
  )
  
  if (is.null(out_full)) return(NULL)
  
  out_df <- as.data.frame(out_full)
  names(out_df)[1:10] <- c(
    "time", "Al", "Ak", "Ar", "Ap", "Afil", "Aust", "Aurine", "Afeces", "Ainput"
  )
  
  t_yr <- out_df$time / 8760
  t_days <- out_df$time / 24
  
  # ----------------------------------------------------------
  # NOTE ON PRENATAL / POSTNATAL ATTRIBUTION
  # ----------------------------------------------------------
  # The current nonlinear renal model does NOT satisfy linear
  # superposition. The following values are therefore retained only as
  # approximate counterfactual ratios for backward compatibility.
  # Final publication analysis should use explicit TOTAL / PRENATAL_ONLY /
  # BREASTFEEDING_ONLY counterfactual simulations.
  # ----------------------------------------------------------
  contribution_target_h <- min(
    SIM_YEARS * 365.25 * 24,
    CONTRIBUTION_MONTHS / 12 * 365.25 * 24
  )
  contribution_times <- c(0, contribution_target_h)
  
  pre_out <- tryCatch(
    run_ode(
      init_full,
      function(t_days) rep(0, length(t_days)),
      contribution_times
    ),
    error = function(e) NULL
  )
  
  post_init <- init_full
  post_init[] <- 0
  
  post_out <- tryCatch(
    run_ode(
      post_init,
      f_C_milk,
      contribution_times
    ),
    error = function(e) NULL
  )
  
  pre_contrib_approx <- NA_real_
  post_contrib_approx <- NA_real_
  
  if (!is.null(pre_out) && !is.null(post_out)) {
    pre_df <- as.data.frame(pre_out)
    post_df <- as.data.frame(post_out)
    
    names(pre_df)[1:10] <- names(out_df)[1:10]
    names(post_df)[1:10] <- names(out_df)[1:10]
    
    idx12 <- which.min(abs(pre_df$time - contribution_target_h))
    
    pre_C12 <- pre_df$Ap[idx12] /
      max(f_V_pl(pre_df$time[idx12] / 8760), MIN_VOLUME_L)
    
    post_C12 <- post_df$Ap[idx12] /
      max(f_V_pl(post_df$time[idx12] / 8760), MIN_VOLUME_L)
    
    if (
      is.finite(pre_C12) &&
      is.finite(post_C12) &&
      pre_C12 + post_C12 > 0
    ) {
      pre_contrib_approx <- pre_C12 / (pre_C12 + post_C12)
      post_contrib_approx <- post_C12 / (pre_C12 + post_C12)
    }
  }
  
  # ----------------------------------------------------------
  # Vectorized diagnostics
  # ----------------------------------------------------------
  BW <- f_BW(t_yr)
  Vpl <- pmax(f_V_pl(t_yr), MIN_VOLUME_L)
  Vli <- pmax(f_V_li(t_yr), MIN_VOLUME_L)
  Vkd <- pmax(f_V_kd(t_yr), MIN_VOLUME_L)
  Vre <- pmax(f_V_re(t_yr), MIN_VOLUME_L)
  Vfil <- pmax(renal_cache$f_Vfil(t_yr), MIN_VOLUME_L)
  
  Qfil_L_h <- pmax(renal_cache$f_Qfil(t_yr), 0)
  Free <- pmax(renal_cache$f_Free(t_yr), 0)
  Tm <- pmax(renal_cache$f_Tm(t_yr), 0)
  Kt <- pmax(renal_cache$f_Kt(t_yr), 1e-18)
  Kust <- pmax(renal_cache$f_kust(t_yr), 0)
  Ku <- renal_cache$f_ku(t_yr)
  Ku_valid <- is.finite(Ku) & Ku >= 0
  Ku[!Ku_valid] <- NA_real_
  
  C_Plasma <- out_df$Ap / Vpl
  C_Liver <- out_df$Al / Vli
  C_Kidney <- out_df$Ak / Vkd
  C_Rest <- out_df$Ar / Vre
  C_Filtrate <- out_df$Afil / Vfil
  
  Milk_C <- f_C_milk(t_days)
  
  R_Filtration <- Qfil_L_h * Free * pmax(C_Plasma, 0)
  
  R_Reabsorption <- ifelse(
    C_Filtrate > 0,
    Tm * C_Filtrate / (Kt + C_Filtrate),
    0
  )
  
  R_Urine_Storage_In <- Kust * pmax(out_df$Afil, 0)
  R_Urine_Excretion <- ifelse(
    Ku_valid,
    Ku * pmax(out_df$Aust, 0),
    NA_real_
  )
  
  Reabsorption_Rate_Ratio <- ifelse(
    R_Filtration > 0,
    R_Reabsorption / R_Filtration,
    NA_real_
  )
  
  Urinary_Clearance_Rate_Ratio <- ifelse(
    R_Filtration > 0 & is.finite(R_Urine_Excretion),
    R_Urine_Excretion / R_Filtration,
    NA_real_
  )
  
  # ----------------------------------------------------------
  # Mass balance
  # ----------------------------------------------------------
  # A_Total_Model includes cumulative urine/feces. Therefore excretion must
  # NOT be subtracted again from the expected total. The exact identity is:
  #
  #   Al + Ak + Ar + Ap + Afil + Aust + Aurine + Afeces
  #       = initial_total + cumulative_external_input
  #
  # Cumulative input is integrated inside the ODE as state Ainput.
  A_Total_Model <-
    out_df$Al +
    out_df$Ak +
    out_df$Ar +
    out_df$Ap +
    out_df$Afil +
    out_df$Aust +
    out_df$Aurine +
    out_df$Afeces
  
  initial_total <-
    out_df$Al[1] +
    out_df$Ak[1] +
    out_df$Ar[1] +
    out_df$Ap[1] +
    out_df$Afil[1] +
    out_df$Aust[1] +
    out_df$Aurine[1] +
    out_df$Afeces[1]
  
  Ingestion_Rate <- f_KMilk(t_yr) * Milk_C * mix_ratio
  
  # Exact LSODA-integrated cumulative external input.
  A_Input_Cum <- out_df$Ainput
  
  # Retained as a separate diagnostic only. It is already part of
  # A_Total_Model and must not be subtracted from Mass_Balance_Expected.
  A_Excretion_Cum <- out_df$Aurine + out_df$Afeces
  
  Mass_Balance_Expected <- initial_total + A_Input_Cum
  
  Mass_Balance_Error <- A_Total_Model - Mass_Balance_Expected
  
  Mass_Balance_Relative_Error <-
    Mass_Balance_Error /
    pmax(abs(Mass_Balance_Expected), 1e-18)
  
  tibble::tibble(
    Time_yr = t_yr,
    Time_days = t_days,
    ID = sub$ID,
    Compound = p_name,
    Scenario = scenario,
    Renal_Model = "Unified_Tm_Kt_UrinaryStorage",
    
    BW_current_kg = BW,
    BW_birth_kg = sub$BW_birth_kg,
    Macrosomia_4kg_Flag = sub$Macrosomia_4kg_Flag,
    Extreme_High_BW_4p5kg_Flag = sub$Extreme_High_BW_4p5kg_Flag,
    Very_Extreme_BW_5kg_Flag = sub$Very_Extreme_BW_5kg_Flag,
    Milk_C = Milk_C,
    
    C_Plasma = C_Plasma,
    C_Liver = C_Liver,
    C_Kidney = C_Kidney,
    C_Rest = C_Rest,
    C_Filtrate = C_Filtrate,
    
    A_Urine_Storage = out_df$Aust,
    A_Urine = out_df$Aurine,
    A_Feces = out_df$Afeces,
    
    A_Total_Model = A_Total_Model,
    # Classical tissue/plasma body burden retained for backward compatibility.
    A_Body_Burden =
      out_df$Al + out_df$Ak + out_df$Ar + out_df$Ap,
    # Systemically accessible burden includes renal filtrate but excludes Aust,
    # which is downstream of the return-to-plasma pathway.
    A_Systemic_Burden =
      out_df$Al + out_df$Ak + out_df$Ar + out_df$Ap + out_df$Afil,
    
    Qfil_L_h = Qfil_L_h,
    R_Filtration = R_Filtration,
    R_Reabsorption = R_Reabsorption,
    R_Urine_Storage_In = R_Urine_Storage_In,
    R_Urine_Excretion = R_Urine_Excretion,
    
    Reabsorption_Rate_Ratio = Reabsorption_Rate_Ratio,
    Urinary_Clearance_Rate_Ratio = Urinary_Clearance_Rate_Ratio,
    
    Tm_ng_h = Tm,
    Tm_ug_h = Tm / 1000,
    Kt_ng_L = Kt,
    Kt_ug_L = Kt / 1000,
    Free_fraction = Free,
    
    Kust_h = Kust,
    Ku_h = Ku,
    Vfil_L = Vfil,
    GFR_L_h = Qfil_L_h,
    
    Ingestion_Rate_ng_h = Ingestion_Rate,
    Ainput_ODE = out_df$Ainput,
    A_Input_Cum = A_Input_Cum,
    A_Excretion_Cum = A_Excretion_Cum,
    
    Mass_Balance_Expected = Mass_Balance_Expected,
    Mass_Balance_Error = Mass_Balance_Error,
    Mass_Balance_Relative_Error = Mass_Balance_Relative_Error,
    
    Pre_Contribution_Approx_12mo = pre_contrib_approx,
    Post_Contribution_Approx_12mo = post_contrib_approx,
    
    Tm_Source = renal_cache$metadata$Tm_source,
    Parameter_Class = renal_cache$metadata$Parameter_Class,
    Kust_Source = renal_cache$metadata$Kust_source,
    Ku_Source = renal_cache$metadata$Kurine_source,
    Ku_Missing = renal_cache$metadata$Ku_missing
  )
}

# ============================================================
# 17. SUBJECT-LEVEL PARALLEL SIMULATION
# ============================================================
run_subject_all_unified <- function(job) {
  sub <- job$sub
  subject_phys <- job$phys
  
  out <- list()
  k <- 1L
  
  for (p_name in pfas_cols) {
    renal_cache <- make_renal_cache_unified(
      p_name,
      subject_phys$phys_df
    )
    
    for (scenario in SCENARIOS) {
      res <- run_single_person_unified(
        sub = sub,
        p_name = p_name,
        scenario = scenario,
        subject_phys = subject_phys,
        renal_cache = renal_cache
      )
      
      if (!is.null(res)) {
        out[[k]] <- res
        k <- k + 1L
      }
    }
  }
  
  if (!length(out)) return(tibble::tibble())
  dplyr::bind_rows(out)
}

# ============================================================

# ============================================================

# ============================================================
# PFAS_RCM_RBM_source_attribution_sensitivity_v1.R
#
# Joint R_CM / R_BM sensitivity of formal two-source Shapley attribution
# plus shared-C0 scaling diagnostic.
#
# Purpose
# -------
# Directly test the robustness of the manuscript's early prenatal-to-postnatal
# source transition to the transfer assumptions used to reconstruct birth burden.
#
# The analysis uses the FINAL FREEZE v1.1 PBPK structure above without changing
# tissue, renal, physiology, solver, milk-trajectory, or cohort parameters.
#
# Primary source framework
# ------------------------
# TOTAL:
#   reconstructed birth burden + postnatal breast-milk input
#
# PRENATAL_ONLY:
#   reconstructed birth burden + zero postnatal input
#
# POSTNATAL_ONLY:
#   zero initial burden + postnatal breast-milk input
#
# For each time t:
#   I(t)          = C_TOTAL(t) - C_PRE(t) - C_POST(t)
#   C_PRE_SH(t)   = C_PRE(t)  + I(t)/2
#   C_POST_SH(t)  = C_POST(t) + I(t)/2
#   F_POST(t)     = C_POST_SH(t) / C_TOTAL(t)
#
# The zero-initial, zero-input baseline is zero, so this is the two-source
# Shapley allocation used in the manuscript.
#
# IMPORTANT STRUCTURAL POINT
# --------------------------
# In the Shanghai primary cohort:
#
#   C_fetal = C_milk * R_CM / R_BM
#
# while measured milk itself drives postnatal intake. Therefore R_CM and R_BM
# affect prenatal reconstruction through their R_CM/R_BM ratio. Grid cells with
# the same multiplier ratio RCM_mult/RBM_mult are algebraically equivalent for
# the prenatal initialization. The script deliberately reports the full joint
# grid but computes each unique ratio only once.
#
# The default R_CM/R_BM multiplier grid is a PRESPECIFIED FACTOR SENSITIVITY,
# not a confidence interval and not a claim that every factor is equally likely.
# Replace the factor grid with compound-specific literature bounds if a harmonized
# evidence table becomes available.
#
# The C0 diagnostic scales the measured colostrum concentration simultaneously
# in BOTH pathways, because the actual model uses the same measured C0 to
# reconstruct birth burden and to drive postnatal input. This specifically tests
# whether source fractions are structurally insensitive to shared exposure scale.
#
# Outputs
# -------
# 1) Individual and compound-level crossover ages
# 2) Formal Shapley postnatal shares at 1, 3, 6, and 12 months
# 3) First crossing and 14-day sustained crossing
# 4) Full R_CM x R_BM joint-grid summaries
# 5) Shared-C0 scaling diagnostic
# 6) Baseline daily Shapley trajectories for final source-attribution plotting
# 7) Numerical and attribution QC
#
# No source conclusion should be interpreted as total real-world exposure:
# only prenatal reconstructed burden and postnatal breast-milk input are included.
# ============================================================

SENSITIVITY_VERSION <- "PFAS_RCM_RBM_source_attribution_sensitivity_v1_EP"

SENSITIVITY_OUTPUT_DIR <- file.path(
  PROJECT_ROOT, "results", "results_RCM_RBM_source_attribution_sensitivity_v1_EP_FAST3"
)

dir.create(SENSITIVITY_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(SENSITIVITY_OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(SENSITIVITY_OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(SENSITIVITY_OUTPUT_DIR, "qc"), recursive = TRUE, showWarnings = FALSE)

# Attribution is evaluated for the literature-informed milk trajectory.
ATTR_SCENARIO <- "S1"

# Daily reporting is used here specifically to avoid over-interpreting the
# 7-day main-model reporting grid when deriving crossover age.
ATTR_STEP_DAYS <- 1
ATTR_END_DAYS <- 365.25

ATTR_DAYS <- unique(sort(c(
  seq(0, floor(ATTR_END_DAYS), by = ATTR_STEP_DAYS),
  ATTR_END_DAYS
)))
ATTR_TIMES_H <- ATTR_DAYS * 24

# Endpoint ages used for source-share reporting.
ENDPOINT_MONTHS <- c(1, 3, 6, 12)
ENDPOINT_DAYS <- ENDPOINT_MONTHS * 365.25 / 12

# First sustained postnatal dominance must persist for this many consecutive days.
SUSTAIN_DAYS <- 14L

# Joint transfer-ratio sensitivity factors.
# These are sensitivity factors, NOT statistical confidence limits.
RCM_MULTIPLIERS <- c(0.50, 0.75, 1.00, 1.25, 1.50)
RBM_MULTIPLIERS <- c(0.50, 0.75, 1.00, 1.25, 1.50)

# Shared measured-colostrum scaling diagnostic.
C0_MULTIPLIERS <- c(0.50, 0.75, 1.00, 1.25, 1.50, 2.00)

# Numerical thresholds.
ATTR_EPS <- 1e-12
SHARE_QC_TOL <- 1e-8
MASS_BALANCE_QC_TOL <- 1e-5

# Save baseline daily trajectories for plotting / Supplement.
SAVE_BASELINE_DAILY_TRAJECTORIES <- TRUE

# Parallel execution.
ATTR_USE_PARALLEL <- USE_PARALLEL
# Windows multisession processes can compete for RAM at very high worker counts.
# Raise this cap if sufficient memory is available; output values do not change.
ATTR_WORKERS <- min(PARALLEL_WORKERS, 6L)

# Enable direct-binding RHS: exactly the same ODE expression without creating
# as.list(parameters) at every LSODA derivative call.
ATTR_FAST_DIRECT_RHS <- TRUE

# ============================================================
# A1. HELPERS
# ============================================================

safe_median_attr <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_q_attr <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probs = p, names = FALSE, na.rm = TRUE))
}

interp_attr <- function(x, y, xout) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L) return(rep(NA_real_, length(xout)))
  approx(
    x = x[keep],
    y = y[keep],
    xout = xout,
    rule = 2,
    ties = mean
  )$y
}

first_crossing_day <- function(days, fpost, threshold = 0.5) {
  ok <- is.finite(days) & is.finite(fpost)
  days <- days[ok]
  fpost <- fpost[ok]
  if (!length(days)) return(NA_real_)

  idx <- which(fpost > threshold)
  if (!length(idx)) return(NA_real_)

  i <- idx[1]

  # Linear interpolation around the first observed threshold crossing.
  if (i > 1L &&
      is.finite(fpost[i - 1L]) &&
      fpost[i - 1L] <= threshold &&
      fpost[i] != fpost[i - 1L]) {

    x0 <- days[i - 1L]
    x1 <- days[i]
    y0 <- fpost[i - 1L]
    y1 <- fpost[i]

    return(
      x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)
    )
  }

  days[i]
}

first_sustained_crossing_day <- function(
    days, fpost, threshold = 0.5, sustain_days = SUSTAIN_DAYS
) {
  ok <- is.finite(days) & is.finite(fpost)
  days <- days[ok]
  fpost <- fpost[ok]
  if (length(days) < sustain_days + 1L) return(NA_real_)

  # Exact equivalent of the original search for the FIRST contiguous
  # (sustain_days + 1)-point above-threshold window. No smoothing.
  n_req <- sustain_days + 1L
  above <- fpost > threshold
  # Prefix sum yields O(n) instead of repeated all(above[i:j]) evaluation.
  prefix <- c(0L, cumsum(as.integer(above)))
  starts <- seq_len(length(above) - n_req + 1L)
  ends <- starts + n_req - 1L
  candidates <- starts[
    (prefix[ends + 1L] - prefix[starts] == n_req) &
      ((days[ends] - days[starts]) >= sustain_days - 1e-8)
  ]
  if (!length(candidates)) return(NA_real_)
  i <- candidates[1L]
  j <- i + n_req - 1L
  ix <- max(1L, i - 1L):j
  first_crossing_day(days[ix], fpost[ix], threshold)
}

# ============================================================
# A2. TRANSFER GRID
# ============================================================

transfer_baseline <- transfer_data %>%
  select(
    Compound,
    R_CM_Baseline = R_CM,
    R_BM_Baseline = R_BM,
    R_CM_Source,
    R_BM_Source,
    R_BM_Evidence_Class
  )

joint_grid <- tidyr::expand_grid(
  RCM_Multiplier = RCM_MULTIPLIERS,
  RBM_Multiplier = RBM_MULTIPLIERS
) %>%
  mutate(
    Prenatal_Ratio_Multiplier =
      RCM_Multiplier /
      RBM_Multiplier,

    Grid_ID =
      paste0(
        "RCM_", format(RCM_Multiplier, trim = TRUE),
        "__RBM_", format(RBM_Multiplier, trim = TRUE)
      )
  )

# QC: every R_CM x R_BM design point must appear exactly once.
stopifnot(
  nrow(joint_grid) == length(RCM_MULTIPLIERS) * length(RBM_MULTIPLIERS),
  !anyDuplicated(joint_grid[c("RCM_Multiplier", "RBM_Multiplier")]),
  all(is.finite(joint_grid$Prenatal_Ratio_Multiplier))
)

# Compute only unique ratio multipliers because R_CM and R_BM enter
# the primary cohort initialization through R_CM/R_BM.
unique_ratio_grid <- joint_grid %>%
  distinct(
    Prenatal_Ratio_Multiplier
  ) %>%
  arrange(
    Prenatal_Ratio_Multiplier
  ) %>%
  mutate(
    Ratio_ID =
      sprintf(
        "RATIO_%02d",
        row_number()
      )
  )

joint_grid <- joint_grid %>%
  left_join(
    unique_ratio_grid,
    by = "Prenatal_Ratio_Multiplier"
  )

c0_grid <- tibble::tibble(
  C0_Multiplier = C0_MULTIPLIERS,
  C0_ID = paste0(
    "C0_",
    format(C0_MULTIPLIERS, trim = TRUE)
  )
)

# ============================================================
# A3. PRECOMPUTE PHYSIOLOGY CACHES
# ============================================================

cat("\n============================================================\n")
cat("RCM/RBM SOURCE-ATTRIBUTION SENSITIVITY v1\n")
cat("============================================================\n")
cat("Analytic cohort N = ", nrow(cleaned_data), "\n", sep = "")
cat("Compounds = ", paste(pfas_cols, collapse = ", "), "\n", sep = "")
cat("Attribution reporting grid = daily\n")
cat("Formal source framework = TOTAL / PRENATAL_ONLY / POSTNATAL_ONLY\n")
cat("Joint grid = ", nrow(joint_grid), " R_CM x R_BM cells; ",
    nrow(unique_ratio_grid), " unique R_CM/R_BM ratios\n", sep = "")

phys_cache_list <- vector("list", nrow(cleaned_data))

for (i in seq_len(nrow(cleaned_data))) {
  if (
    i == 1L ||
    i %% max(1L, floor(nrow(cleaned_data) / 10L)) == 0L ||
    i == nrow(cleaned_data)
  ) {
    cat(sprintf(
      "  physiology cache: %d/%d\n",
      i,
      nrow(cleaned_data)
    ))
  }

  phys_cache_list[[i]] <- tryCatch(
    make_subject_physiology_cache(
      cleaned_data[i, , drop = FALSE]
    ),
    error = function(e) {
      warning(
        "Physiology cache failed for ID=",
        cleaned_data$ID[i],
        ": ",
        e$message
      )
      NULL
    }
  )
}

# ============================================================
# A3B. SAME PBPK EQUATIONS, LESS R INTERPRETER OVERHEAD
# ============================================================
# Derive a direct-binding function from the ORIGINAL frozen RHS expression.
# The original is function(t,state,parameters) with(as.list(parameters), { ... }).
# This factory changes only how names are bound: list2env once per ODE run,
# instead of as.list(parms) on every derivative evaluation. The inner equation
# syntax, LSODA method, tolerances and daily output grid remain untouched.
make_attribution_rhs <- function(pars) {
  rhs <- pbpk_ode_unified_tmkt
  original_body <- body(rhs)
  stopifnot(
    is.call(original_body), identical(original_body[[1L]], as.name("{")),
    is.call(original_body[[2L]]),
    identical(original_body[[2L]][[1L]], as.name("with"))
  )
  equations <- original_body[[2L]][[3L]]
  local_env <- list2env(pars, parent = environment(rhs))
  compiler::cmpfun(eval(call("function", formals(rhs), equations), envir = local_env))
}

# ============================================================
# A4. CORE FORMAL ATTRIBUTION SIMULATOR
# ============================================================

run_formal_attribution_once <- function(
    sub,
    p_name,
    subject_phys,
    renal_cache,
    prenatal_ratio_multiplier = 1,
    c0_multiplier = 1,
    return_daily = FALSE,
    ode_cache = NULL
) {

  c_milk0_base <- as.numeric(sub[[p_name]])

  if (!is.finite(c_milk0_base) || c_milk0_base <= 0) {
    return(NULL)
  }

  if (
    !is.finite(prenatal_ratio_multiplier) ||
    prenatal_ratio_multiplier <= 0
  ) {
    stop("Invalid prenatal_ratio_multiplier.")
  }

  if (!is.finite(c0_multiplier) || c0_multiplier <= 0) {
    stop("Invalid c0_multiplier.")
  }

  mix_ratio <- as.numeric(sub$mix_ratio)
  if (!is.finite(mix_ratio) || mix_ratio < 0 || mix_ratio > 1) {
    mix_ratio <- 1
  }

  c_milk0 <- c_milk0_base * c0_multiplier

  f_V_pl <- subject_phys$fns$f_V_pl
  f_V_li <- subject_phys$fns$f_V_li
  f_V_kd <- subject_phys$fns$f_V_kd
  f_V_re <- subject_phys$fns$f_V_re
  f_Q_li <- subject_phys$fns$f_Q_li
  f_Q_kd <- subject_phys$fns$f_Q_kd
  f_Q_re <- subject_phys$fns$f_Q_re
  f_KMilk <- subject_phys$fns$f_KMilk

  f_C_milk <- make_milk_function(
    p_name,
    c_milk0,
    ATTR_SCENARIO
  )

  pars <- list(
    P_li = as.numeric(pfas_params[[p_name]]["P_li"]),
    P_kd = as.numeric(pfas_params[[p_name]]["P_kd"]),
    P_re = as.numeric(pfas_params[[p_name]]["P_re"]),
    mix_ratio = mix_ratio,

    f_V_pl = f_V_pl,
    f_V_li = f_V_li,
    f_V_kd = f_V_kd,
    f_V_re = f_V_re,
    f_Q_li = f_Q_li,
    f_Q_kd = f_Q_kd,
    f_Q_re = f_Q_re,
    f_KMilk = f_KMilk,

    f_Free = renal_cache$f_Free,
    f_Tm = renal_cache$f_Tm,
    f_Kt = renal_cache$f_Kt,
    f_Qfil = renal_cache$f_Qfil,
    f_kust = renal_cache$f_kust,
    f_ku = renal_cache$f_ku,
    f_V_fil = renal_cache$f_Vfil,

    f_C_milk = function(t_days) f_C_milk(t_days)
  )

  tr <- transfer_data[
    transfer_data$Compound == p_name,
    ,
    drop = FALSE
  ]

  if (nrow(tr) != 1L) {
    stop("Transfer row not found for ", p_name)
  }

  # Baseline:
  # C_fetal = C_milk * R_CM / R_BM
  #
  # Joint multiplier grid:
  # (R_CM*m_rcm)/(R_BM*m_rbm)
  # = baseline ratio * prenatal_ratio_multiplier.
  C_fetal_plasma <-
    c_milk0 *
    (tr$R_CM / tr$R_BM) *
    prenatal_ratio_multiplier

  init_full <- c(
    Al =
      C_fetal_plasma *
      pars$P_li *
      f_V_li(0),

    Ak =
      C_fetal_plasma *
      pars$P_kd *
      f_V_kd(0),

    Ar =
      C_fetal_plasma *
      pars$P_re *
      f_V_re(0),

    Ap =
      C_fetal_plasma *
      f_V_pl(0),

    Afil = 0,
    Aust = 0,
    Aurine = 0,
    Afeces = 0,
    Ainput = 0
  )

  init_zero <- init_full
  init_zero[] <- 0

  run_ode_attr <- function(init, milk_fun) {
    pars_run <- pars
    pars_run$f_C_milk <- function(t_days) milk_fun(t_days)

    deSolve::ode(
      y = init,
      times = ATTR_TIMES_H,
      func = if (ATTR_FAST_DIRECT_RHS) make_attribution_rhs(pars_run)
             else pbpk_ode_unified_tmkt,
      parms = pars_run,
      method = "lsoda",
      rtol = ODE_RTOL,
      atol = ODE_ATOL
    )
  }

  zero_milk <- function(t_days) {
    rep(0, length(t_days))
  }

  # PERFORMANCE: the three counterfactuals are still computed with the
  # original LSODA settings and daily reporting grid. Only the results of
  # mathematically identical pre-only and post-only solves are reused within
  # this subject x compound job. In particular:
  #   PRENATAL_ONLY depends on the actual initial fetal concentration;
  #   POSTNATAL_ONLY depends on the actual milk concentration C0;
  #   TOTAL depends on both and is ALWAYS solved for each setting.
  # Using the calculated concentrations as keys (rather than rounded grid
  # multipliers) avoids unintentionally merging genuinely different inputs.
  cache_key <- function(prefix, value) {
    paste0(prefix, sprintf("%.17g", as.numeric(value)))
  }

  run_cached <- function(key, init, milk_fun) {
    if (!is.null(ode_cache) &&
        exists(key, envir = ode_cache, inherits = FALSE)) {
      return(get(key, envir = ode_cache, inherits = FALSE))
    }
    result <- tryCatch(
      run_ode_attr(init, milk_fun),
      error = function(e) NULL
    )
    # Never cache failed runs, so the original failure handling is preserved.
    if (!is.null(result) && !is.null(ode_cache)) {
      assign(key, as.data.frame(result), envir = ode_cache)
    }
    result
  }

  total_out <- tryCatch(
    run_ode_attr(init_full, f_C_milk),
    error = function(e) NULL
  )
  pre_out <- run_cached(
    cache_key("PRE_", C_fetal_plasma), init_full, zero_milk
  )
  post_out <- run_cached(
    cache_key("POST_", c_milk0), init_zero, f_C_milk
  )

  if (
    is.null(total_out) ||
    is.null(pre_out) ||
    is.null(post_out)
  ) {
    return(NULL)
  }

  total_df <- as.data.frame(total_out)
  pre_df <- as.data.frame(pre_out)
  post_df <- as.data.frame(post_out)

  names(total_df)[1:10] <- c(
    "time", "Al", "Ak", "Ar", "Ap",
    "Afil", "Aust", "Aurine", "Afeces", "Ainput"
  )

  names(pre_df)[1:10] <- names(total_df)[1:10]
  names(post_df)[1:10] <- names(total_df)[1:10]

  days <- total_df$time / 24
  age_yr <- total_df$time / 8760

  Vpl <- pmax(
    f_V_pl(age_yr),
    MIN_VOLUME_L
  )

  C_total <- total_df$Ap / Vpl
  C_pre <- pre_df$Ap / Vpl
  C_post <- post_df$Ap / Vpl

  interaction <- C_total - C_pre - C_post

  C_pre_shapley <- C_pre + interaction / 2
  C_post_shapley <- C_post + interaction / 2

  F_pre <- rep(NA_real_, length(C_total))
  F_post <- rep(NA_real_, length(C_total))

  valid_total <- is.finite(C_total) & C_total > ATTR_EPS

  F_pre[valid_total] <-
    C_pre_shapley[valid_total] /
    C_total[valid_total]

  F_post[valid_total] <-
    C_post_shapley[valid_total] /
    C_total[valid_total]

  share_sum_error <- F_pre + F_post - 1

  first_cross <- first_crossing_day(
    days,
    F_post,
    threshold = 0.5
  )

  sustained_cross <- first_sustained_crossing_day(
    days,
    F_post,
    threshold = 0.5,
    sustain_days = SUSTAIN_DAYS
  )

  endpoint_total <- interp_attr(
    days,
    C_total,
    ENDPOINT_DAYS
  )

  endpoint_pre_share <- interp_attr(
    days,
    F_pre,
    ENDPOINT_DAYS
  )

  endpoint_post_share <- interp_attr(
    days,
    F_post,
    ENDPOINT_DAYS
  )

  endpoint_interaction_rel <- interp_attr(
    days,
    ifelse(
      C_total > ATTR_EPS,
      interaction / C_total,
      NA_real_
    ),
    ENDPOINT_DAYS
  )

  # Exact mass-balance error in TOTAL simulation.
  A_total_model <-
    total_df$Al +
    total_df$Ak +
    total_df$Ar +
    total_df$Ap +
    total_df$Afil +
    total_df$Aust +
    total_df$Aurine +
    total_df$Afeces

  initial_total <-
    sum(
      init_full[
        c(
          "Al", "Ak", "Ar", "Ap",
          "Afil", "Aust"
        )
      ]
    )

  mass_expected <-
    initial_total +
    total_df$Ainput

  mass_rel_error <-
    (
      A_total_model -
      mass_expected
    ) /
    pmax(
      abs(mass_expected),
      ATTR_EPS
    )

  summary_row <- tibble::tibble(
    ID = sub$ID,
    Compound = p_name,
    Scenario = ATTR_SCENARIO,

    C0_Base_ng_L = c_milk0_base,
    C0_Multiplier = c0_multiplier,
    C0_Used_ng_L = c_milk0,

    R_CM_Baseline = tr$R_CM,
    R_BM_Baseline = tr$R_BM,
    Prenatal_Ratio_Multiplier = prenatal_ratio_multiplier,

    First_Crossover_Day = first_cross,
    Sustained14_Crossover_Day = sustained_cross,

    C_Total_1mo_ng_L = endpoint_total[1],
    C_Total_3mo_ng_L = endpoint_total[2],
    C_Total_6mo_ng_L = endpoint_total[3],
    C_Total_12mo_ng_L = endpoint_total[4],

    F_Pre_1mo = endpoint_pre_share[1],
    F_Pre_3mo = endpoint_pre_share[2],
    F_Pre_6mo = endpoint_pre_share[3],
    F_Pre_12mo = endpoint_pre_share[4],

    F_Post_1mo = endpoint_post_share[1],
    F_Post_3mo = endpoint_post_share[2],
    F_Post_6mo = endpoint_post_share[3],
    F_Post_12mo = endpoint_post_share[4],

    Interaction_Rel_1mo = endpoint_interaction_rel[1],
    Interaction_Rel_3mo = endpoint_interaction_rel[2],
    Interaction_Rel_6mo = endpoint_interaction_rel[3],
    Interaction_Rel_12mo = endpoint_interaction_rel[4],

    Min_F_Pre = safe_min(F_pre),
    Max_F_Pre = safe_max(F_pre),
    Min_F_Post = safe_min(F_post),
    Max_F_Post = safe_max(F_post),

    Max_Abs_Share_Sum_Error =
      safe_max(
        abs(
          share_sum_error
        )
      ),

    Any_Negative_Shapley_Share =
      any(
        c(
          F_pre,
          F_post
        ) <
          -SHARE_QC_TOL,
        na.rm = TRUE
      ),

    Any_Share_Above_One =
      any(
        c(
          F_pre,
          F_post
        ) >
          1 + SHARE_QC_TOL,
        na.rm = TRUE
      ),

    Max_Abs_Mass_Balance_Relative_Error =
      safe_max(
        abs(
          mass_rel_error
        )
      )
  )

  if (!return_daily) {
    return(
      list(
        summary = summary_row,
        daily = NULL
      )
    )
  }

  daily_df <- tibble::tibble(
    ID = sub$ID,
    Compound = p_name,
    Scenario = ATTR_SCENARIO,
    Time_days = days,

    C_Total_ng_L = C_total,
    C_PreOnly_ng_L = C_pre,
    C_PostOnly_ng_L = C_post,
    Interaction_ng_L = interaction,

    C_Pre_Shapley_ng_L = C_pre_shapley,
    C_Post_Shapley_ng_L = C_post_shapley,

    F_Pre = F_pre,
    F_Post = F_post,

    Prenatal_Ratio_Multiplier =
      prenatal_ratio_multiplier,

    C0_Multiplier =
      c0_multiplier
  )

  list(
    summary = summary_row,
    daily = daily_df
  )
}

# ============================================================
# A5. ONE SUBJECT x COMPOUND JOB
# ============================================================

run_source_job <- function(job) {

  sub <- job$sub
  subject_phys <- job$phys
  p_name <- job$Compound

  renal_cache <- make_renal_cache_unified(
    p_name,
    subject_phys$phys_df
  )

  # A private cache per subject/compound avoids cross-subject or
  # cross-compound mixing and is discarded after this job finishes.
  ode_cache <- new.env(parent = emptyenv(), hash = TRUE)

  summary_list <- list()
  daily_baseline <- NULL
  k <- 1L

  # ----------------------------------------------------------
  # Joint R_CM / R_BM sensitivity:
  # simulate unique prenatal-ratio multipliers at C0 multiplier = 1.
  # ----------------------------------------------------------
  for (rr in unique_ratio_grid$Prenatal_Ratio_Multiplier) {

    is_baseline <-
      abs(rr - 1) < 1e-12

    ans <- run_formal_attribution_once(
      sub = sub,
      p_name = p_name,
      subject_phys = subject_phys,
      renal_cache = renal_cache,
      prenatal_ratio_multiplier = rr,
      c0_multiplier = 1,
      ode_cache = ode_cache,
      return_daily =
        SAVE_BASELINE_DAILY_TRAJECTORIES &&
        is_baseline
    )

    if (is.null(ans)) next

    summary_list[[k]] <- ans$summary %>%
      mutate(
        Analysis_Block =
          "JOINT_RCM_RBM"
      )

    k <- k + 1L

    if (
      SAVE_BASELINE_DAILY_TRAJECTORIES &&
      is_baseline &&
      !is.null(ans$daily)
    ) {
      daily_baseline <- ans$daily
    }
  }

  # ----------------------------------------------------------
  # Shared-C0 sensitivity:
  # R_CM / R_BM baseline; measured C0 scaled in BOTH pathways.
  #
  # Avoid duplicate baseline run because it is already above.
  # ----------------------------------------------------------
  for (cm in C0_MULTIPLIERS) {

    if (abs(cm - 1) < 1e-12) {
      # Reuse the baseline summary already computed.
      base_idx <- which(
        vapply(
          summary_list,
          function(z) {
            nrow(z) == 1L &&
            z$Analysis_Block[1] == "JOINT_RCM_RBM" &&
            abs(
              z$Prenatal_Ratio_Multiplier[1] -
              1
            ) < 1e-12
          },
          logical(1)
        )
      )

      if (length(base_idx)) {
        summary_list[[k]] <-
          summary_list[[base_idx[1]]] %>%
          mutate(
            Analysis_Block =
              "SHARED_C0",
            C0_Multiplier =
              1,
            C0_Used_ng_L =
              C0_Base_ng_L
          )

        k <- k + 1L
      }

      next
    }

    ans <- run_formal_attribution_once(
      sub = sub,
      p_name = p_name,
      subject_phys = subject_phys,
      renal_cache = renal_cache,
      prenatal_ratio_multiplier = 1,
      c0_multiplier = cm,
      return_daily = FALSE,
      ode_cache = ode_cache
    )

    if (is.null(ans)) next

    summary_list[[k]] <- ans$summary %>%
      mutate(
        Analysis_Block =
          "SHARED_C0"
      )

    k <- k + 1L
  }

  list(
    summary =
      bind_rows(
        summary_list
      ),

    daily_baseline =
      daily_baseline
  )
}

# ============================================================
# A6. BUILD SUBJECT JOBS (same structure as 01_run_main_model.R)
# ============================================================
# Keep each subject and its physiology cache in one future. Run the six PFAS
# locally within that worker. This cuts job scheduling, data transport, and
# returned-list overhead without changing the simulations or their outputs.
subject_jobs_attr <- vector("list", nrow(cleaned_data))
for (i in seq_len(nrow(cleaned_data))) {
  if (is.null(phys_cache_list[[i]])) next
  subject_jobs_attr[[i]] <- list(
    sub = cleaned_data[i, , drop = FALSE],
    phys = phys_cache_list[[i]]
  )
}
subject_jobs_attr <- Filter(Negate(is.null), subject_jobs_attr)

run_subject_attribution_all <- function(subject_job) {
  sub <- subject_job$sub
  phys <- subject_job$phys
  results <- vector("list", length(pfas_cols))
  k <- 0L
  for (compound in pfas_cols) {
    c0 <- as.numeric(sub[[compound]])
    if (!is.finite(c0) || c0 <= 0) next
    k <- k + 1L
    results[[k]] <- run_source_job(list(
      sub = sub, phys = phys, Compound = compound
    ))
  }
  if (!k) return(list(summary = tibble::tibble(),
                      daily_baseline = tibble::tibble()))
  results <- results[seq_len(k)]
  list(
    summary = dplyr::bind_rows(lapply(results, `[[`, "summary")),
    daily_baseline = dplyr::bind_rows(lapply(results, `[[`, "daily_baseline"))
  )
}

cat("Attribution subjects = ", length(subject_jobs_attr),
    "; compounds per subject = ", length(pfas_cols), "\n", sep = "")

# ============================================================
# A7. CHECKPOINTED SUBJECT-LEVEL PARALLELISM (FAST3)
# ============================================================
# A checkpoint is committed only after a WHOLE subject's six-PFAS output has
# returned successfully. Stops/interrupts therefore never corrupt a saved
# subject and completed subjects are not simulated again.
# Checkpoints are separated from the published tables and figures.
ATTR_CHECKPOINT_DIR <- file.path(SENSITIVITY_OUTPUT_DIR, "checkpoints")
dir.create(ATTR_CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

# Set to 2 initially for rapid smoke testing, then use Inf for all 175.
ATTR_MAX_SUBJECTS <- as.numeric(Sys.getenv("PFAS_ATTR_MAX_SUBJECTS", "Inf"))
# Small batches allow visible progress and survive an interrupted future.
ATTR_BATCH_SUBJECTS <- as.integer(Sys.getenv("PFAS_ATTR_BATCH_SUBJECTS", "6"))
ATTR_BATCH_SUBJECTS <- max(1L, ATTR_BATCH_SUBJECTS)
ATTR_RESUME <- TRUE

# Use the actual compound, grids, patient data, and numerical settings to
# detect changes; do not silently reuse incompatible checkpoints.
attr_signature <- function() {
  list(
    cohort = cleaned_data[, c("ID", "infsex", "BW_birth_kg", "mix_ratio", pfas_cols)],
    RCM_MULTIPLIERS = RCM_MULTIPLIERS,
    RBM_MULTIPLIERS = RBM_MULTIPLIERS,
    C0_MULTIPLIERS = C0_MULTIPLIERS,
    ATTR_TIMES_H = ATTR_TIMES_H,
    ATTR_SCENARIO = ATTR_SCENARIO,
    ATTR_EPS = ATTR_EPS,
    SUSTAIN_DAYS = SUSTAIN_DAYS,
    ENDPOINT_DAYS = ENDPOINT_DAYS,
    ODE_RTOL = ODE_RTOL,
    ODE_ATOL = ODE_ATOL,
    TMKT_PARAMS = TMKT_PARAMS,
    transfer_data = transfer_data,
    pfas_params = pfas_params,
    trajectory_params = trajectory_params,
    RENAL_CACHE_STEP_DAYS = RENAL_CACHE_STEP_DAYS,
    BW_INDIVIDUALIZATION_MODE = BW_INDIVIDUALIZATION_MODE,
    BW_BIRTH_EFFECT_REMAINING_AT_12MO = BW_BIRTH_EFFECT_REMAINING_AT_12MO,
    SAVE_BASELINE_DAILY_TRAJECTORIES = SAVE_BASELINE_DAILY_TRAJECTORIES
  )
}

signature_path <- file.path(ATTR_CHECKPOINT_DIR, "run_signature.rds")
current_signature <- attr_signature()
if (file.exists(signature_path)) {
  prev_signature <- readRDS(signature_path)
  if (!identical(prev_signature, current_signature)) {
    stop("FAST3 checkpoints were generated with different settings or data. ",
         "Move/remove the checkpoints directory to start a new analysis.")
  }
} else {
  saveRDS(current_signature, signature_path)
}

subject_jobs_run <- subject_jobs_attr
if (is.finite(ATTR_MAX_SUBJECTS)) {
  if (ATTR_MAX_SUBJECTS < 1) stop("PFAS_ATTR_MAX_SUBJECTS must be >=1")
  subject_jobs_run <- head(subject_jobs_run, as.integer(ATTR_MAX_SUBJECTS))
}

subject_id <- function(job) as.integer(job$sub$ID[[1]])
checkpoint_path <- function(job) file.path(
  ATTR_CHECKPOINT_DIR, sprintf("subject_%05d.rds", subject_id(job)))
read_checkpoint <- function(job) {
  path <- checkpoint_path(job)
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(readRDS(path), error=function(e) NULL)
  if (is.null(x) || !is.list(x) ||
      !all(c("summary", "daily_baseline") %in% names(x))) return(NULL)
  x
}
commit_checkpoint <- function(job, ans) {
  dest <- checkpoint_path(job)
  tmp <- tempfile(pattern="subject_tmp_", tmpdir=ATTR_CHECKPOINT_DIR, fileext=".rds")
  saveRDS(ans, tmp, compress=FALSE)
  # Rename in the same directory to avoid partially written final files.
  ok <- file.rename(tmp, dest)
  if (!ok) {
    unlink(tmp)
    stop("Failed to save checkpoint: ", dest)
  }
  invisible(dest)
}

run_start <- Sys.time()
ready <- vapply(subject_jobs_run, function(j) !is.null(read_checkpoint(j)), logical(1))
todo <- which(!ready)
cat(sprintf("FAST3: %d/%d subjects already saved, %d pending.\n",
            sum(ready), length(subject_jobs_run), length(todo)))

if (length(todo)) {
  batches <- split(todo, ceiling(seq_along(todo) / ATTR_BATCH_SUBJECTS))
  use_parallel <- isTRUE(ATTR_USE_PARALLEL) && ATTR_WORKERS > 1L
  old_plan <- future::plan()
  if (use_parallel) {
    workers <- min(ATTR_WORKERS, ATTR_BATCH_SUBJECTS, length(todo))
    future::plan(future::multisession, workers=workers)
    cat(sprintf("FAST3: %d workers; %d subject(s) per batch; %d batches.\n",
                workers, ATTR_BATCH_SUBJECTS, length(batches)))
  }

  tryCatch({
    for (bi in seq_along(batches)) {
      ix <- batches[[bi]]
      batch_t0 <- Sys.time()
      batch_jobs <- subject_jobs_run[ix]
      batch_ans <- if (use_parallel) {
        future.apply::future_lapply(
          batch_jobs, run_subject_attribution_all,
          future.seed = TRUE,
          future.scheduling = Inf,
          future.packages = c("deSolve", "dplyr", "tidyr", "tibble", "compiler")
        )
      } else {
        lapply(batch_jobs, run_subject_attribution_all)
      }
      for (k in seq_along(ix)) commit_checkpoint(batch_jobs[[k]], batch_ans[[k]])
      completed <- sum(vapply(subject_jobs_run,
                              function(j) file.exists(checkpoint_path(j)), logical(1)))
      elapsed_sec <- as.numeric(difftime(Sys.time(), run_start, units="secs"))
      eta_min <- if (completed > sum(ready)) {
        elapsed_sec / (completed - sum(ready)) * (length(subject_jobs_run) - completed) / 60
      } else NA_real_
      cat(sprintf("FAST3 batch %d/%d: saved %d/%d subjects; batch %.1f min; ETA %.1f min\n",
                  bi, length(batches), completed, length(subject_jobs_run),
                  as.numeric(difftime(Sys.time(), batch_t0, units="mins")), eta_min))
      rm(batch_ans)
      gc(verbose=FALSE)
    }
  }, finally = {
    if (use_parallel) future::plan(old_plan)
  })
}

# Reconstruct outputs in input ID order. No changes to downstream summaries.
job_results <- lapply(subject_jobs_run, read_checkpoint)
if (any(vapply(job_results, is.null, logical(1)))) {
  stop("Attribution run is incomplete. Rerun source() to resume.")
}
elapsed_min <- as.numeric(difftime(Sys.time(), run_start, units="mins"))
cat(sprintf("FAST3 checkpointed attribution complete: %.2f minutes this session.\n", elapsed_min))

# ============================================================
# A8. COMBINE UNIQUE-RATIO AND C0 RESULTS
# ============================================================

summary_unique <- bind_rows(
  lapply(
    job_results,
    `[[`,
    "summary"
  )
)

baseline_daily <- bind_rows(
  lapply(
    job_results,
    `[[`,
    "daily_baseline"
  )
)

# Attach Ratio_ID to unique-ratio results.
joint_unique_individual <- summary_unique %>%
  filter(
    Analysis_Block ==
      "JOINT_RCM_RBM"
  ) %>%
  left_join(
    unique_ratio_grid,
    by = "Prenatal_Ratio_Multiplier"
  )

# Expand unique R_CM/R_BM-ratio simulations back to the full joint grid.
joint_individual <- joint_unique_individual %>%
  select(
    -Ratio_ID
  ) %>%
  inner_join(
    joint_grid,
    by = "Prenatal_Ratio_Multiplier",
    relationship = "many-to-many"
  ) %>%
  left_join(
    transfer_baseline,
    by = "Compound",
    suffix = c("", "_Table")
  ) %>%
  mutate(
    R_CM_Used =
      R_CM_Baseline *
      RCM_Multiplier,

    R_BM_Used =
      R_BM_Baseline *
      RBM_Multiplier,

    RCM_to_RBM_Used =
      R_CM_Used /
      R_BM_Used
  )

# QC: the intentional many-to-many join must expand each simulated
# individual/compound/ratio row to precisely its number of grid points.
expected_joint_rows <- joint_unique_individual %>%
  count(Prenatal_Ratio_Multiplier, name = "n_sim") %>%
  inner_join(
    joint_grid %>% count(Prenatal_Ratio_Multiplier, name = "n_grid"),
    by = "Prenatal_Ratio_Multiplier",
    relationship = "one-to-one"
  ) %>%
  summarise(n = sum(n_sim * n_grid)) %>%
  pull(n)

stopifnot(
  nrow(joint_individual) == expected_joint_rows,
  all(!is.na(joint_individual$Grid_ID))
)

c0_individual <- summary_unique %>%
  filter(
    Analysis_Block ==
      "SHARED_C0"
  ) %>%
  left_join(
    c0_grid,
    by = "C0_Multiplier"
  )

# ============================================================
# A9. COMPOUND-LEVEL JOINT-GRID SUMMARIES
# ============================================================

summarise_attribution_group <- function(df) {

  df %>%
    summarise(
      N_ID =
        n_distinct(
          ID
        ),

      N_First_Crossover =
        sum(
          is.finite(
            First_Crossover_Day
          )
        ),

      N_Sustained14_Crossover =
        sum(
          is.finite(
            Sustained14_Crossover_Day
          )
        ),

      Median_First_Crossover_Day =
        safe_median_attr(
          First_Crossover_Day
        ),

      P025_First_Crossover_Day =
        safe_q_attr(
          First_Crossover_Day,
          0.025
        ),

      P25_First_Crossover_Day =
        safe_q_attr(
          First_Crossover_Day,
          0.25
        ),

      P75_First_Crossover_Day =
        safe_q_attr(
          First_Crossover_Day,
          0.75
        ),

      P975_First_Crossover_Day =
        safe_q_attr(
          First_Crossover_Day,
          0.975
        ),

      Median_Sustained14_Crossover_Day =
        safe_median_attr(
          Sustained14_Crossover_Day
        ),

      P025_Sustained14_Crossover_Day =
        safe_q_attr(
          Sustained14_Crossover_Day,
          0.025
        ),

      P975_Sustained14_Crossover_Day =
        safe_q_attr(
          Sustained14_Crossover_Day,
          0.975
        ),

      Median_F_Post_1mo =
        safe_median_attr(
          F_Post_1mo
        ),

      P025_F_Post_1mo =
        safe_q_attr(
          F_Post_1mo,
          0.025
        ),

      P975_F_Post_1mo =
        safe_q_attr(
          F_Post_1mo,
          0.975
        ),

      Median_F_Post_3mo =
        safe_median_attr(
          F_Post_3mo
        ),

      P025_F_Post_3mo =
        safe_q_attr(
          F_Post_3mo,
          0.025
        ),

      P975_F_Post_3mo =
        safe_q_attr(
          F_Post_3mo,
          0.975
        ),

      Median_F_Post_6mo =
        safe_median_attr(
          F_Post_6mo
        ),

      P025_F_Post_6mo =
        safe_q_attr(
          F_Post_6mo,
          0.025
        ),

      P975_F_Post_6mo =
        safe_q_attr(
          F_Post_6mo,
          0.975
        ),

      Median_F_Post_12mo =
        safe_median_attr(
          F_Post_12mo
        ),

      P025_F_Post_12mo =
        safe_q_attr(
          F_Post_12mo,
          0.025
        ),

      P975_F_Post_12mo =
        safe_q_attr(
          F_Post_12mo,
          0.975
        ),

      Median_Interaction_Rel_1mo =
        safe_median_attr(
          Interaction_Rel_1mo
        ),

      Median_Interaction_Rel_12mo =
        safe_median_attr(
          Interaction_Rel_12mo
        ),

      Fraction_First_Crossover_Within_1_5wk =
        mean(
          First_Crossover_Day >= 7 &
          First_Crossover_Day <= 35,
          na.rm = TRUE
        ),

      Fraction_Sustained_Crossover_Within_1_5wk =
        mean(
          Sustained14_Crossover_Day >= 7 &
          Sustained14_Crossover_Day <= 35,
          na.rm = TRUE
        ),

      Any_Negative_Shapley_Share =
        any(
          Any_Negative_Shapley_Share,
          na.rm = TRUE
        ),

      Any_Share_Above_One =
        any(
          Any_Share_Above_One,
          na.rm = TRUE
        ),

      Max_Abs_Share_Sum_Error =
        safe_max(
          Max_Abs_Share_Sum_Error
        ),

      Max_Abs_Mass_Balance_Relative_Error =
        safe_max(
          Max_Abs_Mass_Balance_Relative_Error
        ),

      .groups =
        "drop"
    )
}

joint_summary <- joint_individual %>%
  group_by(
    Compound,
    Grid_ID,
    RCM_Multiplier,
    RBM_Multiplier,
    Prenatal_Ratio_Multiplier,
    R_CM_Used,
    R_BM_Used
  ) %>%
  summarise_attribution_group()

c0_summary <- c0_individual %>%
  group_by(
    Compound,
    C0_ID,
    C0_Multiplier
  ) %>%
  summarise_attribution_group()

# ============================================================
# A10. BASELINE FORMAL SOURCE-ATTRIBUTION SUMMARY
# ============================================================

baseline_individual <- joint_individual %>%
  filter(
    abs(
      RCM_Multiplier -
      1
    ) <
      1e-12,
    abs(
      RBM_Multiplier -
      1
    ) <
      1e-12
  )

baseline_summary <- baseline_individual %>%
  group_by(
    Compound
  ) %>%
  summarise_attribution_group()

# Baseline daily population trajectories.
baseline_daily_summary <- baseline_daily %>%
  group_by(
    Compound,
    Time_days
  ) %>%
  summarise(
    N_ID =
      n_distinct(
        ID
      ),

    Median_F_Post =
      safe_median_attr(
        F_Post
      ),

    P025_F_Post =
      safe_q_attr(
        F_Post,
        0.025
      ),

    P975_F_Post =
      safe_q_attr(
        F_Post,
        0.975
      ),

    Median_F_Pre =
      safe_median_attr(
        F_Pre
      ),

    P025_F_Pre =
      safe_q_attr(
        F_Pre,
        0.025
      ),

    P975_F_Pre =
      safe_q_attr(
        F_Pre,
        0.975
      ),

    Median_Interaction_Rel =
      safe_median_attr(
        ifelse(
          C_Total_ng_L >
            ATTR_EPS,
          Interaction_ng_L /
            C_Total_ng_L,
          NA_real_
        )
      ),

    .groups =
      "drop"
  )

# ============================================================
# A11. ROBUSTNESS OF THE "EARLY 1-5 WEEK" COMPOUND-LEVEL FINDING
# ============================================================

joint_robustness <- joint_summary %>%
  mutate(
    Compound_Median_First_Crossover_In_1_5wk =
      is.finite(
        Median_First_Crossover_Day
      ) &
      Median_First_Crossover_Day >= 7 &
      Median_First_Crossover_Day <= 35,

    Compound_Median_Sustained_Crossover_In_1_5wk =
      is.finite(
        Median_Sustained14_Crossover_Day
      ) &
      Median_Sustained14_Crossover_Day >= 7 &
      Median_Sustained14_Crossover_Day <= 35
  ) %>%
  group_by(
    Compound
  ) %>%
  summarise(
    N_Grid_Cells =
      n(),

    Fraction_Grid_Cells_Median_First_In_1_5wk =
      mean(
        Compound_Median_First_Crossover_In_1_5wk,
        na.rm = TRUE
      ),

    Fraction_Grid_Cells_Median_Sustained_In_1_5wk =
      mean(
        Compound_Median_Sustained_Crossover_In_1_5wk,
        na.rm = TRUE
      ),

    Min_Median_First_Crossover_Day =
      safe_min(
        Median_First_Crossover_Day
      ),

    Max_Median_First_Crossover_Day =
      safe_max(
        Median_First_Crossover_Day
      ),

    Min_Median_Sustained_Crossover_Day =
      safe_min(
        Median_Sustained14_Crossover_Day
      ),

    Max_Median_Sustained_Crossover_Day =
      safe_max(
        Median_Sustained14_Crossover_Day
      ),

    Min_Median_F_Post_1mo =
      safe_min(
        Median_F_Post_1mo
      ),

    Max_Median_F_Post_1mo =
      safe_max(
        Median_F_Post_1mo
      ),

    Min_Median_F_Post_3mo =
      safe_min(
        Median_F_Post_3mo
      ),

    Max_Median_F_Post_3mo =
      safe_max(
        Median_F_Post_3mo
      ),

    Min_Median_F_Post_6mo =
      safe_min(
        Median_F_Post_6mo
      ),

    Max_Median_F_Post_6mo =
      safe_max(
        Median_F_Post_6mo
      ),

    Min_Median_F_Post_12mo =
      safe_min(
        Median_F_Post_12mo
      ),

    Max_Median_F_Post_12mo =
      safe_max(
        Median_F_Post_12mo
      ),

    .groups =
      "drop"
  )

# All-six compound-level medians within 1-5 weeks for each grid cell.
all_compound_grid_robustness <- joint_summary %>%
  mutate(
    First_In_1_5wk =
      is.finite(
        Median_First_Crossover_Day
      ) &
      Median_First_Crossover_Day >= 7 &
      Median_First_Crossover_Day <= 35,

    Sustained_In_1_5wk =
      is.finite(
        Median_Sustained14_Crossover_Day
      ) &
      Median_Sustained14_Crossover_Day >= 7 &
      Median_Sustained14_Crossover_Day <= 35
  ) %>%
  group_by(
    Grid_ID,
    RCM_Multiplier,
    RBM_Multiplier,
    Prenatal_Ratio_Multiplier
  ) %>%
  summarise(
    N_Compounds =
      n(),

    N_First_In_1_5wk =
      sum(
        First_In_1_5wk,
        na.rm = TRUE
      ),

    All6_First_In_1_5wk =
      N_First_In_1_5wk ==
      6,

    N_Sustained_In_1_5wk =
      sum(
        Sustained_In_1_5wk,
        na.rm = TRUE
      ),

    All6_Sustained_In_1_5wk =
      N_Sustained_In_1_5wk ==
      6,

    .groups =
      "drop"
  )

# ============================================================
# A12. C0 STRUCTURAL-DIAGNOSTIC SUMMARY
# ============================================================

c0_robustness <- c0_summary %>%
  group_by(
    Compound
  ) %>%
  summarise(
    Min_Median_First_Crossover_Day =
      safe_min(
        Median_First_Crossover_Day
      ),

    Max_Median_First_Crossover_Day =
      safe_max(
        Median_First_Crossover_Day
      ),

    Range_Median_First_Crossover_Day =
      Max_Median_First_Crossover_Day -
      Min_Median_First_Crossover_Day,

    Min_Median_F_Post_1mo =
      safe_min(
        Median_F_Post_1mo
      ),

    Max_Median_F_Post_1mo =
      safe_max(
        Median_F_Post_1mo
      ),

    Range_Median_F_Post_1mo =
      Max_Median_F_Post_1mo -
      Min_Median_F_Post_1mo,

    Min_Median_F_Post_12mo =
      safe_min(
        Median_F_Post_12mo
      ),

    Max_Median_F_Post_12mo =
      safe_max(
        Median_F_Post_12mo
      ),

    Range_Median_F_Post_12mo =
      Max_Median_F_Post_12mo -
      Min_Median_F_Post_12mo,

    .groups =
      "drop"
  )

# ============================================================
# A13. ATTRIBUTION / NUMERICAL QC
# ============================================================

attribution_qc <- bind_rows(
  joint_individual %>%
    transmute(
      Analysis_Block =
        "JOINT_RCM_RBM",
      Compound,
      ID,
      Any_Negative_Shapley_Share,
      Any_Share_Above_One,
      Max_Abs_Share_Sum_Error,
      Max_Abs_Mass_Balance_Relative_Error
    ),

  c0_individual %>%
    transmute(
      Analysis_Block =
        "SHARED_C0",
      Compound,
      ID,
      Any_Negative_Shapley_Share,
      Any_Share_Above_One,
      Max_Abs_Share_Sum_Error,
      Max_Abs_Mass_Balance_Relative_Error
    )
) %>%
  group_by(
    Analysis_Block,
    Compound
  ) %>%
  summarise(
    N_Rows =
      n(),

    Any_Negative_Shapley_Share =
      any(
        Any_Negative_Shapley_Share,
        na.rm = TRUE
      ),

    Any_Share_Above_One =
      any(
        Any_Share_Above_One,
        na.rm = TRUE
      ),

    Max_Abs_Share_Sum_Error =
      safe_max(
        Max_Abs_Share_Sum_Error
      ),

    Max_Abs_Mass_Balance_Relative_Error =
      safe_max(
        Max_Abs_Mass_Balance_Relative_Error
      ),

    Mass_Balance_Pass_1e5 =
      Max_Abs_Mass_Balance_Relative_Error <=
      MASS_BALANCE_QC_TOL,

    .groups =
      "drop"
  )


# ============================================================
# A13B. ENVIRONMENTAL POLLUTION / ELSEVIER FIGURE SETTINGS
# ============================================================

# Publication defaults for all figures in this analysis:
#   * English text only
#   * Arial font
#   * Enlarged text for readability after journal scaling
#   * Clean white background
#   * PDF vector output + 600-dpi TIFF + 600-dpi PNG preview
#   * Short axis/legend labels; detailed interpretation belongs in captions

EP_FONT_FAMILY <- "Arial"
EP_BASE_SIZE <- 15

EP_THEME <- function() {
  theme_bw(
    base_size = EP_BASE_SIZE,
    base_family = EP_FONT_FAMILY
  ) +
    theme(
      plot.title = element_text(
        size = 17,
        face = "bold",
        hjust = 0,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_blank(),
      axis.title = element_text(
        size = 16
      ),
      axis.text = element_text(
        size = 13
      ),
      strip.text = element_text(
        size = 14,
        face = "bold"
      ),
      legend.title = element_text(
        size = 14
      ),
      legend.text = element_text(
        size = 13
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.25
      ),
      panel.border = element_rect(
        linewidth = 0.6
      ),
      plot.margin = margin(
        10, 12, 10, 10
      )
    )
}

save_ep_figure <- function(
    plot_object,
    filename_stem,
    width_in = 7.5,
    height_in = 6.5
) {

  png_file <- file.path(
    SENSITIVITY_OUTPUT_DIR,
    "plots",
    paste0(filename_stem, ".png")
  )

  pdf_file <- file.path(
    SENSITIVITY_OUTPUT_DIR,
    "plots",
    paste0(filename_stem, ".pdf")
  )

  tiff_file <- file.path(
    SENSITIVITY_OUTPUT_DIR,
    "plots",
    paste0(filename_stem, ".tiff")
  )

  ggsave(
    png_file,
    plot_object,
    width = width_in,
    height = height_in,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    pdf_file,
    plot_object,
    width = width_in,
    height = height_in,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  ggsave(
    tiff_file,
    plot_object,
    width = width_in,
    height = height_in,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )

  invisible(
    c(
      PNG = png_file,
      PDF = pdf_file,
      TIFF = tiff_file
    )
  )
}


# ============================================================
# A14. FIGURES — EP SUBMISSION FORMAT
# ============================================================

# Figure 1. Joint R_CM/R_BM sensitivity of crossover age
p_cross_heat <- ggplot(
  joint_summary,
  aes(
    x = factor(RBM_Multiplier),
    y = factor(RCM_Multiplier),
    fill = Median_First_Crossover_Day
  )
) +
  geom_tile(
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = ifelse(
        is.finite(Median_First_Crossover_Day),
        sprintf("%.1f", Median_First_Crossover_Day),
        "NA"
      )
    ),
    size = 4.0,
    family = EP_FONT_FAMILY
  ) +
  facet_wrap(
    ~ Compound,
    ncol = 3
  ) +
  labs(
    title = "Transfer-ratio sensitivity of source crossover age",
    x = expression(R[BM]~"multiplier"),
    y = expression(R[CM]~"multiplier"),
    fill = "Crossover age\n(days)"
  ) +
  EP_THEME() +
  theme(
    legend.position = "right"
  )

save_ep_figure(
  p_cross_heat,
  "Figure_1_RCM_RBM_crossover_heatmap_EP",
  width_in = 7.5,
  height_in = 6.8
)

# Figure 2. Joint R_CM/R_BM sensitivity of the 12-month postnatal share
p_share_heat <- ggplot(
  joint_summary,
  aes(
    x = factor(RBM_Multiplier),
    y = factor(RCM_Multiplier),
    fill = Median_F_Post_12mo
  )
) +
  geom_tile(
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = ifelse(
        is.finite(Median_F_Post_12mo),
        sprintf("%.3f", Median_F_Post_12mo),
        "NA"
      )
    ),
    size = 4.0,
    family = EP_FONT_FAMILY
  ) +
  facet_wrap(
    ~ Compound,
    ncol = 3
  ) +
  labs(
    title = "Transfer-ratio sensitivity of postnatal contribution",
    x = expression(R[BM]~"multiplier"),
    y = expression(R[CM]~"multiplier"),
    fill = "Postnatal fraction\nat 12 months"
  ) +
  EP_THEME() +
  theme(
    legend.position = "right"
  )

save_ep_figure(
  p_share_heat,
  "Figure_2_RCM_RBM_postnatal_share_12mo_heatmap_EP",
  width_in = 7.5,
  height_in = 6.8
)

# Figure 3. Shared-colostrum scaling diagnostic
p_c0 <- ggplot(
  c0_summary,
  aes(
    x = C0_Multiplier,
    y = Median_First_Crossover_Day
  )
) +
  geom_line(
    linewidth = 1.0
  ) +
  geom_point(
    size = 3.0
  ) +
  facet_wrap(
    ~ Compound,
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    title = "Shared colostrum-scaling diagnostic",
    x = expression(C[0]~"multiplier"),
    y = "Median crossover age (days)"
  ) +
  EP_THEME()

save_ep_figure(
  p_c0,
  "Figure_3_shared_C0_crossover_sensitivity_EP",
  width_in = 7.5,
  height_in = 6.5
)

# Figure 4. Baseline formal Shapley source-attribution trajectories
if (
  SAVE_BASELINE_DAILY_TRAJECTORIES &&
  nrow(baseline_daily_summary) > 0
) {

  p_source <- ggplot(
    baseline_daily_summary,
    aes(
      x = Time_days / 30.4375,
      y = Median_F_Post
    )
  ) +
    geom_ribbon(
      aes(
        ymin = P025_F_Post,
        ymax = P975_F_Post
      ),
      alpha = 0.18
    ) +
    geom_line(
      linewidth = 1.05
    ) +
    geom_hline(
      yintercept = 0.5,
      linetype = 2,
      linewidth = 0.75
    ) +
    facet_wrap(
      ~ Compound,
      ncol = 3
    ) +
    coord_cartesian(
      ylim = c(0, 1)
    ) +
    scale_x_continuous(
      breaks = seq(0, 12, by = 2),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = "Prenatal-to-postnatal source transition during infancy",
      x = "Postnatal age (months)",
      y = "Postnatal source contribution"
    ) +
    EP_THEME()

  save_ep_figure(
    p_source,
    "Figure_4_baseline_formal_source_attribution_EP",
    width_in = 7.5,
    height_in = 6.8
  )
}

# Figure 5. Range of compound-level median crossover ages across the joint grid
robustness_plot_df <- joint_robustness %>%
  mutate(
    Compound = factor(
      Compound,
      levels = pfas_cols
    ),
    Baseline_Crossover_Day =
      baseline_summary$Median_First_Crossover_Day[
        match(
          Compound,
          baseline_summary$Compound
        )
      ]
  )

p_robust <- ggplot(
  robustness_plot_df,
  aes(
    x = Compound,
    y = Baseline_Crossover_Day,
    ymin = Min_Median_First_Crossover_Day,
    ymax = Max_Median_First_Crossover_Day
  )
) +
  geom_errorbar(
    width = 0.18,
    linewidth = 0.95
  ) +
  geom_point(
    size = 3.2
  ) +
  geom_hline(
    yintercept = c(7, 35),
    linetype = 2,
    linewidth = 0.65
  ) +
  labs(
    title = "Robustness of source crossover age",
    x = NULL,
    y = "Median crossover age (days)"
  ) +
  EP_THEME() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

save_ep_figure(
  p_robust,
  "Figure_5_crossover_robustness_range_EP",
  width_in = 7.2,
  height_in = 5.6
)

# ============================================================
# A14B. FIGURE-CAPTION TEXT FOR MANUSCRIPT / SUPPLEMENT
# ============================================================

figure_captions <- tibble::tribble(
  ~Figure,
  ~Caption,

  "Figure 1",
  paste0(
    "Joint sensitivity of the prenatal-to-postnatal source crossover age to R_CM and R_BM. ",
    "Each cell shows the compound-level median first age at which the formal two-source Shapley postnatal contribution exceeded 0.5. ",
    "R_CM and R_BM were varied as multiplicative sensitivity factors around their baseline values. ",
    "The factors represent prespecified sensitivity bounds rather than statistical confidence limits."
  ),

  "Figure 2",
  paste0(
    "Joint sensitivity of the 12-month postnatal source contribution to R_CM and R_BM. ",
    "Values are compound-level medians of the formal two-source Shapley postnatal fraction. ",
    "Contributions refer only to the two modeled sources: reconstructed prenatal burden and postnatal breast-milk input."
  ),

  "Figure 3",
  paste0(
    "Shared-colostrum scaling diagnostic for source crossover age. ",
    "The measured colostrum concentration C0 was multiplied simultaneously in the reconstructed prenatal burden and the postnatal milk-input pathway. ",
    "Points show compound-level median crossover ages."
  ),

  "Figure 4",
  paste0(
    "Formal postnatal source contribution during the first year of life under baseline transfer parameters. ",
    "Solid lines show medians across simulated infants and shaded bands show the 2.5th-97.5th percentiles across individuals. ",
    "The horizontal dashed line denotes equal prenatal and postnatal contribution (0.5). ",
    "The shaded bands describe between-infant variation and are not confidence intervals."
  ),

  "Figure 5",
  paste0(
    "Robustness of compound-level median source crossover age across the joint R_CM/R_BM sensitivity grid. ",
    "Points indicate the baseline median crossover age and vertical ranges indicate the minimum and maximum compound-level median crossover ages across the tested grid. ",
    "Horizontal dashed lines mark 1 and 5 postnatal weeks."
  )
)

write.csv(
  figure_captions,
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "plots",
    "EP_figure_captions.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# A15. OUTPUT WORKBOOK
# ============================================================

notes <- tibble::tribble(
  ~Item,
  ~Definition,

  "Purpose",
  "Jointly perturb R_CM and R_BM and quantify whether the manuscript's early prenatal-to-postnatal source transition remains stable.",

  "Formal attribution",
  "TOTAL, PRENATAL_ONLY and POSTNATAL_ONLY simulations are combined using the two-source Shapley allocation I=Ctotal-Cpre-Cpost; Cpre*=Cpre+I/2; Cpost*=Cpost+I/2.",

  "Source scope",
  "Fractions refer only to reconstructed prenatal burden and postnatal breast-milk input. They do not represent total real-world postnatal exposure.",

  "R_CM/R_BM structure",
  "In the Shanghai primary cohort R_CM and R_BM enter prenatal initialization through the ratio R_CM/R_BM. Grid cells with the same multiplier ratio are algebraically equivalent for birth-burden reconstruction.",

  "Joint factor grid",
  paste0(
    "R_CM and R_BM multipliers: ",
    paste(
      RCM_MULTIPLIERS,
      collapse = ", "
    ),
    ". These are prespecified sensitivity factors, not confidence limits."
  ),

  "Shared C0 diagnostic",
  paste0(
    "Measured colostrum C0 is scaled simultaneously in prenatal reconstruction and postnatal milk input by factors ",
    paste(
      C0_MULTIPLIERS,
      collapse = ", "
    ),
    "."
  ),

  "Crossover",
  "First crossover is the first age at which formal Shapley Fpost exceeds 0.5, with linear interpolation between daily output points.",

  "Sustained crossover",
  paste0(
    "A robustness endpoint also requires Fpost > 0.5 for at least ",
    SUSTAIN_DAYS,
    " consecutive days."
  ),

  "Intervals",
  "2.5th-97.5th percentile ranges across infants represent between-infant distributions, not parameter confidence intervals.",

  "Daily grid",
  "Attribution trajectories are reported daily for crossover extraction; LSODA remains adaptive internally.",

  "Interpretation",
  "If crossover age changes strongly across plausible transfer-factor perturbations, the 1-5 week conclusion should be presented as conditional rather than as a sharply identified biological threshold."
)

out_xlsx <- file.path(
  SENSITIVITY_OUTPUT_DIR,
  "tables",
  "PFAS_RCM_RBM_source_attribution_sensitivity_v1_EP.xlsx"
)

safe_write_xlsx(
  list(
    Notes =
      notes,

    Baseline_Transfer =
      transfer_baseline,

    Joint_Grid =
      joint_grid,

    Unique_Ratio_Grid =
      unique_ratio_grid,

    Baseline_Formal_Summary =
      baseline_summary,

    Joint_RCM_RBM_Summary =
      joint_summary,

    Joint_Robustness =
      joint_robustness,

    All6_Grid_Robustness =
      all_compound_grid_robustness,

    Shared_C0_Summary =
      c0_summary,

    Shared_C0_Robustness =
      c0_robustness,

    Attribution_QC =
      attribution_qc,

    Joint_RCM_RBM_Individual =
      joint_individual,

    Shared_C0_Individual =
      c0_individual,

    Baseline_Daily_Summary =
      baseline_daily_summary,

    Figure_Captions =
      figure_captions
  ),
  out_xlsx
)

# Save raw compact results and optional individual daily baseline trajectories.
saveRDS(
  list(
    joint_individual =
      joint_individual,
    joint_summary =
      joint_summary,
    c0_individual =
      c0_individual,
    c0_summary =
      c0_summary,
    baseline_individual =
      baseline_individual,
    baseline_daily =
      baseline_daily,
    baseline_daily_summary =
      baseline_daily_summary,
    attribution_qc =
      attribution_qc
  ),
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "tables",
    "PFAS_RCM_RBM_source_attribution_sensitivity_v1_EP.rds"
  )
)

# ============================================================
# A16. CONSOLE SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("SOURCE-ATTRIBUTION SENSITIVITY COMPLETE\n")
cat("============================================================\n")

cat("\nBaseline formal Shapley source-attribution summary:\n")
print(
  baseline_summary %>%
    select(
      Compound,
      N_ID,
      Median_First_Crossover_Day,
      P025_First_Crossover_Day,
      P975_First_Crossover_Day,
      Median_Sustained14_Crossover_Day,
      Median_F_Post_1mo,
      Median_F_Post_3mo,
      Median_F_Post_6mo,
      Median_F_Post_12mo,
      Median_Interaction_Rel_1mo,
      Median_Interaction_Rel_12mo
    ),
  n = Inf,
  width = Inf
)

cat("\nJoint R_CM/R_BM robustness by compound:\n")
print(
  joint_robustness,
  n = Inf,
  width = Inf
)

cat("\nShared C0 scaling robustness:\n")
print(
  c0_robustness,
  n = Inf,
  width = Inf
)

cat("\nAll-six-compound joint-grid robustness:\n")
print(
  all_compound_grid_robustness,
  n = Inf,
  width = Inf
)

cat("\nAttribution / numerical QC:\n")
print(
  attribution_qc,
  n = Inf,
  width = Inf
)

cat(
  "\nWorkbook:\n",
  out_xlsx,
  "\n",
  sep = ""
)

cat(
  "\nPlots:\n",
  file.path(
    SENSITIVITY_OUTPUT_DIR,
    "plots"
  ),
  "\n",
  sep = ""
)

cat(
  "\nInterpretation reminder:\n"
)
cat(
  "- R_CM and R_BM are not separately identifiable from source attribution in this primary-cohort initialization; their ratio controls reconstructed prenatal burden.\n"
)
cat(
  "- Shared C0 scaling tests the structural coupling created by using the same colostrum measurement in both prenatal reconstruction and postnatal input.\n"
)
cat(
  "- The 1-5 week conclusion should be retained as a robust summary only if it remains stable across the joint transfer-factor grid.\n"
)

cat("\nDone.\n")

# ============================================================
# END OF PFAS_RCM_RBM_source_attribution_sensitivity_v1_EP.R
# ============================================================
