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
# GitHub/reproducible path configuration.
# Run this script from the repository root unless PFAS_PBPK_ROOT is set.
PROJECT_ROOT <- Sys.getenv("PFAS_PBPK_ROOT", unset = getwd())
PROJECT_ROOT <- normalizePath(
  PROJECT_ROOT,
  winslash = "/",
  mustWork = FALSE
)

file_input <- Sys.getenv(
  "PFAS_PBPK_INPUT",
  unset = file.path(
    PROJECT_ROOT,
    "data",
    "input-PBPK.xlsx"
  )
)

# Retained only for provenance/documentation. The primary infant simulation
# does not read the Abraham workbook directly.
file_abraham_reference <- file.path(
  PROJECT_ROOT,
  "data",
  "external",
  "1-s2.0-S0160412024006330-mmc2.xlsx"
)

BASE_RESULTS_DIR <- file.path(
  PROJECT_ROOT,
  "results"
)

MODEL_FREEZE_VERSION <- "PFAS_infant_PBPK_final_freeze_v1_1_GitHub"
MODEL_FREEZE_DATE <- "2026-10-15"

OUTPUT_DIR <- file.path(
  BASE_RESULTS_DIR,
  MODEL_FREEZE_VERSION
)

INPUT_MILK_UNIT <- "ng/L"
INPUT_BW_UNIT <- "g"

SIM_YEARS <- 1
DAYS_PER_YEAR <- 365.25
HOURS_PER_YEAR <- DAYS_PER_YEAR * 24
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

input_sheets <- readxl::excel_sheets(file_input)

if (!"Sheet1" %in% input_sheets) {
  stop(
    "Input workbook must contain a sheet named 'Sheet1'. Available sheets: ",
    paste(input_sheets, collapse = ", ")
  )
}

raw_data <- read_excel(
  file_input,
  sheet = "Sheet1"
)
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
SIM_END_HOURS <- SIM_YEARS * HOURS_PER_YEAR

# Weekly reporting grid plus the exact 365.25-day endpoint.
# LSODA remains adaptive internally; this only controls requested output times.
times <- sort(
  unique(
    c(
      seq(
        0,
        SIM_END_HOURS,
        by = TIME_STEP_HOURS
      ),
      SIM_END_HOURS
    )
  )
)

# Use the same 365.25-day year definition as SIM_END_HOURS.
# This keeps the exact 365.25-day endpoint at Time_yr = 1.0.
age_seq_yr <- times / HOURS_PER_YEAR

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
    age_yr <- t / HOURS_PER_YEAR
    
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
  
  t_yr <- out_df$time / HOURS_PER_YEAR
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
    SIM_YEARS * HOURS_PER_YEAR,
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
      max(f_V_pl(pre_df$time[idx12] / HOURS_PER_YEAR), MIN_VOLUME_L)
    
    post_C12 <- post_df$Ap[idx12] /
      max(f_V_pl(post_df$time[idx12] / HOURS_PER_YEAR), MIN_VOLUME_L)
    
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
# 18. PRECOMPUTE SUBJECT PHYSIOLOGY CACHES
# ============================================================
cat("\n============================================================\n")
cat("PRECOMPUTING SUBJECT PHYSIOLOGY CACHES\n")
cat("============================================================\n")

cache_start <- Sys.time()
subject_jobs <- vector("list", nrow(cleaned_data))

for (i in seq_len(nrow(cleaned_data))) {
  if (
    i == 1 ||
    i %% max(1, floor(nrow(cleaned_data) / 10)) == 0 ||
    i == nrow(cleaned_data)
  ) {
    cat(sprintf(
      "  physiology cache: %d/%d\n",
      i,
      nrow(cleaned_data)
    ))
  }
  
  phys_i <- tryCatch(
    make_subject_physiology_cache(cleaned_data[i, , drop = FALSE]),
    error = function(e) {
      warning(
        "Physiology cache failed for ID ",
        cleaned_data$ID[i],
        ": ",
        e$message
      )
      NULL
    }
  )
  
  if (!is.null(phys_i)) {
    subject_jobs[[i]] <- list(
      sub = cleaned_data[i, , drop = FALSE],
      phys = phys_i
    )
  }
}

subject_jobs <- subject_jobs[
  !vapply(subject_jobs, is.null, logical(1))
]

cache_elapsed <- as.numeric(
  difftime(Sys.time(), cache_start, units = "secs")
)

cat(sprintf(
  "Physiology cache complete in %.2f seconds for %d subjects.\n",
  cache_elapsed,
  length(subject_jobs)
))

# ============================================================
# 19. RUN ALL SIMULATIONS
# ============================================================
sim_start <- Sys.time()

if (
  USE_PARALLEL &&
  length(subject_jobs) > 1L &&
  PARALLEL_WORKERS > 1L
) {
  workers <- min(PARALLEL_WORKERS, length(subject_jobs))
  
  cat(
    "Running parallel subject simulations with ",
    workers,
    " workers...\n",
    sep = ""
  )
  
  old_plan <- future::plan()
  future::plan(future::multisession, workers = workers)
  
  results_subject <- tryCatch(
    future.apply::future_lapply(
      subject_jobs,
      run_subject_all_unified,
      future.seed = TRUE,
      future.packages = c("deSolve", "dplyr", "tibble")
    ),
    finally = future::plan(old_plan)
  )
} else {
  cat("Running serial subject simulations...\n")
  results_subject <- lapply(
    subject_jobs,
    run_subject_all_unified
  )
}

results_all <- dplyr::bind_rows(results_subject)

sim_elapsed <- as.numeric(
  difftime(Sys.time(), sim_start, units = "secs")
)

cat(sprintf(
  "Simulation completed in %.2f seconds (%.2f minutes).\n",
  sim_elapsed,
  sim_elapsed / 60
))

if (nrow(results_all) == 0) {
  stop("No PBPK simulation results were generated.")
}

write_xlsx(
  list(PBPK_Results = results_all),
  path = file.path(
    OUTPUT_DIR,
    "tables",
    "PBPK_all_S0_S1_unified_TmKt_UrinaryStorage_results.xlsx"
  )
)

# ============================================================
# 20. POST-RUN RENAL QC
# ============================================================
renal_qc_summary <- results_all %>%
  group_by(Scenario, Compound, Renal_Model) %>%
  summarise(
    Median_Reabsorption_Rate_Ratio =
      safe_median(Reabsorption_Rate_Ratio),
    Max_Reabsorption_Rate_Ratio =
      safe_max(Reabsorption_Rate_Ratio),
    Median_Urinary_Clearance_Rate_Ratio =
      safe_median(Urinary_Clearance_Rate_Ratio),
    Max_Urinary_Clearance_Rate_Ratio =
      safe_max(Urinary_Clearance_Rate_Ratio),
    Max_Abs_Mass_Balance_Error =
      safe_max(abs(Mass_Balance_Error)),
    Max_Abs_Mass_Balance_Relative_Error =
      safe_max(abs(Mass_Balance_Relative_Error)),
    BW_Min_kg = safe_min(BW_current_kg),
    BW_Max_kg = safe_max(BW_current_kg),
    Ku_Missing_Fraction = safe_mean(as.numeric(Ku_Missing)),
    N = n(),
    .groups = "drop"
  )

renal_qc_flags <- renal_qc_summary %>%
  mutate(
    Reabsorption_Ratio_Flag =
      is.finite(Max_Reabsorption_Rate_Ratio) &
      Max_Reabsorption_Rate_Ratio > REABSORPTION_RATIO_WARNING,
    Mass_Balance_Flag =
      is.finite(Max_Abs_Mass_Balance_Relative_Error) &
      Max_Abs_Mass_Balance_Relative_Error > MASS_BALANCE_REL_TOL,
    Missing_Ku_Flag =
      is.finite(Ku_Missing_Fraction) & Ku_Missing_Fraction > 0
  )

write_xlsx(
  list(
    Renal_QC_Summary = renal_qc_summary,
    Renal_QC_Flags = renal_qc_flags
  ),
  path = file.path(
    OUTPUT_DIR,
    "qc",
    "renal_QC_summary_v3.xlsx"
  )
)

cat("\n========== MASS BALANCE DIAGNOSTIC ==========\n")
mb_diag <- results_all %>%
  group_by(Scenario, Compound) %>%
  summarise(
    Max_Rel_Error = safe_max(abs(Mass_Balance_Relative_Error)),
    Median_Rel_Error = safe_median(abs(Mass_Balance_Relative_Error)),
    Max_Abs_Error_ng = safe_max(abs(Mass_Balance_Error)),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_Rel_Error))
print(mb_diag, n = Inf)

write_xlsx(
  list(Mass_Balance_Diagnostic = mb_diag),
  path = file.path(OUTPUT_DIR, "qc", "mass_balance_diagnostic_v3.xlsx")
)

if (any(renal_qc_flags$Mass_Balance_Flag, na.rm = TRUE)) {
  warning(
    "Mass-balance relative error exceeded tolerance in at least one group. ",
    "Inspect qc/renal_QC_summary_v3.xlsx and qc/mass_balance_diagnostic_v3.xlsx."
  )
}

# ============================================================
# 20B. EXTREME BIRTH-WEIGHT INFLUENCE / SENSITIVITY ANALYSIS
# ============================================================
# Primary PBPK results retain every valid subject. This secondary analysis asks
# whether subjects with birth weight >=4.5 kg materially drive population-level
# 12-month predictions. No parameter is refitted.

endpoint_all <- results_all %>%
  group_by(Scenario, ID, Compound, Renal_Model) %>%
  slice_max(order_by = Time_days, n = 1, with_ties = FALSE) %>%
  ungroup()

endpoint_non_extreme <- endpoint_all %>%
  filter(!Extreme_High_BW_4p5kg_Flag)

summarise_endpoint_population <- function(df, population_label) {
  df %>%
    group_by(Scenario, Compound, Renal_Model) %>%
    summarise(
      Population = population_label,
      N_ID = n_distinct(ID),
      Median_C_Plasma_12mo = median(C_Plasma, na.rm = TRUE),
      P025_C_Plasma_12mo = quantile(C_Plasma, 0.025, na.rm = TRUE, names = FALSE),
      P975_C_Plasma_12mo = quantile(C_Plasma, 0.975, na.rm = TRUE, names = FALSE),
      Median_A_Body_Burden_12mo = median(A_Body_Burden, na.rm = TRUE),
      .groups = "drop"
    )
}

extreme_bw_endpoint_sensitivity <- bind_rows(
  summarise_endpoint_population(endpoint_all, "Primary: all valid birth weights"),
  summarise_endpoint_population(endpoint_non_extreme, "Sensitivity: exclude BW >=4.5 kg")
)

extreme_bw_influence_comparison <- extreme_bw_endpoint_sensitivity %>%
  select(Scenario, Compound, Renal_Model, Population, N_ID, Median_C_Plasma_12mo) %>%
  pivot_wider(
    names_from = Population,
    values_from = c(N_ID, Median_C_Plasma_12mo)
  )

write_xlsx(
  list(
    Endpoint_Sensitivity = extreme_bw_endpoint_sensitivity,
    Influence_Comparison = extreme_bw_influence_comparison,
    Extreme_Subjects = cleaned_data %>%
      filter(Extreme_High_BW_4p5kg_Flag) %>%
      select(ID, infsex, BW_birth_kg, Macrosomia_4kg_Flag, Extreme_High_BW_4p5kg_Flag, Very_Extreme_BW_5kg_Flag)
  ),
  path = file.path(OUTPUT_DIR, "qc", "extreme_birth_weight_sensitivity_v3.xlsx")
)

# ============================================================
# 21. RENAL MECHANISTIC DIAGNOSTICS
# ============================================================
renal_mechanistic_summary <- results_all %>%
  group_by(Scenario, Compound, Time_days) %>%
  summarise(
    Median_GFR_L_h = median(GFR_L_h, na.rm = TRUE),
    Lower95_GFR_L_h = quantile(GFR_L_h, 0.025, na.rm = TRUE, names = FALSE),
    Upper95_GFR_L_h = quantile(GFR_L_h, 0.975, na.rm = TRUE, names = FALSE),
    
    Median_Tm_ug_h = median(Tm_ug_h, na.rm = TRUE),
    Lower95_Tm_ug_h = quantile(Tm_ug_h, 0.025, na.rm = TRUE, names = FALSE),
    Upper95_Tm_ug_h = quantile(Tm_ug_h, 0.975, na.rm = TRUE, names = FALSE),
    
    Median_Kt_ug_L = median(Kt_ug_L, na.rm = TRUE),
    Median_Free = median(Free_fraction, na.rm = TRUE),
    
    Median_R_Filtration = median(R_Filtration, na.rm = TRUE),
    Median_R_Reabsorption = median(R_Reabsorption, na.rm = TRUE),
    Median_R_Urine_Excretion = safe_median(R_Urine_Excretion),
    
    Median_Reabsorption_Rate_Ratio = safe_median(Reabsorption_Rate_Ratio),
    Median_Urinary_Clearance_Rate_Ratio = safe_median(Urinary_Clearance_Rate_Ratio),
    
    Median_A_Urine_Storage = median(A_Urine_Storage, na.rm = TRUE),
    .groups = "drop"
  )

write_xlsx(
  list(Renal_Mechanistic_Summary = renal_mechanistic_summary),
  path = file.path(
    OUTPUT_DIR,
    "qc",
    "renal_mechanistic_diagnostics_v3.xlsx"
  )
)

# ============================================================
# 22. MILK TRAJECTORY OUTPUT
# ============================================================
milk_baseline <- cleaned_data %>%
  select(ID, all_of(pfas_cols)) %>%
  pivot_longer(
    cols = all_of(pfas_cols),
    names_to = "Compound",
    values_to = "C0"
  ) %>%
  mutate(C0 = suppressWarnings(as.numeric(C0))) %>%
  filter(is.finite(C0), C0 > 0)

validation_times <- sort(
  unique(
    c(
      seq(
        0,
        SIM_YEARS * DAYS_PER_YEAR,
        by = 7
      ),
      SIM_YEARS * DAYS_PER_YEAR
    )
  )
)

milk_validation <- milk_baseline %>%
  crossing(
    Scenario = SCENARIOS,
    Time_days = validation_times
  ) %>%
  rowwise() %>%
  mutate(
    Milk_C = {
      f <- make_milk_function(
        Compound,
        C0,
        Scenario
      )
      f(Time_days)[1]
    }
  ) %>%
  ungroup() %>%
  filter(is.finite(Milk_C), Milk_C >= 0)

write_xlsx(
  list(Milk_Trajectories = milk_validation),
  path = file.path(
    OUTPUT_DIR,
    "milk_trajectories",
    "milk_validation_subject_scenario_time_v3.xlsx"
  )
)

# ============================================================
# 23. METRICS
# ============================================================
trapz <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2) return(NA_real_)
  sum(
    diff(x) *
      (y[-1] + y[-length(y)]) / 2
  )
}

calc_metrics <- function(res_df) {
  res_df %>%
    group_by(
      Scenario,
      ID,
      Compound,
      Renal_Model
    ) %>%
    arrange(Time_yr, .by_group = TRUE) %>%
    summarise(
      C0_Plasma = C_Plasma[1],
      Cmax_Plasma = max(C_Plasma, na.rm = TRUE),
      Tmax_Plasma_days = Time_days[which.max(C_Plasma)][1],
      
      AUC_Plasma_0_3mo = trapz(
        Time_yr[Time_yr <= 0.25],
        C_Plasma[Time_yr <= 0.25]
      ),
      AUC_Plasma_0_6mo = trapz(
        Time_yr[Time_yr <= 0.5],
        C_Plasma[Time_yr <= 0.5]
      ),
      AUC_Plasma_0_12mo = trapz(
        Time_yr[Time_yr <= 1],
        C_Plasma[Time_yr <= 1]
      ),
      
      C12mo_Plasma =
        C_Plasma[which.min(abs(Time_yr - 1))][1],
      
      T50_birth_burden_days = {
        burden <- A_Body_Burden
        initial_burden <- burden[1]
        idx_drop <- which(
          burden <= 0.5 * initial_burden
        )
        if (length(idx_drop)) {
          Time_days[idx_drop[1]]
        } else {
          NA_real_
        }
      },
      
      Reabsorption_Rate_Ratio_12mo =
        Reabsorption_Rate_Ratio[
          which.min(abs(Time_yr - 1))
        ][1],
      
      Urinary_Clearance_Rate_Ratio_12mo =
        Urinary_Clearance_Rate_Ratio[
          which.min(abs(Time_yr - 1))
        ][1],
      
      Max_Abs_Mass_Balance_Relative_Error = max(
        abs(Mass_Balance_Relative_Error),
        na.rm = TRUE
      ),
      
      Pre_Contribution_Approx_12mo =
        Pre_Contribution_Approx_12mo[
          which.min(abs(Time_yr - 1))
        ][1],
      
      Post_Contribution_Approx_12mo =
        Post_Contribution_Approx_12mo[
          which.min(abs(Time_yr - 1))
        ][1],
      
      Ku_Missing = first(Ku_Missing),
      .groups = "drop"
    )
}

metrics <- calc_metrics(results_all)

write_xlsx(
  list(Individual_Metrics = metrics),
  path = file.path(
    OUTPUT_DIR,
    "tables",
    "individual_PBPK_metrics_unified_TmKt_UrinaryStorage.xlsx"
  )
)

metric_summary <- metrics %>%
  group_by(Scenario, Compound, Renal_Model) %>%
  summarise(
    C0_median = median(C0_Plasma, na.rm = TRUE),
    Cmax_median = median(Cmax_Plasma, na.rm = TRUE),
    Tmax_median_days = median(Tmax_Plasma_days, na.rm = TRUE),
    AUC0_3_median = median(AUC_Plasma_0_3mo, na.rm = TRUE),
    AUC0_6_median = median(AUC_Plasma_0_6mo, na.rm = TRUE),
    AUC0_12_median = median(AUC_Plasma_0_12mo, na.rm = TRUE),
    C12_median = median(C12mo_Plasma, na.rm = TRUE),
    T50_median_days = median(T50_birth_burden_days, na.rm = TRUE),
    
    Reabsorption_Rate_Ratio_12mo_median = median(
      Reabsorption_Rate_Ratio_12mo,
      na.rm = TRUE
    ),
    
    Urinary_Clearance_Rate_Ratio_12mo_median = median(
      Urinary_Clearance_Rate_Ratio_12mo,
      na.rm = TRUE
    ),
    
    Max_Mass_Balance_Rel_Error = max(
      Max_Abs_Mass_Balance_Relative_Error,
      na.rm = TRUE
    ),
    
    Pre_contrib_approx_median = median(
      Pre_Contribution_Approx_12mo,
      na.rm = TRUE
    ),
    Post_contrib_approx_median = median(
      Post_Contribution_Approx_12mo,
      na.rm = TRUE
    ),
    
    Ku_Missing_Fraction = mean(Ku_Missing, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

write_xlsx(
  list(Population_Metric_Summary = metric_summary),
  path = file.path(
    OUTPUT_DIR,
    "tables",
    "population_metric_summary_unified_TmKt_UrinaryStorage.xlsx"
  )
)

# ============================================================
# 24. PARAMETER / MODEL DOCUMENTATION
# ============================================================
renal_parameter_documentation <- TMKT_PARAMS %>%
  transmute(
    Compound = Compound,
    Renal_Model = "Unified_Tm_Kt_UrinaryStorage",
    Tmc_ug_h_bw075 = Tmc_ug_h_bw075,
    Tmc_derivation = Tm_basis,
    Tm_70kg_ug_h = Tmc_ug_h_bw075 * 70^0.75,
    Kt_ug_L = Kt_ug_L,
    Free_fraction = Free_fraction,
    KurineC_h_bwminus025 = KurineC_h_bwminus025,
    Kurine_source = Kurine_source,
    Tm_source = Tm_source,
    Parameter_class = Parameter_class
  )

write_xlsx(
  list(Renal_Parameters = renal_parameter_documentation),
  path = file.path(
    OUTPUT_DIR,
    "tables",
    "unified_renal_parameter_documentation_v3.xlsx"
  )
)

model_documentation <- tibble::tribble(
  ~Item, ~Definition,
  
  "Renal structure",
  "All six PFAS use the same plasma -> filtrate -> saturable Tm/Kt reabsorption -> urinary-storage -> urine architecture.",
  
  "Glomerular filtration",
  "R_filtration = Qfil × Free × C_plasma; Qfil is the age-dependent infant GFR in the current implementation.",
  
  "Saturable reabsorption",
  "R_reabsorption = Tm × C_filtrate / (Kt + C_filtrate).",
  
  "Filtrate-to-urinary-storage transfer",
  "R_ust = Kust × A_filtrate; Kust is a physiological residence-rate approximation Qfil/Vfil and is not the chemical-specific Ku.",
  
  "Urinary excretion",
  "R_urine = Ku × A_urinary_storage; Ku = KurineC × BW^(-0.25) when chemical-specific KurineC is available.",
  
  "PFOA",
  "EFSA/PFAS PBK family: Tmc 6000 ug/h/kg^0.75; Kt 55 ug/L; Free 0.020; KurineC 0.0003.",
  
  "PFOS",
  "EFSA/PFAS PBK family: Tmc 3500 ug/h/kg^0.75; Kt 23 ug/L; Free 0.025; KurineC 0.001.",
  
  "PFNA",
  "Fragki/Iulini implementation: Tmc 7900 ug/h/kg^0.75; Kt 55 ug/L; Free 0.020; KurineC 0.0003.",
  
  "PFHxS",
  "Fragki/Iulini implementation: Tmc 7000 ug/h/kg^0.75; Kt 23 ug/L; Free 0.025; KurineC 0.001.",
  
  "PFDA",
  "Fàbrega-informed Kt/Free with effective Tmc = 14.700113 ug/h/kg^0.75 carried forward from Abraham-2024 adult-v5 minimal Tmc calibration; KurineC remains non-identifiable from plasma.",
  
  "PFUdA",
  "Fàbrega-informed Kt/Free with effective Tmc = 260.217197 ug/h/kg^0.75 carried forward from Abraham-2024 adult-v5 minimal Tmc calibration; KurineC remains non-identifiable from plasma.",
  
  "Infant GFR",
  "GFR (mL/min) = -17.74 + 99.054*BSA - 6.1604*BSA^2. v3 removes the erroneous extra ×0.06 factor present in v2.",
  
  "Placental transfer R_CM",
  "Median cord-to-maternal whole-blood ratios from Liu et al. 2024 Table 3: PFOA 0.290, PFNA 0.230, PFDA 0.170, PFUdA/PFUnDA 0.230, PFHxS 0.220, PFOS 0.090.",
  
  "Breast-milk transfer R_BM",
  "Chemical-specific milk/maternal-serum ratios retained from Zheng et al. 2022 parameterization: PFOA 0.053, PFNA 0.103, PFDA 0.023, PFUdA 0.022, PFHxS 0.047, PFOS 0.042. Liu et al. 2024 cites Zheng et al. as the matched maternal-serum/cord-serum/breast-milk transfer study.",
  
  "BW input",
  "Input BW is converted to kilograms before physiological scaling.",
  
  "First-order plasma clearance",
  "Not used. Elimination is represented by mechanistic renal filtration, saturable reabsorption, urinary-storage transfer, and urinary excretion.",
  
  "Prenatal/postnatal attribution",
  "Current output retains approximate counterfactual ratios for backward compatibility; final analysis should use explicit total, prenatal-only, and breastfeeding-only nonlinear simulations.",
  
  "Infant transporter maturation",
  "No direct infant-specific OAT/PTC maturation parameterization is assumed.",
  
  "PFDA/PFUdA Ku",
  "No GFR/Vfil surrogate is used as chemical-specific Ku. Adult v5 demonstrated KurineC is structurally non-identifiable from plasma in this model; urinary outputs for PFDA/PFUdA remain provisional until urine/fecal data inform Ku.",
  
  "Quality control",
  "Mass balance, renal flow, Tm/Kt, reabsorption-rate ratio, urinary-clearance ratio and missing-Ku flags are generated automatically."
)

write_xlsx(
  list(Model_Documentation = model_documentation),
  path = file.path(
    OUTPUT_DIR,
    "tables",
    "model_documentation_unified_TmKt_UrinaryStorage_v3.xlsx"
  )
)

# ============================================================
# 25. POPULATION SUMMARIES + PBPK PLOTS
# ============================================================
compartments <- c(
  "C_Plasma",
  "C_Liver",
  "C_Kidney",
  "C_Rest",
  "C_Filtrate"
)

summarise_interval <- function(df, value_col) {
  df %>%
    group_by(
      Scenario,
      Compound,
      Renal_Model,
      Time_yr
    ) %>%
    summarise(
      Median = safe_median(
        .data[[value_col]]
      ),
      Lower95 = safe_quantile(
        .data[[value_col]],
        0.025
      ),
      Upper95 = safe_quantile(
        .data[[value_col]],
        0.975
      ),
      Mean = safe_mean(
        .data[[value_col]]
      ),
      N = sum(
        is.finite(
          .data[[value_col]]
        )
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Compartment = value_col
    )
}

summary_long <- bind_rows(
  lapply(
    compartments,
    function(x) summarise_interval(results_all, x)
  )
)

write_xlsx(
  list(Population_Summary = summary_long),
  path = file.path(
    OUTPUT_DIR,
    "tables",
    "PBPK_population_summary_95_UI_unified_TmKt_UrinaryStorage.xlsx"
  )
)

make_compound_plot <- function(
    compound,
    compartment = "C_Plasma"
) {
  plot_df <- summary_long %>%
    filter(
      Compound == compound,
      Compartment == compartment
    )
  
  if (!nrow(plot_df)) return(NULL)
  
  ggplot(
    plot_df,
    aes(
      x = Time_yr,
      y = Median,
      linetype = Scenario,
      group = Scenario
    )
  ) +
    geom_ribbon(
      aes(
        ymin = Lower95,
        ymax = Upper95,
        fill = Scenario
      ),
      alpha = 0.16,
      colour = NA
    ) +
    geom_line(linewidth = 1.0) +
    labs(
      title = paste0(
        compound,
        " — ",
        gsub("C_", "", compartment)
      ),
      subtitle = "Unified renal Tm/Kt PBPK model with urinary-storage compartment",
      x = "Postnatal age (years)",
      y = "Concentration (ng/L)",
      linetype = "Scenario",
      fill = "Scenario"
    ) +
    scale_x_continuous(
      limits = c(0, SIM_YEARS),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
}

for (comp in pfas_cols) {
  for (compartment in compartments) {
    p <- make_compound_plot(comp, compartment)
    if (is.null(p)) next
    
    base_name <- paste0(
      comp,
      "_",
      gsub("C_", "", compartment),
      "_unified_TmKt_UrinaryStorage_S0_S1"
    )
    
    ggsave(
      filename = file.path(
        OUTPUT_DIR,
        "plots",
        paste0(base_name, ".png")
      ),
      plot = p,
      width = PLOT_WIDTH,
      height = PLOT_HEIGHT,
      dpi = PLOT_DPI,
      units = "in"
    )
    
    ggsave(
      filename = file.path(
        OUTPUT_DIR,
        "plots",
        paste0(base_name, ".pdf")
      ),
      plot = p,
      width = PLOT_WIDTH,
      height = PLOT_HEIGHT,
      units = "in",
      device = cairo_pdf
    )
  }
}

# ============================================================
# 26. MILK TRAJECTORY SUMMARY / PLOTS
# ============================================================
milk_plot_data <- milk_validation %>%
  group_by(
    Scenario,
    Compound,
    Time_days
  ) %>%
  summarise(
    Median_Milk_C = median(Milk_C, na.rm = TRUE),
    Lower95 = quantile(
      Milk_C,
      0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    Upper95 = quantile(
      Milk_C,
      0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  )

write_xlsx(
  list(Milk_Trajectories = milk_plot_data),
  path = file.path(
    OUTPUT_DIR,
    "tables",
    "milk_trajectory_summary_unified_TmKt_UrinaryStorage.xlsx"
  )
)

for (comp in pfas_cols) {
  df <- milk_plot_data %>%
    filter(Compound == comp)
  
  if (!nrow(df)) next
  
  p <- ggplot(
    df,
    aes(
      x = Time_days,
      y = Median_Milk_C,
      linetype = Scenario,
      group = Scenario
    )
  ) +
    geom_ribbon(
      aes(
        ymin = Lower95,
        ymax = Upper95,
        fill = Scenario
      ),
      alpha = 0.14,
      colour = NA
    ) +
    geom_line(linewidth = 1.1) +
    labs(
      title = paste0(
        comp,
        " breast-milk concentration trajectory"
      ),
      subtitle = "Unified renal Tm/Kt PBPK model input; S0 and S1",
      x = "Postnatal age (days)",
      y = "Breast-milk concentration (ng/L)"
    ) +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    file.path(
      OUTPUT_DIR,
      "milk_trajectories",
      paste0(
        comp,
        "_milk_unified_TmKt_UrinaryStorage_S0_S1.png"
      )
    ),
    p,
    width = PLOT_WIDTH,
    height = PLOT_HEIGHT,
    dpi = PLOT_DPI,
    units = "in"
  )
  
  ggsave(
    file.path(
      OUTPUT_DIR,
      "milk_trajectories",
      paste0(
        comp,
        "_milk_unified_TmKt_UrinaryStorage_S0_S1.pdf"
      )
    ),
    p,
    width = PLOT_WIDTH,
    height = PLOT_HEIGHT,
    units = "in",
    device = cairo_pdf
  )
}

# ============================================================
# 27. FINAL REPRODUCIBILITY SUMMARY
# ============================================================
fmt_chr <- function(x) {
  if (length(x) == 0 || is.null(x)) return(NA_character_)
  as.character(x)
}

fmt_num <- function(x, digits = 6) {
  z <- suppressWarnings(as.numeric(x))[1]
  if (!is.finite(z)) return(NA_character_)
  formatC(z, format = "fg", digits = digits, flag = "#")
}

repro_summary <- tibble::tribble(
  ~Item, ~Value,
  "Input file", fmt_chr(file_input),
  "Abraham reference file", fmt_chr(file_abraham_reference),
  "Simulation years", fmt_num(SIM_YEARS, 6),
  "Output time step (h)", fmt_num(TIME_STEP_HOURS, 6),
  "Renal cache step (days)", fmt_num(RENAL_CACHE_STEP_DAYS, 6),
  "TM scaling", "Tm = Tmc × BW^(0.75)",
  "Ku scaling", "Ku = KurineC × BW^(-0.25)",
  "Kust definition", "Kust = Qfil / Vfil",
  "Filtration flow", "Qfil = age-dependent infant GFR",
  "Renal storage", "Afiltrate -> urinary storage -> cumulative urine",
  "Infant GFR equation", "GFR = -17.74 + 99.054*BSA - 6.1604*BSA^2 mL/min (no ×0.06)",
  "R_CM source", "Liu et al. 2024 J Hazard Mater Table 3 median cord/maternal whole-blood ratios",
  "R_BM source", "Zheng et al. 2022 matched maternal serum / breast milk transfer; chemical-specific values retained",
  "PFDA Tmc", "14.700113 ug/h/kg^0.75; Fàbrega-informed, Abraham-2024 adult-v5 calibrated",
  "PFUdA Tmc", "260.217197 ug/h/kg^0.75; Fàbrega-informed, Abraham-2024 adult-v5 calibrated",
  "PFDA Ku", "Not identifiable from plasma; no chemical-specific value assigned",
  "PFUdA Ku", "Not identifiable from plasma; no chemical-specific value assigned",
  "Fecal route", "Absent in the current base infant model",
  "Prenatal/postnatal attribution", "Approximate backward-compatible ratio only; final analysis should use nonlinear counterfactual simulations",
  "Independent validation", "Not performed in this infant simulation script"
)

write_xlsx(
  list(Reproducibility_Summary = repro_summary),
  path = file.path(
    OUTPUT_DIR,
    "tables",
    "reproducibility_summary_v3.xlsx"
  )
)

# ============================================================
# 28. FINAL REPORT
# ============================================================
cat("\n\n============================================================\n")
cat("UNIFIED PFAS INFANT PBPK MODEL v3 COMPLETED\n")
cat("============================================================\n")
cat(
  "Scenarios: ",
  paste(SCENARIOS, collapse = ", "),
  "\n",
  sep = ""
)
cat(
  "Compounds: ",
  paste(pfas_cols, collapse = ", "),
  "\n",
  sep = ""
)
cat("Renal structure: Plasma -> Filtrate -> Tm/Kt reabsorption + urinary storage -> urine\n")
cat("Filtration flow: literature GFR = -17.74 + 99.054*BSA - 6.1604*BSA^2 mL/min\n")
cat("Kust: Qfil/Vfil residence-rate approximation\n")
cat("Ku: chemical-specific KurineC × BW^(-0.25); PFDA/PFUdA Ku remain unidentified from plasma\n")
cat("Output directory: ", OUTPUT_DIR, "\n", sep = "")
cat("Subjects simulated: ", length(subject_jobs), "\n", sep = "")
cat("Simulation time: ", round(sim_elapsed / 60, 2), " min\n", sep = "")

cat("\nParameter classes:\n")
print(
  TMKT_PARAMS %>%
    select(
      Compound,
      Tm_source,
      Kurine_source,
      Parameter_class
    ),
  row.names = FALSE
)

cat("\nKey QC files:\n")
cat(
  file.path(
    OUTPUT_DIR,
    "qc",
    "renal_QC_summary_v3.xlsx"
  ),
  "\n"
)
cat(
  file.path(
    OUTPUT_DIR,
    "qc",
    "renal_mechanistic_diagnostics_v3.xlsx"
  ),
  "\n"
)

cat("\nDone.\n")

# ============================================================
# REFERENCES / MODEL BASIS
# ============================================================
# 1. Andersen ME, Clewell HJ 3rd, Tan YM, Butenhoff JL, Olsen GW.
#    Pharmacokinetic modeling of saturable, renal resorption of
#    perfluoroalkylacids in monkeys—probing the determinants of long
#    plasma half-lives. Toxicology. 2006;227(1-2):156-164.
#    DOI: 10.1016/j.tox.2006.08.004.
#
# 2. Loccisano AE, Campbell JL Jr, Andersen ME, Clewell HJ 3rd.
#    Evaluation and prediction of pharmacokinetics of PFOA and PFOS in
#    the monkey and human using a PBPK model.
#    Regul Toxicol Pharmacol. 2011;59(1):157-175.
#    DOI: 10.1016/j.yrtph.2010.12.004.
#
# 3. Loccisano AE, Campbell JL Jr, Butenhoff JL, Andersen ME, Clewell HJ 3rd.
#    Evaluation of the transfer of PFOA and PFOS through the placenta,
#    maternal milk and developmental exposure in rats using a PBPK model.
#    Reprod Toxicol. 2012;33(4):468-490.
#
# 4. EFSA CONTAM Panel.
#    Risk to human health related to the presence of PFOS and PFOA in food.
#    EFSA Journal. 2018;16(12):5194.
#    DOI: 10.2903/j.efsa.2018.5194.
#
# 5. EFSA CONTAM Panel.
#    Risk to human health related to the presence of perfluoroalkyl
#    substances in food. EFSA Journal. 2020;18(9):6223.
#    DOI: 10.2903/j.efsa.2020.6223.
#
# 6. Bernstein K, Kapraun DF, Schlosser PM et al.
#    Model Template Approach for Rapid Evaluation and Application of
#    Physiologically Based Pharmacokinetic Models for Use in Human Health
#    Risk Assessments: A Case Study on PFAS.
#    Toxicological Sciences. 2021;182(2):215-228.
#
# 7. Kim SJ, Choi EJ, Choi GW, Lee YB, Cho HY.
#    Exploring sex differences in human health risk assessment for PFNA
#    and PFDA using a PBPK model.
#    Archives of Toxicology. 2019;93:311-330.
#    DOI: 10.1007/s00204-018-2365-y.
#
# 8. Fragki S, Louisse J, Bokkers B, Luijten M.
#    New approach methodologies: a quantitative in vitro to in vivo
#    extrapolation case study with PFASs.
#    Food and Chemical Toxicology. 2023;172:113559.
#    DOI: 10.1016/j.fct.2022.113559.
#
# 9. Iulini M et al.
#    Advancing PFAS risk assessment: Integrative approaches using
#    agent-based modelling and physiologically-based kinetic for
#    environmental and health safety.
#    Computational and Structural Biotechnology Journal. 2024;23:2763-2778.
#    DOI: 10.1016/j.csbj.2024.06.036.
#
# 10. Fàbrega F et al. Human PBPK modeling of multiple PFAS including
#     PFDA and PFUnDA/PFUdA (2015).
#     PFDA and PFUnDA Tm/Kt/Free values are retained as literature-derived
#     or derived quantities and are explicitly flagged uncertain.
#
# 11. Blomberg et al. 2023.
#     Longitudinal PFAS trajectories in breast milk, including
#     PFAS-specific early-lactation and monthly changes.
#
# 12. Liu L, Yan P, Liu X, et al. Profiles and transplacental transfer of
#     per- and polyfluoroalkyl substances in maternal and umbilical cord
#     blood: A birth cohort study in Zhoushan, Zhejiang Province, China.
#     J Hazard Mater. 2024;466:133501. DOI:10.1016/j.jhazmat.2024.133501.
#     Used for v3 R_CM median cord/maternal whole-blood ratios.
#
# 13. Zheng P, Liu Y, An Q, et al. Prenatal and postnatal exposure to
#     emerging and legacy per-/polyfluoroalkyl substances: levels and
#     transfer in maternal serum, cord serum, and breast milk.
#     Sci Total Environ. 2022;812:152446.
#     DOI:10.1016/j.scitotenv.2021.152446.
#     Retained as the prior R_BM provenance for PFOA/PFDA/PFUdA/PFOS.
#     PFNA/PFHxS final-freeze R_BM values were updated to independent
#     pooled-literature milk:serum ratios after targeted transfer audit.
#
# ============================================================
# END OF SCRIPT
# ============================================================


# ============================================================
# FINAL-FREEZE QUALIFICATION SUMMARY
# ============================================================
#
# This block summarizes only information produced by THIS frozen main run.
# It does not fabricate external-validation metrics. Study-specific Varsi
# re-evaluation should be run separately with the frozen transfer_data values.
# ============================================================

freeze_qualification_run <- renal_qc_flags %>%
  transmute(
    Scenario,
    Compound,
    Max_Abs_Mass_Balance_Relative_Error,
    Mass_Balance_Pass =
      !Mass_Balance_Flag,
    Reabsorption_Diagnostic_Flag =
      Reabsorption_Ratio_Flag,
    Missing_Ku_Flag
  )

freeze_endpoint_12m <- results_all %>%
  filter(
    abs(
      Time_days - max(Time_days, na.rm = TRUE)
    ) < 1e-8
  ) %>%
  group_by(
    Scenario,
    Compound
  ) %>%
  summarise(
    N_ID = n_distinct(ID),
    Median_Plasma_ng_L =
      safe_median(C_Plasma),
    P025_Plasma_ng_L =
      safe_quantile(C_Plasma, 0.025),
    P975_Plasma_ng_L =
      safe_quantile(C_Plasma, 0.975),
    Median_Body_Burden_ng =
      safe_median(A_Body_Burden),
    .groups = "drop"
  )

final_transfer_change_log <- tibble::tribble(
  ~Compound, ~Previous_v3_R_BM, ~Final_Freeze_R_BM, ~Relative_to_Previous, ~Reason,

  "PFNA",
  0.103,
  0.0300,
  0.0300 / 0.103,
  "Independent pooled literature plus targeted transfer audit; not fitted to Varsi infant outcome.",

  "PFHxS",
  0.047,
  0.0133,
  0.0133 / 0.047,
  "Independent pooled literature plus targeted transfer audit; not fitted to Varsi infant outcome."
)

safe_write_xlsx(
  list(
    Freeze_Manifest = model_freeze_manifest,
    Final_Transfer_Parameters = transfer_data,
    Transfer_Change_Log = final_transfer_change_log,
    Main_Run_QC = freeze_qualification_run,
    Endpoint_12m = freeze_endpoint_12m
  ),
  file.path(
    OUTPUT_DIR,
    "tables",
    "PFAS_infant_PBPK_final_freeze_v1_summary.xlsx"
  )
)

cat("\n============================================================\n")
cat("FINAL MODEL FREEZE COMPLETED\n")
cat("============================================================\n")
cat("Version: ", MODEL_FREEZE_VERSION, "\n", sep = "")
cat("Freeze date: ", MODEL_FREEZE_DATE, "\n", sep = "")
cat("\nFinal PFNA R_BM: 0.0300\n")
cat("Final PFHxS R_BM: 0.0133\n")
cat("\nMain freeze summary:\n")
cat(
  file.path(
    OUTPUT_DIR,
    "tables",
    "PFAS_infant_PBPK_final_freeze_v1_summary.xlsx"
  ),
  "\n",
  sep = ""
)
cat(
  "\nIMPORTANT: Do not tune frozen parameters to subsequent validation outcomes.\n"
)


# ============================================================
# REPRODUCIBILITY RECORD
# ============================================================
writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)

cat(
  "Session information written to: ",
  file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  ),
  "\n",
  sep = ""
)
