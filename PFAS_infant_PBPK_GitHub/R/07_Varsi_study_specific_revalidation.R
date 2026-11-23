# REVISED EXTERNAL VALIDATION: independent Varsi-specific maternal-to-milk exposure reconstruction.
# Kidney uses tissue/total-plasma P_kd; primary runs S1 only.
# IMPORTANT: Varsi-specific exposure trajectories follow the observed maternal series;
# they do NOT directly reuse the study's infant serum outcomes to set predictions.
# The revised code is derived from historical v1 and REQUIRES execution to refresh results.
# ============================================================
# PFAS infant PBPK model — VARSI 2022 FINAL FREEZE VALIDATION v1
# Uses the frozen final infant model parameterization; no validation-driven tuning.
# Derived directly from final infant PBPK v3-growth
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
# 3) Breast-milk/maternal-serum R_BM values are retained from the Zheng et al.
#    2022 matched maternal-serum/cord-serum/breast-milk transfer study, cited
#    as reference 35 in Liu et al. 2024. Liu 2024 does not itself tabulate R_BM.
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
PROJECT_ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT",getwd()),winslash="/",mustWork=FALSE)
file_input <- Sys.getenv("PFAS_PBPK_INPUT",file.path(PROJECT_ROOT,"data","input-PBPK.xlsx"))
file_abraham_reference <- file.path(PROJECT_ROOT,"data","external","1-s2.0-S0160412024006330-mmc2.xlsx")

BASE_RESULTS_DIR <- file.path(PROJECT_ROOT,"results")
MODEL_FREEZE_VERSION <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
VALIDATION_VERSION <- "Varsi_revised"
VALIDATION_DATE <- "2026-11-03"
OUTPUT_DIR <- file.path(BASE_RESULTS_DIR, MODEL_FREEZE_VERSION, VALIDATION_VERSION)

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

SCENARIOS <- "S1"
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
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

# ============================================================
# 2. LITERATURE-INFORMED MILK TRAJECTORIES
# ============================================================
trajectory_params <- tibble::tribble(
  ~Compound, ~supported_by_blomberg, ~early_change_mean, ~early_ci_low, ~early_ci_high,
  ~monthly_change_mean, ~monthly_ci_low, ~monthly_ci_high, ~source_note,
  "PFOA",  TRUE,  -0.167, -0.280, -0.035, -0.120, -0.220, -0.015,
  "Blomberg 2023 Tables 3–4",
  "PFNA",  FALSE, 0.000, NA, NA, -0.120, -0.220, -0.015,
  "No supported longitudinal trajectory estimate in Blomberg 2023; flat scenario fallback",
  "PFDA",  FALSE, 0.000, NA, NA, -0.120, -0.220, -0.015,
  "No supported longitudinal trajectory estimate in Blomberg 2023; flat scenario fallback",
  "PFUdA", FALSE, 0.000, NA, NA, -0.120, -0.220, -0.015,
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
# R_BM: retained from the existing parameterization attributed to Zheng et al.
# 2022 (Sci Total Environ 812:152446), which is cited as reference 35 by
# Liu et al. 2024 and studied matched maternal serum, cord serum and breast milk.
# Liu et al. 2024 itself does NOT report the six chemical-specific R_BM values,
# so the R_BM provenance is explicitly separated from the directly observed R_CM.
# FINAL-FROZEN transfer parameters.
# PFNA and PFHxS R_BM were updated before model freeze using independent
# literature-pooled milk:serum ratios, following a targeted transfer audit.
# These values are NOT fitted to Varsi infant outcomes.
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
    "Prior v3 chemical-specific value retained",
    "FINAL FREEZE: independent pooled-literature milk:serum ratio",
    "Prior v3 chemical-specific value retained",
    "Prior v3 chemical-specific value retained",
    "FINAL FREEZE: independent pooled-literature milk:serum ratio",
    "Prior v3 chemical-specific value retained"
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

  Freeze_Status = c(
    "FROZEN",
    "FROZEN_UPDATED",
    "FROZEN",
    "FROZEN",
    "FROZEN_UPDATED",
    "FROZEN"
  ),

  Matrix_Compatibility_Note = c(
    "R_CM and R_BM originate from different studies/matrices and are treated as effective transfer factors.",
    "PFNA final R_BM is literature-pooled and frozen before this final Varsi validation.",
    "R_CM and R_BM originate from different studies/matrices and are treated as effective transfer factors.",
    "R_CM and R_BM originate from different studies/matrices and are treated as effective transfer factors.",
    "PFHxS final R_BM is literature-pooled and frozen before this final Varsi validation.",
    "R_CM and R_BM originate from different studies/matrices and are treated as effective transfer factors."
  ),

  stringsAsFactors = FALSE
)

write_xlsx(
  list(Transfer_Parameters = transfer_data),
  path = file.path(OUTPUT_DIR, "tables", "maternal_cord_milk_transfer_parameters_v3.xlsx")
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
    
    dAk <- Q_kd * (Cp - Ck / P_kd)
    
    dAr <- Q_re * (Cp - Cr / P_re)
    
    # ----------------------------------------------------------
    # Plasma mass balance
    # ----------------------------------------------------------
    dAp <-
      Q_li * Cl / P_li +
      Q_kd * Ck / P_kd +
      Q_re * Cr / P_re -
      (Q_li + Q_kd + Q_re) * Cp -
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
# 18. VARSI 2022 STUDY-SPECIFIC EXTERNAL EVALUATION v1
# ============================================================
#
# Source:
#   Varsi K et al. Environmental Research 204 (2022) 112392
#   DOI: 10.1016/j.envres.2021.112392
#
# Purpose:
#   Perform the FINAL study-specific Varsi 2022 external evaluation using the
#   structurally frozen infant PBPK model.
#
#   Reconstruct Varsi-specific maternal -> fetal -> milk exposure using only
#   published maternal serum summaries and frozen transfer/renal parameters,
#   then predict infant plasma at six months.
#
# FINAL-FREEZE RULE:
#   Infant 6-month PFAS measurements are NEVER used to fit, choose, or revise
#   any parameter in this script. All model parameters are frozen before the
#   observed infant outcomes are compared with predictions.
#
# Evaluation class:
#   Study-specific POPULATION-LEVEL external evaluation because published
#   summary trajectories are available rather than individual paired
#   mother-infant longitudinal records.
# ============================================================

VARSI_N <- 114L
VARSI_N_MALE <- 60L
VARSI_N_FEMALE <- 54L

VARSI_GEST_AGE_WK <- 39.9
VARSI_BIRTH_BW_MEAN_KG <- 3.573
VARSI_BIRTH_BW_SD_KG <- 0.418
VARSI_6M_BW_MEAN_KG <- 7.969
VARSI_6M_BW_SD_KG <- 0.987

VARSI_EXCLUSIVE_BF_MONTHS <- 3.8
VARSI_BF_PREVALENCE_6M <- 0.86

DAYS_PER_MONTH <- 30.4375
VARSI_6M_DAYS <- 6 * DAYS_PER_MONTH
VARSI_4M_DAYS <- 4 * DAYS_PER_MONTH
VARSI_6WK_DAYS <- 42
VARSI_EBF_END_DAYS <- VARSI_EXCLUSIVE_BF_MONTHS * DAYS_PER_MONTH

VARSI_TRAJECTORY_LEVELS <- c("Q25", "MEDIAN", "Q75")

VARSI_FEEDING_SCENARIOS <- c(
  "FULL_BF",
  "COHORT_PREVALENCE_APPROX",
  "PARTIAL_50_SENSITIVITY"
)

VARSI_PRIMARY_TRAJECTORY <- "MEDIAN"
VARSI_PRIMARY_FEEDING <- "COHORT_PREVALENCE_APPROX"

VARSI_OUTPUT_DAYS <- sort(unique(c(
  seq(0, floor(VARSI_6M_DAYS), by = 1),
  VARSI_6M_DAYS
)))
VARSI_OUTPUT_H <- VARSI_OUTPUT_DAYS * 24

VARSI_DIR <- OUTPUT_DIR
VARSI_TABLE_DIR <- file.path(VARSI_DIR, "tables")
VARSI_PLOT_DIR <- file.path(VARSI_DIR, "plots")
VARSI_QC_DIR <- file.path(VARSI_DIR, "qc")

for (d in c(VARSI_DIR, VARSI_TABLE_DIR, VARSI_PLOT_DIR, VARSI_QC_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# 19. SAFE HELPERS
# ============================================================

v_safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

v_safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x)
}

v_safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

v_safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

v_safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

v_safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(
    quantile(
      x,
      probs = p,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )
}

v_safe_divide <- function(num, den) {
  ifelse(
    is.finite(num) & is.finite(den) & abs(den) > 1e-30,
    num / den,
    NA_real_
  )
}

v_safe_write_xlsx <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  excel_limit <- 1048576L
  nr <- vapply(
    x,
    function(z) if (is.data.frame(z)) nrow(z) else 0L,
    integer(1)
  )

  too_large <- names(nr)[nr + 1L > excel_limit]

  if (length(too_large)) {
    stop(
      "Excel row limit exceeded in sheet(s): ",
      paste(too_large, collapse = ", ")
    )
  }

  ok <- tryCatch(
    {
      writexl::write_xlsx(x = x, path = path)
      TRUE
    },
    error = function(e) {
      message(
        "write_xlsx failed for: ", path,
        "\nReason: ", conditionMessage(e)
      )
      FALSE
    }
  )

  if (ok) return(path)

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  alt <- file.path(
    dirname(path),
    paste0(
      tools::file_path_sans_ext(basename(path)),
      "_", stamp, ".xlsx"
    )
  )

  writexl::write_xlsx(x = x, path = alt)
  warning("Saved to fallback path:\n", alt)
  alt
}

# ============================================================
# 20. PUBLISHED VARSI DATA
# ============================================================
#
# Serum concentration unit: ng/mL.
# PFUnDA is mapped to PBPK compound PFUdA.
# Sum PFHxS and Sum PFOS are mapped to PFHxS and PFOS.
# ============================================================

varsi_maternal <- tibble::tribble(
  ~Compound, ~Time_Label, ~Time_days_relative_birth, ~Median, ~Q25, ~Q75,

  "PFOA",  "GW18",  (18 - VARSI_GEST_AGE_WK) * 7, 1.15, 0.79, 1.49,
  "PFOA",  "GW28",  (28 - VARSI_GEST_AGE_WK) * 7, 1.04, 0.77, 1.35,
  "PFOA",  "GW36",  (36 - VARSI_GEST_AGE_WK) * 7, 0.96, 0.70, 1.24,
  "PFOA",  "PP6WK", VARSI_6WK_DAYS,             0.93, 0.67, 1.11,
  "PFOA",  "PP4M",  VARSI_4M_DAYS,              0.76, 0.56, 0.96,
  "PFOA",  "PP6M",  VARSI_6M_DAYS,              0.67, 0.50, 0.87,

  "PFNA",  "GW18",  (18 - VARSI_GEST_AGE_WK) * 7, 0.45, 0.37, 0.55,
  "PFNA",  "GW28",  (28 - VARSI_GEST_AGE_WK) * 7, 0.41, 0.33, 0.52,
  "PFNA",  "GW36",  (36 - VARSI_GEST_AGE_WK) * 7, 0.38, 0.31, 0.48,
  "PFNA",  "PP6WK", VARSI_6WK_DAYS,             0.39, 0.30, 0.48,
  "PFNA",  "PP4M",  VARSI_4M_DAYS,              0.36, 0.28, 0.48,
  "PFNA",  "PP6M",  VARSI_6M_DAYS,              0.35, 0.27, 0.46,

  "PFDA",  "GW18",  (18 - VARSI_GEST_AGE_WK) * 7, 0.22, 0.17, 0.29,
  "PFDA",  "GW28",  (28 - VARSI_GEST_AGE_WK) * 7, 0.21, 0.16, 0.26,
  "PFDA",  "GW36",  (36 - VARSI_GEST_AGE_WK) * 7, 0.19, 0.15, 0.26,
  "PFDA",  "PP6WK", VARSI_6WK_DAYS,             0.21, 0.16, 0.25,
  "PFDA",  "PP4M",  VARSI_4M_DAYS,              0.19, 0.15, 0.24,
  "PFDA",  "PP6M",  VARSI_6M_DAYS,              0.20, 0.15, 0.25,

  "PFUdA", "GW18",  (18 - VARSI_GEST_AGE_WK) * 7, 0.26, 0.18, 0.34,
  "PFUdA", "GW28",  (28 - VARSI_GEST_AGE_WK) * 7, 0.26, 0.18, 0.35,
  "PFUdA", "GW36",  (36 - VARSI_GEST_AGE_WK) * 7, 0.24, 0.18, 0.34,
  "PFUdA", "PP6WK", VARSI_6WK_DAYS,             0.22, 0.14, 0.28,
  "PFUdA", "PP4M",  VARSI_4M_DAYS,              0.21, 0.14, 0.30,
  "PFUdA", "PP6M",  VARSI_6M_DAYS,              0.22, 0.14, 0.30,

  "PFHxS", "GW18",  (18 - VARSI_GEST_AGE_WK) * 7, 0.54, 0.39, 0.67,
  "PFHxS", "GW28",  (28 - VARSI_GEST_AGE_WK) * 7, 0.51, 0.37, 0.65,
  "PFHxS", "GW36",  (36 - VARSI_GEST_AGE_WK) * 7, 0.54, 0.38, 0.69,
  "PFHxS", "PP6WK", VARSI_6WK_DAYS,             0.52, 0.38, 0.67,
  "PFHxS", "PP4M",  VARSI_4M_DAYS,              0.44, 0.36, 0.66,
  "PFHxS", "PP6M",  VARSI_6M_DAYS,              0.45, 0.35, 0.66,

  "PFOS",  "GW18",  (18 - VARSI_GEST_AGE_WK) * 7, 4.30, 3.23, 5.94,
  "PFOS",  "GW28",  (28 - VARSI_GEST_AGE_WK) * 7, 4.00, 3.03, 5.58,
  "PFOS",  "GW36",  (36 - VARSI_GEST_AGE_WK) * 7, 3.86, 2.96, 5.21,
  "PFOS",  "PP6WK", VARSI_6WK_DAYS,             3.76, 2.70, 5.00,
  "PFOS",  "PP4M",  VARSI_4M_DAYS,              3.55, 2.64, 4.68,
  "PFOS",  "PP6M",  VARSI_6M_DAYS,              3.53, 2.42, 4.37
)

varsi_infant_observed <- tibble::tribble(
  ~Compound, ~Observed_Median_ng_mL, ~Observed_Q25_ng_mL, ~Observed_Q75_ng_mL,
  "PFOA",  3.27, 2.14, 4.23,
  "PFNA",  0.71, 0.56, 0.91,
  "PFDA",  0.23, 0.16, 0.30,
  "PFUdA", 0.16, 0.10, 0.21,
  "PFHxS", 0.68, 0.53, 0.92,
  "PFOS",  4.71, 3.53, 6.23
) %>%
  mutate(
    N_Infant = 94L,
    Age_days = VARSI_6M_DAYS,
    Source = "Varsi et al. 2022, Environmental Research 204:112392"
  )

# ============================================================
# 21. MATERNAL TRAJECTORY RECONSTRUCTION
# ============================================================

get_varsi_level_value <- function(df, level) {
  if (level == "Q25") return(df$Q25)
  if (level == "MEDIAN") return(df$Median)
  if (level == "Q75") return(df$Q75)
  stop("Unknown level: ", level)
}

log_interp_two <- function(x1, y1, x2, y2, xout) {
  if (
    !all(is.finite(c(x1, y1, x2, y2, xout))) ||
    y1 <= 0 || y2 <= 0 || x1 == x2
  ) return(NA_real_)

  exp(
    log(y1) +
      (xout - x1) / (x2 - x1) *
      (log(y2) - log(y1))
  )
}

get_varsi_birth_maternal_ng_mL <- function(compound, level) {

  d <- varsi_maternal %>%
    filter(
      Compound == compound,
      Time_Label %in% c("GW36", "PP6WK")
    ) %>%
    arrange(Time_days_relative_birth)

  if (nrow(d) != 2L) {
    stop("Birth maternal reconstruction failed for ", compound)
  }

  y <- get_varsi_level_value(d, level)

  log_interp_two(
    d$Time_days_relative_birth[1], y[1],
    d$Time_days_relative_birth[2], y[2],
    0
  )
}

make_varsi_maternal_postpartum_function <- function(compound, level) {

  d <- varsi_maternal %>%
    filter(
      Compound == compound,
      Time_Label %in% c("PP6WK", "PP4M", "PP6M")
    ) %>%
    arrange(Time_days_relative_birth)

  birth_val <- get_varsi_birth_maternal_ng_mL(
    compound,
    level
  )

  y_post <- get_varsi_level_value(d, level)

  x <- c(0, d$Time_days_relative_birth)
  y <- c(birth_val, y_post)

  f <- stats::approxfun(
    x = x,
    y = log(y),
    method = "linear",
    rule = 2
  )

  function(t_days) {
    exp(f(as.numeric(t_days)))
  }
}

# ============================================================
# 22. STUDY-SPECIFIC MILK RECONSTRUCTION
# ============================================================

get_transfer_row <- function(compound) {
  z <- transfer_data %>%
    filter(Compound == compound)

  if (nrow(z) != 1L) {
    stop("Transfer row not uniquely found for ", compound)
  }

  z
}

make_varsi_milk_function <- function(
    compound,
    trajectory_level,
    feeding_scenario
) {

  maternal_fun <- make_varsi_maternal_postpartum_function(
    compound,
    trajectory_level
  )

  tr <- get_transfer_row(compound)
  rbm <- as.numeric(tr$R_BM)

  feeding_fraction <- function(t_days) {

    t <- pmax(as.numeric(t_days), 0)

    if (feeding_scenario == "FULL_BF") {
      return(rep(1, length(t)))
    }

    target_6m <- if (
      feeding_scenario == "COHORT_PREVALENCE_APPROX"
    ) {
      VARSI_BF_PREVALENCE_6M
    } else if (
      feeding_scenario == "PARTIAL_50_SENSITIVITY"
    ) {
      0.50
    } else {
      stop("Unknown feeding scenario: ", feeding_scenario)
    }

    out <- rep(1, length(t))
    after <- t > VARSI_EBF_END_DAYS

    if (any(after)) {
      frac <- (
        t[after] - VARSI_EBF_END_DAYS
      ) / (
        VARSI_6M_DAYS - VARSI_EBF_END_DAYS
      )

      frac <- pmin(pmax(frac, 0), 1)

      out[after] <- 1 - frac * (1 - target_6m)
    }

    out
  }

  list(
    maternal_serum_ng_mL = function(t_days) {
      maternal_fun(t_days)
    },

    reconstructed_milk_ng_L = function(t_days) {
      maternal_fun(t_days) * 1000 * rbm
    },

    feeding_fraction = feeding_fraction,

    effective_milk_ng_L = function(t_days) {
      maternal_fun(t_days) *
        1000 *
        rbm *
        feeding_fraction(t_days)
    }
  )
}

# ============================================================
# 23. VARSI VIRTUAL COHORT
# ============================================================

bw_probs <- (seq_len(VARSI_N) - 0.5) / VARSI_N

bw_virtual <- qnorm(
  bw_probs,
  mean = VARSI_BIRTH_BW_MEAN_KG,
  sd = VARSI_BIRTH_BW_SD_KG
)

bw_virtual <- pmin(
  pmax(bw_virtual, 2.0),
  5.5
)

set.seed(20260908)

sex_vec <- c(
  rep(1, VARSI_N_MALE),
  rep(2, VARSI_N_FEMALE)
)

sex_vec <- sample(
  sex_vec,
  size = length(sex_vec),
  replace = FALSE
)

varsi_virtual <- tibble::tibble(
  ID = seq_len(VARSI_N),
  infsex = sex_vec,
  BW_birth_kg = bw_virtual,
  mix_ratio = 1
)

varsi_virtual_qc <- tibble::tibble(
  Metric = c(
    "N",
    "Male_N",
    "Female_N",
    "Birth_weight_mean_kg",
    "Birth_weight_SD_kg"
  ),
  Published = c(
    VARSI_N,
    VARSI_N_MALE,
    VARSI_N_FEMALE,
    VARSI_BIRTH_BW_MEAN_KG,
    VARSI_BIRTH_BW_SD_KG
  ),
  Virtual = c(
    nrow(varsi_virtual),
    sum(varsi_virtual$infsex == 1),
    sum(varsi_virtual$infsex == 2),
    mean(varsi_virtual$BW_birth_kg),
    sd(varsi_virtual$BW_birth_kg)
  )
) %>%
  mutate(Difference = Virtual - Published)

# ============================================================
# 24. PRECOMPUTE VARSI PHYSIOLOGY
# ============================================================

cat("\n============================================================\n")
cat("PRECOMPUTING VARSI VIRTUAL COHORT PHYSIOLOGY\n")
cat("============================================================\n")

varsi_jobs <- vector("list", VARSI_N)

for (i in seq_len(VARSI_N)) {

  if (
    i == 1L ||
    i %% 20L == 0L ||
    i == VARSI_N
  ) {
    cat(sprintf(
      "  physiology: %d/%d\n",
      i,
      VARSI_N
    ))
  }

  phys <- make_subject_physiology_cache(
    varsi_virtual[i, , drop = FALSE]
  )

  varsi_jobs[[i]] <- list(
    sub = varsi_virtual[i, , drop = FALSE],
    phys = phys
  )
}

varsi_6m_phys <- bind_rows(
  lapply(
    varsi_jobs,
    function(job) {
      a_yr <- VARSI_6M_DAYS / 365.25

      tibble::tibble(
        ID = job$sub$ID,
        Model_6m_BW_kg =
          job$phys$fns$f_BW(a_yr)
      )
    }
  )
)

varsi_physiology_validation <- tibble::tibble(
  Metric = c(
    "Mean_6m_weight_kg",
    "SD_6m_weight_kg"
  ),
  Observed = c(
    VARSI_6M_BW_MEAN_KG,
    VARSI_6M_BW_SD_KG
  ),
  Model = c(
    mean(varsi_6m_phys$Model_6m_BW_kg),
    sd(varsi_6m_phys$Model_6m_BW_kg)
  )
) %>%
  mutate(
    Difference = Model - Observed,
    Percent_Difference = 100 * Difference / Observed
  )

# ============================================================
# 25. DIRECT MATERNAL-SERUM BIRTH INITIALIZATION
# ============================================================

get_varsi_initial_state <- function(
    compound,
    trajectory_level,
    subject_phys
) {

  tr <- get_transfer_row(compound)

  c_mat_birth_ng_mL <-
    get_varsi_birth_maternal_ng_mL(
      compound,
      trajectory_level
    )

  c_fetal_birth_ng_L <-
    c_mat_birth_ng_mL *
    1000 *
    as.numeric(tr$R_CM)

  p <- pfas_params[[compound]]

  c(
    Al =
      c_fetal_birth_ng_L *
      as.numeric(p["P_li"]) *
      subject_phys$fns$f_V_li(0),

    Ak =
      c_fetal_birth_ng_L *
      as.numeric(p["P_kd"]) *
      subject_phys$fns$f_V_kd(0),

    Ar =
      c_fetal_birth_ng_L *
      as.numeric(p["P_re"]) *
      subject_phys$fns$f_V_re(0),

    Ap =
      c_fetal_birth_ng_L *
      subject_phys$fns$f_V_pl(0),

    Afil = 0,
    Aust = 0,
    Aurine = 0,
    Afeces = 0,
    Ainput = 0
  )
}

# ============================================================
# 26. ONE VARSI SIMULATION
# ============================================================

run_varsi_single <- function(
    job,
    compound,
    trajectory_level,
    feeding_scenario
) {

  sub <- job$sub
  subject_phys <- job$phys

  renal_cache <- make_renal_cache_unified(
    compound,
    subject_phys$phys_df
  )

  milk <- make_varsi_milk_function(
    compound,
    trajectory_level,
    feeding_scenario
  )

  p <- pfas_params[[compound]]

  pars <- list(
    P_li = as.numeric(p["P_li"]),
    P_kd = as.numeric(p["P_kd"]),
    P_re = as.numeric(p["P_re"]),
    mix_ratio = 1,

    f_V_pl = subject_phys$fns$f_V_pl,
    f_V_li = subject_phys$fns$f_V_li,
    f_V_kd = subject_phys$fns$f_V_kd,
    f_V_re = subject_phys$fns$f_V_re,

    f_Q_li = subject_phys$fns$f_Q_li,
    f_Q_kd = subject_phys$fns$f_Q_kd,
    f_Q_re = subject_phys$fns$f_Q_re,
    f_KMilk = subject_phys$fns$f_KMilk,

    f_Free = renal_cache$f_Free,
    f_Tm = renal_cache$f_Tm,
    f_Kt = renal_cache$f_Kt,
    f_Qfil = renal_cache$f_Qfil,
    f_kust = renal_cache$f_kust,
    f_ku = renal_cache$f_ku,
    f_V_fil = renal_cache$f_Vfil,

    f_C_milk = function(t_days) {
      milk$effective_milk_ng_L(t_days)
    }
  )

  init <- get_varsi_initial_state(
    compound,
    trajectory_level,
    subject_phys
  )

  out <- tryCatch(
    deSolve::ode(
      y = init,
      times = VARSI_OUTPUT_H,
      func = pbpk_ode_unified_tmkt,
      parms = pars,
      method = "lsoda",
      rtol = ODE_RTOL,
      atol = ODE_ATOL
    ),
    error = function(e) {
      warning(
        "Varsi ODE failed: ID=", sub$ID,
        ", compound=", compound,
        ", trajectory=", trajectory_level,
        ", feeding=", feeding_scenario,
        ": ", conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(out)) return(NULL)

  d <- as.data.frame(out)

  names(d)[1:10] <- c(
    "time",
    "Al", "Ak", "Ar", "Ap",
    "Afil", "Aust",
    "Aurine", "Afeces",
    "Ainput"
  )

  t_days <- d$time / 24
  t_yr <- d$time / 8760

  Vpl <- pmax(
    subject_phys$fns$f_V_pl(t_yr),
    MIN_VOLUME_L
  )

  Vfil <- pmax(
    renal_cache$f_Vfil(t_yr),
    MIN_VOLUME_L
  )

  BW <- subject_phys$fns$f_BW(t_yr)
  Qfil <- pmax(renal_cache$f_Qfil(t_yr), 0)
  Free <- pmax(renal_cache$f_Free(t_yr), 0)
  Tm <- pmax(renal_cache$f_Tm(t_yr), 0)
  Kt <- pmax(renal_cache$f_Kt(t_yr), 1e-18)

  C_Plasma <- d$Ap / Vpl
  C_Filtrate <- d$Afil / Vfil

  R_Filtration <- Qfil * Free * pmax(C_Plasma, 0)

  R_Reabsorption <- ifelse(
    C_Filtrate > 0,
    Tm * C_Filtrate / (Kt + C_Filtrate),
    0
  )

  maternal_serum <- milk$maternal_serum_ng_mL(t_days)
  milk_raw <- milk$reconstructed_milk_ng_L(t_days)
  feeding_fraction <- milk$feeding_fraction(t_days)
  milk_effective <- milk$effective_milk_ng_L(t_days)

  KMilk <- pmax(subject_phys$fns$f_KMilk(t_yr), 0)
  Ingestion_Rate <- KMilk * milk_effective

  A_Systemic <- d$Al + d$Ak + d$Ar + d$Ap

  A_Total <-
    A_Systemic +
    d$Afil +
    d$Aust +
    d$Aurine +
    d$Afeces

  initial_total <- A_Total[1]
  Mass_Balance_Expected <- initial_total + d$Ainput
  Mass_Balance_Error <- A_Total - Mass_Balance_Expected

  Mass_Balance_Relative_Error <-
    Mass_Balance_Error /
    pmax(abs(Mass_Balance_Expected), 1e-18)

  tibble::tibble(
    ID = sub$ID,
    Compound = compound,
    Trajectory_Level = trajectory_level,
    Feeding_Scenario = feeding_scenario,

    Primary_Evaluation =
      trajectory_level == VARSI_PRIMARY_TRAJECTORY &
      feeding_scenario == VARSI_PRIMARY_FEEDING,

    Time_days = t_days,
    BW_current_kg = BW,

    Maternal_Serum_ng_mL = maternal_serum,
    Milk_Reconstructed_ng_L = milk_raw,
    Feeding_Fraction = feeding_fraction,
    Milk_Effective_ng_L = milk_effective,
    Ingestion_Rate_ng_h = Ingestion_Rate,
    A_Input_Cum_ng = d$Ainput,

    C_Plasma_ng_L = C_Plasma,
    C_Plasma_ng_mL = C_Plasma / 1000,
    A_Systemic_Burden_ng = A_Systemic,

    GFR_L_h = Qfil,
    R_Filtration_ng_h = R_Filtration,
    R_Reabsorption_ng_h = R_Reabsorption,

    Reabsorption_Rate_Ratio =
      ifelse(
        R_Filtration > 0,
        R_Reabsorption / R_Filtration,
        NA_real_
      ),

    Mass_Balance_Relative_Error =
      Mass_Balance_Relative_Error
  )
}

# ============================================================
# 27. RUN ALL VARSI SIMULATIONS
# ============================================================

run_varsi_subject <- function(job) {

  out <- list()
  k <- 1L

  for (compound in pfas_cols) {
    for (trajectory_level in VARSI_TRAJECTORY_LEVELS) {
      for (feeding_scenario in VARSI_FEEDING_SCENARIOS) {

        z <- run_varsi_single(
          job = job,
          compound = compound,
          trajectory_level = trajectory_level,
          feeding_scenario = feeding_scenario
        )

        if (!is.null(z)) {
          out[[k]] <- z
          k <- k + 1L
        }
      }
    }
  }

  if (!length(out)) return(tibble::tibble())
  bind_rows(out)
}

expected_runs <-
  length(varsi_jobs) *
  length(pfas_cols) *
  length(VARSI_TRAJECTORY_LEVELS) *
  length(VARSI_FEEDING_SCENARIOS)

cat("\n============================================================\n")
cat("RUNNING VARSI 2022 STUDY-SPECIFIC EXTERNAL EVALUATION\n")
cat("============================================================\n")
cat("Virtual infants: ", length(varsi_jobs), "\n", sep = "")
cat("PFAS: ", length(pfas_cols), "\n", sep = "")
cat("Expected ODE simulations: ", expected_runs, "\n", sep = "")

v_start <- Sys.time()

if (
  USE_PARALLEL &&
  length(varsi_jobs) > 1L &&
  PARALLEL_WORKERS > 1L
) {

  workers <- min(PARALLEL_WORKERS, length(varsi_jobs))

  cat("Parallel execution with ", workers, " workers...\n", sep = "")

  old_plan <- future::plan()

  future::plan(
    future::multisession,
    workers = workers
  )

  varsi_subject_results <- tryCatch(
    future.apply::future_lapply(
      varsi_jobs,
      run_varsi_subject,
      future.seed = TRUE,
      future.packages = c(
        "deSolve",
        "dplyr",
        "tibble"
      )
    ),
    finally = future::plan(old_plan)
  )

} else {

  varsi_subject_results <- vector(
    "list",
    length(varsi_jobs)
  )

  for (i in seq_along(varsi_jobs)) {

    if (
      i == 1L ||
      i %% 10L == 0L ||
      i == length(varsi_jobs)
    ) {
      cat(sprintf(
        "  subjects: %d/%d (%.1f%%)\n",
        i,
        length(varsi_jobs),
        100 * i / length(varsi_jobs)
      ))
    }

    varsi_subject_results[[i]] <-
      run_varsi_subject(
        varsi_jobs[[i]]
      )
  }
}

varsi_all <- bind_rows(varsi_subject_results)

if (!nrow(varsi_all)) {
  stop("No Varsi study-specific simulations were generated.")
}

v_elapsed <- as.numeric(
  difftime(Sys.time(), v_start, units = "secs")
)

cat(sprintf(
  "Varsi simulations completed in %.2f sec (%.2f min).\n",
  v_elapsed,
  v_elapsed / 60
))

# ============================================================
# 28. QC
# ============================================================

varsi_run_qc <- varsi_all %>%
  group_by(
    Compound,
    Trajectory_Level,
    Feeding_Scenario
  ) %>%
  summarise(
    N_ID = n_distinct(ID),
    N_Timepoints = n_distinct(Time_days),

    Max_Abs_Mass_Balance_Relative_Error =
      v_safe_max(
        abs(Mass_Balance_Relative_Error)
      ),

    Min_C_Plasma_ng_L =
      v_safe_min(C_Plasma_ng_L),

    Max_Reabsorption_Rate_Ratio =
      v_safe_max(Reabsorption_Rate_Ratio),

    .groups = "drop"
  ) %>%
  mutate(
    Missing_ID_Flag = N_ID != VARSI_N,

    Missing_Timepoint_Flag =
      N_Timepoints != length(VARSI_OUTPUT_DAYS),

    Mass_Balance_Flag =
      Max_Abs_Mass_Balance_Relative_Error >
      MASS_BALANCE_REL_TOL,

    Negative_Plasma_Flag =
      Min_C_Plasma_ng_L < -1e-12
  )

# ============================================================
# 29. POPULATION TIMECOURSES
# ============================================================

varsi_curve <- varsi_all %>%
  group_by(
    Compound,
    Trajectory_Level,
    Feeding_Scenario,
    Primary_Evaluation,
    Time_days
  ) %>%
  summarise(
    N_ID = n_distinct(ID),

    Median_C_Plasma_ng_mL =
      v_safe_median(C_Plasma_ng_mL),

    P025_C_Plasma_ng_mL =
      v_safe_quantile(C_Plasma_ng_mL, 0.025),

    P975_C_Plasma_ng_mL =
      v_safe_quantile(C_Plasma_ng_mL, 0.975),

    Median_Maternal_Serum_ng_mL =
      v_safe_median(Maternal_Serum_ng_mL),

    Median_Reconstructed_Milk_ng_L =
      v_safe_median(Milk_Reconstructed_ng_L),

    Median_Feeding_Fraction =
      v_safe_median(Feeding_Fraction),

    Median_Effective_Milk_ng_L =
      v_safe_median(Milk_Effective_ng_L),

    Median_A_Input_Cum_ng =
      v_safe_median(A_Input_Cum_ng),

    Median_BW_current_kg =
      v_safe_median(BW_current_kg),

    .groups = "drop"
  )

# ============================================================
# 30. SIX-MONTH EVALUATION
# ============================================================

varsi_6m <- varsi_all %>%
  filter(
    abs(Time_days - VARSI_6M_DAYS) < 1e-8
  ) %>%
  group_by(
    Compound,
    Trajectory_Level,
    Feeding_Scenario,
    Primary_Evaluation
  ) %>%
  summarise(
    N_ID = n_distinct(ID),

    Predicted_Median_ng_mL =
      v_safe_median(C_Plasma_ng_mL),

    Predicted_P025_ng_mL =
      v_safe_quantile(C_Plasma_ng_mL, 0.025),

    Predicted_P975_ng_mL =
      v_safe_quantile(C_Plasma_ng_mL, 0.975),

    Predicted_Mean_ng_mL =
      v_safe_mean(C_Plasma_ng_mL),

    Predicted_SD_ng_mL =
      v_safe_sd(C_Plasma_ng_mL),

    .groups = "drop"
  ) %>%
  left_join(
    varsi_infant_observed,
    by = "Compound"
  ) %>%
  mutate(
    Predicted_to_Observed_Ratio =
      v_safe_divide(
        Predicted_Median_ng_mL,
        Observed_Median_ng_mL
      ),

    Fold_Error =
      ifelse(
        is.finite(Predicted_to_Observed_Ratio) &
          Predicted_to_Observed_Ratio > 0,
        exp(abs(log(Predicted_to_Observed_Ratio))),
        NA_real_
      ),

    Within_2Fold =
      Predicted_to_Observed_Ratio >= 0.5 &
      Predicted_to_Observed_Ratio <= 2,

    Observed_Median_Within_Physiology_95PI =
      Observed_Median_ng_mL >= Predicted_P025_ng_mL &
      Observed_Median_ng_mL <= Predicted_P975_ng_mL,

    Evaluation_Class =
      "Study-specific population-level external evaluation"
  )

# ============================================================
# 31. EXPOSURE-IQR ENVELOPE
# ============================================================

varsi_exposure_envelope <- varsi_6m %>%
  filter(
    Feeding_Scenario == VARSI_PRIMARY_FEEDING
  ) %>%
  group_by(Compound) %>%
  summarise(
    Envelope_Lower_ng_mL =
      v_safe_min(
        Predicted_P025_ng_mL[
          Trajectory_Level == "Q25"
        ]
      ),

    Median_Trajectory_Prediction_ng_mL =
      v_safe_median(
        Predicted_Median_ng_mL[
          Trajectory_Level == "MEDIAN"
        ]
      ),

    Envelope_Upper_ng_mL =
      v_safe_max(
        Predicted_P975_ng_mL[
          Trajectory_Level == "Q75"
        ]
      ),

    .groups = "drop"
  ) %>%
  left_join(
    varsi_infant_observed,
    by = "Compound"
  ) %>%
  mutate(
    Observed_Median_Inside_Exposure_Envelope =
      Observed_Median_ng_mL >= Envelope_Lower_ng_mL &
      Observed_Median_ng_mL <= Envelope_Upper_ng_mL,

    Observed_IQR_Overlap =
      Observed_Q75_ng_mL >= Envelope_Lower_ng_mL &
      Observed_Q25_ng_mL <= Envelope_Upper_ng_mL
  )

# ============================================================
# 32. METRICS
# ============================================================

calc_varsi_metrics <- function(df) {

  d <- df %>%
    filter(
      is.finite(Predicted_Median_ng_mL),
      Predicted_Median_ng_mL > 0,
      is.finite(Observed_Median_ng_mL),
      Observed_Median_ng_mL > 0
    ) %>%
    mutate(
      ratio =
        Predicted_Median_ng_mL /
        Observed_Median_ng_mL,
      log_ratio = log(ratio)
    )

  if (!nrow(d)) {
    return(tibble::tibble(
      N = 0L,
      MFE = NA_real_,
      AFE = NA_real_,
      RMSE_ng_mL = NA_real_,
      R2_linear = NA_real_,
      R2_log10 = NA_real_,
      TwoFold_Accuracy = NA_real_
    ))
  }

  r2_linear <- if (nrow(d) >= 2L) {
    cor(
      d$Predicted_Median_ng_mL,
      d$Observed_Median_ng_mL
    )^2
  } else NA_real_

  r2_log10 <- if (nrow(d) >= 2L) {
    cor(
      log10(d$Predicted_Median_ng_mL),
      log10(d$Observed_Median_ng_mL)
    )^2
  } else NA_real_

  tibble::tibble(
    N = nrow(d),

    MFE =
      exp(mean(d$log_ratio)),

    AFE =
      exp(mean(abs(d$log_ratio))),

    RMSE_ng_mL =
      sqrt(
        mean(
          (
            d$Predicted_Median_ng_mL -
              d$Observed_Median_ng_mL
          )^2
        )
      ),

    R2_linear = r2_linear,
    R2_log10 = r2_log10,

    TwoFold_Accuracy =
      mean(
        d$ratio >= 0.5 &
        d$ratio <= 2
      )
  )
}

varsi_metrics <- varsi_6m %>%
  group_by(
    Trajectory_Level,
    Feeding_Scenario,
    Primary_Evaluation
  ) %>%
  group_modify(
    ~ calc_varsi_metrics(.x)
  ) %>%
  ungroup()

primary_compound <- varsi_6m %>%
  filter(Primary_Evaluation)

primary_metrics <- varsi_metrics %>%
  filter(Primary_Evaluation)

# ============================================================
# 33. BIRTH / MILK EXPOSURE AUDIT
# ============================================================

varsi_birth_audit <- tidyr::crossing(
  Compound = pfas_cols,
  Trajectory_Level = VARSI_TRAJECTORY_LEVELS
) %>%
  rowwise() %>%
  mutate(
    Maternal_Birth_ng_mL =
      get_varsi_birth_maternal_ng_mL(
        Compound,
        Trajectory_Level
      ),

    R_CM =
      as.numeric(
        get_transfer_row(Compound)$R_CM
      ),

    R_BM =
      as.numeric(
        get_transfer_row(Compound)$R_BM
      ),

    Predicted_Fetal_Plasma_Birth_ng_mL =
      Maternal_Birth_ng_mL * R_CM,

    Reconstructed_Milk_Birth_ng_mL =
      Maternal_Birth_ng_mL * R_BM
  ) %>%
  ungroup()

varsi_milk_nodes <- tidyr::crossing(
  Compound = pfas_cols,
  Trajectory_Level = VARSI_TRAJECTORY_LEVELS,
  Feeding_Scenario = VARSI_FEEDING_SCENARIOS,
  Time_days = c(
    0,
    VARSI_6WK_DAYS,
    VARSI_4M_DAYS,
    VARSI_6M_DAYS
  )
) %>%
  rowwise() %>%
  mutate(
    Maternal_Serum_ng_mL = {
      f <- make_varsi_milk_function(
        Compound,
        Trajectory_Level,
        Feeding_Scenario
      )
      f$maternal_serum_ng_mL(Time_days)
    },

    Reconstructed_Milk_ng_L = {
      f <- make_varsi_milk_function(
        Compound,
        Trajectory_Level,
        Feeding_Scenario
      )
      f$reconstructed_milk_ng_L(Time_days)
    },

    Feeding_Fraction = {
      f <- make_varsi_milk_function(
        Compound,
        Trajectory_Level,
        Feeding_Scenario
      )
      f$feeding_fraction(Time_days)
    },

    Effective_Milk_ng_L = {
      f <- make_varsi_milk_function(
        Compound,
        Trajectory_Level,
        Feeding_Scenario
      )
      f$effective_milk_ng_L(Time_days)
    }
  ) %>%
  ungroup()

# ============================================================
# 34. PLOTS
# ============================================================

p_mat <- ggplot(
  varsi_maternal %>%
    filter(
      Time_Label %in% c(
        "GW36",
        "PP6WK",
        "PP4M",
        "PP6M"
      )
    ),
  aes(
    x = Time_days_relative_birth,
    y = Median
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  geom_vline(
    xintercept = 0,
    linetype = 3
  ) +
  facet_wrap(
    ~Compound,
    scales = "free_y"
  ) +
  labs(
    title = "Varsi 2022 maternal serum PFAS trajectory",
    subtitle = "Published medians; delivery reconstructed between GW36 and postpartum week 6",
    x = "Days relative to delivery",
    y = "Maternal serum PFAS (ng/mL)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(
    VARSI_PLOT_DIR,
    "01_varsi_maternal_serum_trajectory_v1.png"
  ),
  p_mat,
  width = 11,
  height = 7,
  dpi = PLOT_DPI
)

primary_curve <- varsi_curve %>%
  filter(Primary_Evaluation)

p_inf <- ggplot(
  primary_curve,
  aes(
    x = Time_days,
    y = Median_C_Plasma_ng_mL
  )
) +
  geom_ribbon(
    aes(
      ymin = P025_C_Plasma_ng_mL,
      ymax = P975_C_Plasma_ng_mL
    ),
    alpha = 0.15
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(
    data = varsi_infant_observed,
    aes(
      x = Age_days,
      y = Observed_Median_ng_mL
    ),
    inherit.aes = FALSE,
    size = 2.8
  ) +
  facet_wrap(
    ~Compound,
    scales = "free_y"
  ) +
  labs(
    title = "Varsi study-specific infant PFAS prediction",
    subtitle = "Primary: maternal median trajectory + cohort-prevalence breastfeeding approximation",
    x = "Infant age (days)",
    y = "Infant plasma/serum PFAS (ng/mL)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(
    VARSI_PLOT_DIR,
    "02_varsi_primary_infant_prediction_v1.png"
  ),
  p_inf,
  width = 11,
  height = 7,
  dpi = PLOT_DPI
)

p_pvo <- ggplot(
  primary_compound,
  aes(
    x = Observed_Median_ng_mL,
    y = Predicted_Median_ng_mL,
    label = Compound
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2
  ) +
  geom_abline(
    intercept = 0,
    slope = 2,
    linetype = 3
  ) +
  geom_abline(
    intercept = 0,
    slope = 0.5,
    linetype = 3
  ) +
  geom_point(size = 3) +
  geom_text(
    nudge_x = 0.05,
    nudge_y = 0.05,
    check_overlap = TRUE
  ) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal() +
  labs(
    title = "Varsi 6-month predicted vs observed",
    subtitle = "1:1 and two-fold reference lines",
    x = "Observed infant serum median (ng/mL)",
    y = "Predicted infant plasma median (ng/mL)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    VARSI_PLOT_DIR,
    "03_varsi_predicted_vs_observed_v1.png"
  ),
  p_pvo,
  width = 7,
  height = 7,
  dpi = PLOT_DPI
)

p_feed <- varsi_6m %>%
  filter(
    Trajectory_Level == "MEDIAN"
  ) %>%
  ggplot(
    aes(
      x = Feeding_Scenario,
      y = Predicted_to_Observed_Ratio
    )
  ) +
  geom_hline(
    yintercept = 1,
    linetype = 2
  ) +
  geom_hline(
    yintercept = c(0.5, 2),
    linetype = 3
  ) +
  geom_point(size = 2.5) +
  facet_wrap(
    ~Compound,
    scales = "free_y"
  ) +
  labs(
    title = "Sensitivity to breastfeeding assumption",
    subtitle = "Maternal median trajectory held fixed",
    x = NULL,
    y = "Predicted / observed infant median"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(
    VARSI_PLOT_DIR,
    "04_varsi_feeding_sensitivity_v1.png"
  ),
  p_feed,
  width = 11,
  height = 7,
  dpi = PLOT_DPI
)

# ============================================================
# 35. OUTPUT WORKBOOKS
# ============================================================

varsi_notes <- tibble::tribble(
  ~Item, ~Definition,

  "Evaluation class",
  "Study-specific population-level external evaluation using published Varsi maternal serum summaries. Not subject-level validation.",

  "No outcome fitting",
  "Infant 6-month PFAS data are used only after simulation for evaluation; no PBPK parameter is fitted to them.",

  "Birth initialization",
  "Maternal serum at delivery is log-linearly interpolated between gestational week 36 and postpartum week 6; fetal plasma at birth equals maternal delivery serum × frozen R_CM.",

  "Postnatal exposure",
  "Maternal postpartum serum is piecewise log-linearly interpolated through delivery, postpartum week 6, month 4 and month 6; reconstructed milk equals maternal serum × frozen R_BM.",

  "Primary feeding assumption",
  "Full model milk intake through mean exclusive breastfeeding duration 3.8 months, then linear population weighting to 0.86 at six months because 86% remained breastfed. This is an approximation, not measured individual milk volume.",

  "FULL_BF",
  "Sensitivity scenario with full modeled milk intake through six months.",

  "PARTIAL_50_SENSITIVITY",
  "Sensitivity scenario in which the post-3.8-month feeding multiplier declines linearly to 0.50 at six months; not fitted.",

  "Q25/MEDIAN/Q75",
  "Published maternal serum IQR trajectories used as exposure-distribution sensitivity scenarios; they are not matched longitudinal individuals.",

  "PFUnDA mapping",
  "Varsi PFUnDA is mapped to PBPK PFUdA/PFUnDA.",

  "Two-fold",
  "Reported descriptively; not used as a hard model-acceptance rule."
)

main_xlsx <- file.path(
  VARSI_TABLE_DIR,
  "Varsi2022_study_specific_validation_summary_v1.xlsx"
)

v_safe_write_xlsx(
  list(
    Notes = varsi_notes,
    Primary_6m_Evaluation = primary_compound,
    Primary_Metrics = primary_metrics,
    All_6m_Evaluations = varsi_6m,
    Metrics_All_Scenarios = varsi_metrics,
    Exposure_IQR_Envelope = varsi_exposure_envelope,
    Maternal_Published = varsi_maternal,
    Infant_Observed = varsi_infant_observed,
    Birth_Exposure_Audit = varsi_birth_audit,
    Milk_Reconstruction = varsi_milk_nodes,
    Virtual_Cohort_QC = varsi_virtual_qc,
    Physiology_6m_Validation = varsi_physiology_validation,
    Run_QC = varsi_run_qc
  ),
  main_xlsx
)

v_safe_write_xlsx(
  list(
    Population_Timecourses = varsi_curve
  ),
  file.path(
    VARSI_TABLE_DIR,
    "Varsi2022_study_specific_timecourses_v1.xlsx"
  )
)

v_safe_write_xlsx(
  list(
    Run_QC = varsi_run_qc,
    Virtual_Cohort_QC = varsi_virtual_qc,
    Physiology_6m_Validation = varsi_physiology_validation
  ),
  file.path(
    VARSI_QC_DIR,
    "Varsi2022_study_specific_QC_v1.xlsx"
  )
)

# ============================================================
# 36. CONSOLE SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("VARSI PHYSIOLOGY CHECK\n")
cat("============================================================\n")
print(varsi_physiology_validation, n = Inf)

cat("\n============================================================\n")
cat("PRIMARY 6-MONTH EXTERNAL EVALUATION\n")
cat("============================================================\n")

print(
  primary_compound %>%
    select(
      Compound,
      Predicted_Median_ng_mL,
      Observed_Median_ng_mL,
      Predicted_to_Observed_Ratio,
      Fold_Error,
      Within_2Fold,
      Observed_Median_Within_Physiology_95PI
    ),
  n = Inf
)

cat("\nPrimary metrics:\n")
print(primary_metrics, n = Inf)

cat("\n============================================================\n")
cat("EXPOSURE-IQR ENVELOPE\n")
cat("============================================================\n")
print(varsi_exposure_envelope, n = Inf)

cat("\n============================================================\n")
cat("RUN QC\n")
cat("============================================================\n")
print(varsi_run_qc, n = Inf)

cat("\n============================================================\n")
cat("VARSI 2022 STUDY-SPECIFIC EXTERNAL EVALUATION COMPLETED\n")
cat("============================================================\n")
cat("Main workbook:\n", main_xlsx, "\n", sep = "")
cat(
  "Timecourse workbook:\n",
  file.path(
    VARSI_TABLE_DIR,
    "Varsi2022_study_specific_timecourses_v1.xlsx"
  ),
  "\n",
  sep = ""
)
cat("Plots:\n", VARSI_PLOT_DIR, "\n", sep = "")
cat("QC:\n", VARSI_QC_DIR, "\n", sep = "")
cat("\nDone.\n")


# ============================================================
# 37. FINAL-FREEZE MANUSCRIPT-READY VALIDATION SUMMARY
# ============================================================

final_validation_compound <- primary_compound %>%
  transmute(
    Model_Version = MODEL_FREEZE_VERSION,
    Validation_Version = VALIDATION_VERSION,
    Compound,

    Predicted_Median_ng_mL,
    Predicted_P025_ng_mL,
    Predicted_P975_ng_mL,

    Observed_Median_ng_mL,
    Observed_Q25_ng_mL,
    Observed_Q75_ng_mL,

    Predicted_to_Observed_Ratio,
    Fold_Error,
    Within_2Fold,

    Observed_Median_Within_Physiology_95PI,

    Interpretation = dplyr::case_when(
      Predicted_to_Observed_Ratio < 0.5 ~
        "More than two-fold underprediction",

      Predicted_to_Observed_Ratio > 2 ~
        "More than two-fold overprediction",

      TRUE ~
        "Within two-fold"
    )
  )

final_validation_metrics <- primary_metrics %>%
  mutate(
    Model_Version = MODEL_FREEZE_VERSION,
    Validation_Version = VALIDATION_VERSION,
    Validation_Date = VALIDATION_DATE,
    Evaluation_Class =
      "Study-specific population-level external evaluation; no outcome fitting"
  ) %>%
  select(
    Model_Version,
    Validation_Version,
    Validation_Date,
    Evaluation_Class,
    everything()
  )

# Exposure-envelope summary uses frozen R_BM and the pre-specified primary
# breastfeeding approximation. This is NOT a formal 95% prediction interval.
final_exposure_envelope <- varsi_exposure_envelope %>%
  mutate(
    Model_Version = MODEL_FREEZE_VERSION,
    Validation_Version = VALIDATION_VERSION,
    Interval_Type =
      "Maternal-exposure IQR envelope; not formal prediction interval"
  ) %>%
  select(
    Model_Version,
    Validation_Version,
    Compound,
    Envelope_Lower_ng_mL,
    Median_Trajectory_Prediction_ng_mL,
    Envelope_Upper_ng_mL,
    Observed_Median_ng_mL,
    Observed_Q25_ng_mL,
    Observed_Q75_ng_mL,
    Observed_Median_Inside_Exposure_Envelope,
    Observed_IQR_Overlap,
    Interval_Type,
    everything()
  )

# Final validation QC gate.
final_validation_qc_gate <- tibble::tibble(
  Check = c(
    "All simulation groups have expected virtual-infant count",
    "All simulation groups have expected timepoint count",
    "All mass-balance checks pass",
    "No negative plasma concentrations",
    "Final PFNA R_BM frozen",
    "Final PFHxS R_BM frozen"
  ),

  Pass = c(
    !any(varsi_run_qc$Missing_ID_Flag, na.rm = TRUE),
    !any(varsi_run_qc$Missing_Timepoint_Flag, na.rm = TRUE),
    !any(varsi_run_qc$Mass_Balance_Flag, na.rm = TRUE),
    !any(varsi_run_qc$Negative_Plasma_Flag, na.rm = TRUE),
    abs(
      transfer_data$R_BM[
        transfer_data$Compound == "PFNA"
      ] - 0.0300
    ) < 1e-12,
    abs(
      transfer_data$R_BM[
        transfer_data$Compound == "PFHxS"
      ] - 0.0133
    ) < 1e-12
  ),

  Detail = c(
    paste0(
      "Expected N=", VARSI_N,
      " in every Compound × Exposure × Feeding group"
    ),
    paste0(
      "Expected timepoints=", length(VARSI_OUTPUT_DAYS)
    ),
    paste0(
      "Tolerance=", MASS_BALANCE_REL_TOL
    ),
    "Minimum simulated plasma must be >= -1e-12 ng/L",
    "PFNA R_BM = 0.0300",
    "PFHxS R_BM = 0.0133"
  )
) %>%
  mutate(
    Overall_Final_Validation_QC_Pass =
      all(Pass)
  )

# Frozen parameter table for direct manuscript/Supplement use.
final_frozen_transfer_table <- transfer_data %>%
  transmute(
    Model_Version = MODEL_FREEZE_VERSION,
    Compound,
    R_CM,
    R_BM,
    R_CM_Source,
    R_BM_Source,
    R_BM_Evidence_Class,
    R_CM_Matrix,
    R_BM_Matrix,
    Freeze_Status
  )

final_notes <- tibble::tribble(
  ~Item, ~Definition,

  "Model status",
  "Revised S1 kidney exchange total-plasma model; independently evaluated with Varsi cohort-specific exposures. Not the old frozen v1.1 predictions.",

  "PFNA final R_BM",
  "0.0300; independent pooled-literature milk:serum ratio selected before this final validation.",

  "PFHxS final R_BM",
  "0.0133; independent pooled-literature milk:serum ratio selected before this final validation.",

  "Birth initialization",
  "Varsi maternal delivery serum is reconstructed independently and multiplied by frozen R_CM. Therefore final R_BM affects postnatal milk exposure but not Varsi prenatal initialization.",

  "Postnatal exposure",
  "Maternal postpartum serum is reconstructed independently from Varsi serial maternal measurements; milk PFAS equals maternal serum × frozen R_BM.",

  "Primary feeding scenario",
  "Full model milk intake through mean exclusive breastfeeding duration of 3.8 months, followed by linear population weighting to 0.86 at six months. This is a cohort-level approximation rather than measured individual milk volume.",

  "Observed infant outcome use",
  "Observed infant six-month serum concentrations are used only after simulation to calculate external-evaluation metrics.",

  "Prediction interval caution",
  "The physiology-only 2.5th-97.5th percentile interval reflects virtual-cohort physiology heterogeneity and is not a full parameter/exposure uncertainty interval.",

  "Exposure envelope",
  "Maternal Q25/median/Q75 exposure trajectories provide a separate exposure-IQR envelope and are not a formal 95% prediction interval.",

  "Acceptance language",
  "Two-fold accuracy is reported descriptively. Model qualification is judged from the complete evidence chain rather than a single hard acceptance threshold."
)

final_xlsx <- file.path(
  VARSI_TABLE_DIR,
  "Varsi2022_final_freeze_validation_summary_v1.xlsx"
)

v_safe_write_xlsx(
  list(
    Final_Notes = final_notes,
    Final_Primary_Validation = final_validation_compound,
    Final_Metrics = final_validation_metrics,
    Final_Exposure_Envelope = final_exposure_envelope,
    Final_QC_Gate = final_validation_qc_gate,
    Frozen_Transfer_Parameters = final_frozen_transfer_table,
    Physiology_6m_Check = varsi_physiology_validation,
    Full_Run_QC = varsi_run_qc,
    Maternal_Published = varsi_maternal,
    Infant_Observed = varsi_infant_observed
  ),
  final_xlsx
)

# ------------------------------------------------------------
# Manuscript-ready predicted vs observed figure
# ------------------------------------------------------------

p_final_pvo <- ggplot(
  final_validation_compound,
  aes(
    x = Observed_Median_ng_mL,
    y = Predicted_Median_ng_mL,
    label = Compound
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.8
  ) +
  geom_abline(
    intercept = 0,
    slope = 2,
    linetype = 2,
    linewidth = 0.6
  ) +
  geom_abline(
    intercept = 0,
    slope = 0.5,
    linetype = 2,
    linewidth = 0.6
  ) +
  geom_point(size = 3.2) +
  geom_text(
    nudge_x = 0.04,
    nudge_y = 0.04,
    check_overlap = TRUE
  ) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal() +
  labs(
    title = "Final frozen-model external evaluation: Varsi 2022",
    subtitle = "Study-specific exposure reconstruction; parameters frozen before infant outcome comparison",
    x = "Observed infant serum median at 6 months (ng/mL)",
    y = "Predicted infant plasma median at 6 months (ng/mL)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    VARSI_PLOT_DIR,
    "05_final_freeze_predicted_vs_observed_v1.png"
  ),
  p_final_pvo,
  width = 7.5,
  height = 7,
  dpi = PLOT_DPI
)

ggsave(
  file.path(
    VARSI_PLOT_DIR,
    "05_final_freeze_predicted_vs_observed_v1.pdf"
  ),
  p_final_pvo,
  width = 7.5,
  height = 7
)

# ------------------------------------------------------------
# Final console gate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL FROZEN-MODEL VARSI VALIDATION\n")
cat("============================================================\n")

print(
  final_validation_compound %>%
    select(
      Compound,
      Predicted_Median_ng_mL,
      Observed_Median_ng_mL,
      Predicted_to_Observed_Ratio,
      Fold_Error,
      Within_2Fold,
      Interpretation
    ),
  n = Inf
)

cat("\nFinal metrics:\n")
print(final_validation_metrics, n = Inf)

cat("\nFinal QC gate:\n")
print(final_validation_qc_gate, n = Inf)

if (!all(final_validation_qc_gate$Pass)) {
  warning(
    "FINAL VALIDATION QC GATE DID NOT FULLY PASS. ",
    "Do not interpret external-evaluation metrics until QC is resolved."
  )
} else {
  cat(
    "\nFINAL VALIDATION QC GATE: PASS\n"
  )
}

cat("\nFinal validation workbook:\n")
cat(final_xlsx, "\n", sep = "")

cat(
  "\nMODEL-FREEZE RULE: validation outcomes above must not be used ",
  "to retune the frozen model.\n",
  sep = ""
)

cat("\nDone.\n")
