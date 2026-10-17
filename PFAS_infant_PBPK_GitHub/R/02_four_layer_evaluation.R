# ============================================================
# PFAS infant PBPK — FOUR-LAYER EVALUATION / QUALIFICATION MODULE — FROZEN PIPELINE (graphics warning fix 2026-10-16)
#
# Purpose
#   Upgrade the current two-layer evaluation (adult + 6-month infant)
#   to a four-layer evidence framework:
#
#   Layer 1. Adult renal qualification (Abraham 2024)
#   Layer 2. Neonatal renal-clearance external benchmark (Yao et al. 2023)
#   Layer 3. 6-month infant population-level external evaluation (Varsi 2022)
#   Layer 4. Chinese childhood longitudinal plausibility benchmark (Wu et al. 2026)
#
# IMPORTANT SCIENTIFIC INTERPRETATION
#   Layer 1 = model qualification, not infant validation.
#   Layer 2 = mechanistic external benchmark of net renal clearance.
#   Layer 3 = population-median external evaluation, not individual-level validation.
#   Layer 4 = longitudinal plausibility benchmark, not formal validation,
#             because the current PBPK model stops at 12 months whereas Wu et al.
#             followed children to 2, 4 and 7 years and includes post-weaning
#             environmental/dietary exposure.
#
# Expected PBPK result columns from the current unified Tm/Kt model:
#   Time_yr, Time_days, ID, Compound, Scenario,
#   BW_current_kg, C_Plasma,
#   R_Filtration, R_Reabsorption, R_Urine_Excretion
#
# Units expected:
#   C_Plasma             ng/L
#   R_*                  ng/h
#   BW_current_kg        kg
#
# Primary output:
#   four_layer_evaluation_summary.xlsx
#   plus publication-quality PNG/PDF figures
# ============================================================

options(stringsAsFactors = FALSE)

# ============================================================
# 0. PACKAGES
# ============================================================

required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "ggplot2",
  "writexl", "purrr", "stringr", "scales"
)

new_pkgs <- required_pkgs[
  !(required_pkgs %in% rownames(installed.packages()))
]

if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(writexl)
  library(purrr)
  library(stringr)
  library(scales)
})

# ============================================================
# 1. USER SETTINGS
# ============================================================

# Freeze-compatible paths. Run from repository root or set PFAS_PBPK_ROOT.
PROJECT_ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", unset = getwd()),
                              winslash = "/", mustWork = FALSE)
MODEL_FREEZE_VERSION <- "PFAS_infant_PBPK_final_freeze_v1_1_GitHub"
MODEL_RESULTS_DIR <- file.path(PROJECT_ROOT, "results", MODEL_FREEZE_VERSION)
PBPK_RESULTS_XLSX <- Sys.getenv("PFAS_PBPK_RESULTS", unset = file.path(
  MODEL_RESULTS_DIR, "tables", "PBPK_all_S0_S1_unified_TmKt_UrinaryStorage_results.xlsx"))
# Optional supporting dataset; not required for the main four-layer comparison.
WU_SUPP_XLSX <- file.path(PROJECT_ROOT, "data", "external",
                          "1-s2.0-S0160412026002813-mmc2.xlsx")
OUTPUT_DIR <- file.path(MODEL_RESULTS_DIR, "four_layer_evaluation")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
PLOT_DIR <- file.path(OUTPUT_DIR, "figures")
LOG_DIR <- file.path(OUTPUT_DIR, "provenance")
for (d in c(TABLE_DIR, PLOT_DIR, LOG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Main infant scenario used for mechanistic / longitudinal evaluation.
PRIMARY_SCENARIO <- "S1"

# Neonatal ages used for CLrenal comparison.
# Current PBPK output is usually every 7 days.
NEONATAL_COMPARE_DAYS <- c(7, 14, 28)

# For matching a requested age to the closest model output.
MAX_DAY_MATCH_DIFF <- 4

# Descriptive 2-fold agreement only; not a universal pass/fail criterion.
TWO_FOLD_LOW  <- 0.5
TWO_FOLD_HIGH <- 2.0

# ============================================================
# 2. HELPERS
# ============================================================

safe_read_pbpk <- function(path) {
  if (!file.exists(path)) {
    stop("PBPK result file not found:\n", path)
  }

  sheets <- readxl::excel_sheets(path)

  preferred <- c(
    "PBPK_Results",
    "Results",
    "results",
    "Sheet1"
  )

  sheet_use <- preferred[preferred %in% sheets][1]

  if (length(sheet_use) == 0 || is.na(sheet_use)) {
    sheet_use <- sheets[1]
  }

  message("Reading PBPK results from sheet: ", sheet_use)

  readxl::read_excel(path, sheet = sheet_use)
}


safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
}


safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x, na.rm = TRUE)
}


safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x, na.rm = TRUE)
}


normalize_compound <- function(x) {
  x <- as.character(x)
  x[x %in% c("PFUnDA", "PFUnA")] <- "PFUdA"
  x
}


calc_fold_metrics <- function(df,
                              pred_col = "Predicted",
                              obs_col  = "Observed") {

  pred <- as.numeric(df[[pred_col]])
  obs  <- as.numeric(df[[obs_col]])

  keep <- is.finite(pred) &
    is.finite(obs) &
    pred > 0 &
    obs > 0

  pred <- pred[keep]
  obs  <- obs[keep]

  if (!length(pred)) {
    return(tibble::tibble(
      N = 0L,
      MFE = NA_real_,
      AFE = NA_real_,
      RMSE = NA_real_,
      R2_linear = NA_real_,
      R2_log10 = NA_real_,
      TwoFold_Accuracy = NA_real_,
      Median_PO = NA_real_,
      Min_PO = NA_real_,
      Max_PO = NA_real_
    ))
  }

  ratio <- pred / obs

  r2_linear <- if (
    length(pred) >= 2 &&
    sd(pred) > 0 &&
    sd(obs) > 0
  ) {
    cor(pred, obs)^2
  } else {
    NA_real_
  }

  r2_log <- if (
    length(pred) >= 2 &&
    sd(log10(pred)) > 0 &&
    sd(log10(obs)) > 0
  ) {
    cor(log10(pred), log10(obs))^2
  } else {
    NA_real_
  }

  tibble::tibble(
    N = length(pred),
    MFE = exp(mean(log(ratio))),
    AFE = exp(mean(abs(log(ratio)))),
    RMSE = sqrt(mean((pred - obs)^2)),
    R2_linear = r2_linear,
    R2_log10 = r2_log,
    TwoFold_Accuracy = mean(
      ratio >= TWO_FOLD_LOW &
        ratio <= TWO_FOLD_HIGH
    ),
    Median_PO = median(ratio),
    Min_PO = min(ratio),
    Max_PO = max(ratio)
  )
}


nearest_age_rows <- function(df, target_day) {

  df %>%
    group_by(ID, Compound, Scenario) %>%
    mutate(
      Day_Diff = abs(Time_days - target_day)
    ) %>%
    filter(Day_Diff == min(Day_Diff, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup() %>%
    filter(Day_Diff <= MAX_DAY_MATCH_DIFF)
}


save_plot_both <- function(
    p, filename, width = 9, height = 6
) {
  ggsave(
    file.path(PLOT_DIR, paste0(filename, ".png")),
    p,
    width = width,
    height = height,
    dpi = 600
  )

  ggsave(
    file.path(PLOT_DIR, paste0(filename, ".pdf")),
    p,
    width = width,
    height = height
  )
}


# ============================================================
# 3. LOAD CURRENT FINAL PBPK RESULTS
# ============================================================

pbpk <- safe_read_pbpk(PBPK_RESULTS_XLSX)
message("Using frozen main-model output: ", normalizePath(PBPK_RESULTS_XLSX, winslash="/"))

required_cols <- c(
  "Time_yr", "Time_days", "ID",
  "Compound", "Scenario",
  "BW_current_kg", "C_Plasma",
  "R_Filtration", "R_Reabsorption",
  "R_Urine_Excretion"
)

missing_cols <- setdiff(required_cols, names(pbpk))

if (length(missing_cols) > 0) {
  stop(
    "PBPK result file is missing required columns:\n",
    paste(missing_cols, collapse = ", ")
  )
}

if (!PRIMARY_SCENARIO %in% as.character(pbpk$Scenario))
  stop("Required scenario ", PRIMARY_SCENARIO, " absent from main-model results")
if (anyDuplicated(names(pbpk))) stop("Duplicate input column names found")

pbpk <- pbpk %>%
  mutate(
    Compound = normalize_compound(Compound),
    Scenario = as.character(Scenario),
    Time_days = as.numeric(Time_days),
    Time_yr = as.numeric(Time_yr),
    BW_current_kg = as.numeric(BW_current_kg),
    C_Plasma = as.numeric(C_Plasma),
    R_Filtration = as.numeric(R_Filtration),
    R_Reabsorption = as.numeric(R_Reabsorption),
    R_Urine_Excretion = as.numeric(R_Urine_Excretion)
  )

cat("\n============================================================\n")
cat("FOUR-LAYER EVALUATION — GITHUB FREEZE\n")
cat("============================================================\n")
cat("Rows in PBPK results: ", nrow(pbpk), "\n", sep = "")
cat("Compounds: ",
    paste(sort(unique(pbpk$Compound)), collapse = ", "),
    "\n",
    sep = "")
cat("Scenarios: ",
    paste(sort(unique(pbpk$Scenario)), collapse = ", "),
    "\n\n",
    sep = "")


# ============================================================
# 4. LAYER 1 — ADULT RENAL QUALIFICATION
# ============================================================
#
# These are the final values reported in the current manuscript.
#
# Interpretation:
#   PFOA/PFNA/PFHxS/PFOS:
#     parameterizations were retained without fitting to Abraham data.
#
#   PFDA/PFUdA:
#     Tmc was calibrated to the Abraham adult terminal kinetics.
#     Therefore model/observed = 1.00 is a calibration result, not
#     independent predictive performance.
# ============================================================

adult_layer <- tibble::tribble(
  ~Compound, ~Observed_HalfLife_d, ~Predicted_HalfLife_d,
  ~Parameter_Use, ~Evaluation_Class,

  "PFOA", 2348.1, 1392.9,
  "Fixed before Abraham benchmark",
  "External adult qualification",

  "PFNA", 1581.3, 1670.6,
  "Fixed before Abraham benchmark",
  "External adult qualification",

  "PFDA", 1122.6, 1122.6,
  "Tmc calibrated to Abraham terminal kinetics",
  "Adult calibration only",

  "PFUdA", 568.2, 568.2,
  "Tmc calibrated to Abraham terminal kinetics",
  "Adult calibration only",

  "PFHxS", 1885.2, 3756.5,
  "Fixed before Abraham benchmark",
  "External adult qualification",

  "PFOS", 1420.8, 2412.8,
  "Fixed before Abraham benchmark",
  "External adult qualification"
) %>%
  mutate(
    PO_Ratio = Predicted_HalfLife_d / Observed_HalfLife_d,
    Within_2Fold =
      PO_Ratio >= TWO_FOLD_LOW &
      PO_Ratio <= TWO_FOLD_HIGH,
    Layer = "Layer 1",
    Layer_Name = "Adult renal qualification"
  )

adult_fixed_only <- adult_layer %>%
  filter(Evaluation_Class == "External adult qualification")

adult_metrics <- calc_fold_metrics(
  adult_fixed_only %>%
    transmute(
      Predicted = Predicted_HalfLife_d,
      Observed = Observed_HalfLife_d
    )
)


# ============================================================
# 5. LAYER 2 — NEONATAL RENAL-CLEARANCE EXTERNAL BENCHMARK
# ============================================================
#
# Yao et al. 2023 neonatal observations:
#   PFOA  = 0.049 mL/kg BW/day
#   PFNA  = 0.047
#   PFDA  = 0.022
#   PFHxS = 0.270
#
# PFOS:
#   PFOS was not detected in neonatal urine in that study;
#   therefore no infant CLrenal benchmark is entered here.
#
# Model-derived renal clearance:
#
#   CLrenal_excretion
#      = R_Urine_Excretion / C_Plasma
#
#   R_Urine_Excretion: ng/h
#   C_Plasma: ng/L
#   => L/h
#   => ×1000 mL/L ×24 h/day / BW(kg)
#
# A second diagnostic is also reported:
#
#   CLrenal_net_entry
#      = max(R_Filtration - R_Reabsorption, 0) / C_Plasma
#
# This is NOT the primary comparison with Yao, but it helps diagnose
# urinary-storage transient effects in the current model.
# ============================================================

yao_obs <- tibble::tribble(
  ~Compound, ~Observed_CLrenal_mL_kg_day,
  ~Source_Note,

  "PFOA",  0.049,
  "Yao et al. 2023; paired neonatal cord blood and urine",

  "PFNA",  0.047,
  "Yao et al. 2023; paired neonatal cord blood and urine",

  "PFDA",  0.022,
  "Yao et al. 2023; paired neonatal cord blood and urine",

  "PFHxS", 0.270,
  "Yao et al. 2023; paired neonatal cord blood and urine"
)


pbpk_cl <- pbpk %>%
  filter(Scenario == PRIMARY_SCENARIO) %>%
  mutate(
    CLrenal_Excretion_mL_kg_day = ifelse(
      is.finite(R_Urine_Excretion) &
        is.finite(C_Plasma) &
        C_Plasma > 0 &
        is.finite(BW_current_kg) &
        BW_current_kg > 0,

      (R_Urine_Excretion / C_Plasma) *
        1000 * 24 / BW_current_kg,

      NA_real_
    ),

    Net_Renal_Entry_ng_h = pmax(
      R_Filtration - R_Reabsorption,
      0
    ),

    CLrenal_NetEntry_mL_kg_day = ifelse(
      is.finite(Net_Renal_Entry_ng_h) &
        is.finite(C_Plasma) &
        C_Plasma > 0 &
        is.finite(BW_current_kg) &
        BW_current_kg > 0,

      (Net_Renal_Entry_ng_h / C_Plasma) *
        1000 * 24 / BW_current_kg,

      NA_real_
    )
  )


layer2_by_age <- purrr::map_dfr(
  NEONATAL_COMPARE_DAYS,
  function(day_target) {

    nearest_age_rows(pbpk_cl, day_target) %>%
      filter(Compound %in% yao_obs$Compound) %>%
      group_by(Compound) %>%
      summarise(
        Target_Day = day_target,
        Actual_Model_Day =
          safe_median(Time_days),

        N_Individuals =
          sum(is.finite(CLrenal_Excretion_mL_kg_day)),

        Pred_CLrenal_Excretion_Median =
          safe_median(CLrenal_Excretion_mL_kg_day),

        Pred_CLrenal_Excretion_P025 =
          safe_quantile(CLrenal_Excretion_mL_kg_day, 0.025),

        Pred_CLrenal_Excretion_P975 =
          safe_quantile(CLrenal_Excretion_mL_kg_day, 0.975),

        Pred_CLrenal_NetEntry_Median =
          safe_median(CLrenal_NetEntry_mL_kg_day),

        Pred_CLrenal_NetEntry_P025 =
          safe_quantile(CLrenal_NetEntry_mL_kg_day, 0.025),

        Pred_CLrenal_NetEntry_P975 =
          safe_quantile(CLrenal_NetEntry_mL_kg_day, 0.975),

        .groups = "drop"
      ) %>%
      left_join(
        yao_obs,
        by = "Compound"
      ) %>%
      mutate(
        PO_Ratio_Excretion =
          Pred_CLrenal_Excretion_Median /
          Observed_CLrenal_mL_kg_day,

        PO_Ratio_NetEntry =
          Pred_CLrenal_NetEntry_Median /
          Observed_CLrenal_mL_kg_day,

        Within_2Fold_Excretion =
          PO_Ratio_Excretion >= TWO_FOLD_LOW &
          PO_Ratio_Excretion <= TWO_FOLD_HIGH,
        Clearance_Status = ifelse(
          is.finite(Pred_CLrenal_Excretion_Median),
          "Estimable", "Not estimable: urinary elimination parameter/output missing"
        ),

        Layer = "Layer 2",
        Layer_Name =
          "Neonatal renal-clearance external benchmark"
      )
  }
)


# Day 7 is primary because the source study used neonatal urine collected
# during the first postnatal week and the current model reports every 7 days.
layer2_primary <- layer2_by_age %>%
  filter(Target_Day == 7)


layer2_metrics <- calc_fold_metrics(
  layer2_primary %>%
    transmute(
      Predicted = Pred_CLrenal_Excretion_Median,
      Observed = Observed_CLrenal_mL_kg_day
    )
)


# ============================================================
# 6. LAYER 3 — VARSI 2022 6-MONTH INFANT EXTERNAL EVALUATION
# ============================================================
#
# Final frozen study-specific predictions and observed population medians
# reported in the current manuscript.
#
# This is a POPULATION-MEDIAN external evaluation.
# It must not be described as individual-level prediction validation.
# ============================================================

varsi_layer <- tibble::tribble(
  ~Compound, ~Predicted_ng_mL, ~Observed_ng_mL,

  "PFOA",  2.853, 3.270,
  "PFNA",  0.834, 0.710,
  "PFDA",  0.200, 0.230,
  "PFUdA", 0.161, 0.160,
  "PFHxS", 0.434, 0.680,
  "PFOS",  6.745, 4.710
) %>%
  mutate(
    Prediction_Provenance = "Previously frozen Varsi-specific population predictions; not recomputed from the main cohort",
    PO_Ratio = Predicted_ng_mL / Observed_ng_mL,
    Within_2Fold =
      PO_Ratio >= TWO_FOLD_LOW &
      PO_Ratio <= TWO_FOLD_HIGH,
    Layer = "Layer 3",
    Layer_Name =
      "6-month infant population-level external evaluation"
  )


varsi_metrics <- calc_fold_metrics(
  varsi_layer %>%
    transmute(
      Predicted = Predicted_ng_mL,
      Observed = Observed_ng_mL
    )
)


# ============================================================
# 7. LAYER 4 — WU 2026 LONGITUDINAL PLAUSIBILITY BENCHMARK
# ============================================================
#
# Wu et al. 2026 Shanghai Birth Cohort:
# adjusted annual mean percent changes from birth to 7 years:
#
#   PFOA   +9.6%
#   PFNA  +11.4%
#   PFDA  +10.8%
#   PFUdA  +6.7%
#   PFHxS +31.1%
#   PFOS   +9.6%
#
# Current model only covers birth to 12 months and does not include later
# dietary/drinking-water/dust exposure. Therefore, the correct comparison
# is direction / qualitative accumulation plausibility, not absolute
# concentration validation.
# ============================================================

wu_trend <- tibble::tribble(
  ~Compound, ~Wu_Adjusted_Annual_Change_pct,

  "PFOA",   9.6,
  "PFNA",  11.4,
  "PFDA",  10.8,
  "PFUdA",  6.7,
  "PFHxS", 31.1,
  "PFOS",   9.6
)


# Subject-specific birth and ~12-month plasma concentration.
model_birth_12 <- pbpk %>%
  filter(
    Scenario == PRIMARY_SCENARIO,
    Compound %in% wu_trend$Compound
  ) %>%
  group_by(ID, Compound, Scenario) %>%
  arrange(Time_days, .by_group = TRUE) %>%
  summarise(
    Birth_Day =
      Time_days[which.min(abs(Time_days - 0))][1],

    C_Birth_ng_L =
      C_Plasma[which.min(abs(Time_days - 0))][1],

    Month12_Day =
      Time_days[which.min(abs(Time_days - 365.25))][1],

    C_12mo_ng_L =
      C_Plasma[which.min(abs(Time_days - 365.25))][1],

    .groups = "drop"
  ) %>%
  mutate(
    Model_12mo_vs_Birth_pct = ifelse(
      is.finite(C_Birth_ng_L) &
        C_Birth_ng_L > 0 &
        is.finite(C_12mo_ng_L),

      100 * (C_12mo_ng_L / C_Birth_ng_L - 1),

      NA_real_
    ),

    Model_Log_Annual_Change_pct = ifelse(
      is.finite(C_Birth_ng_L) &
        C_Birth_ng_L > 0 &
        is.finite(C_12mo_ng_L) &
        C_12mo_ng_L > 0,

      100 * (exp(log(C_12mo_ng_L / C_Birth_ng_L)) - 1),

      NA_real_
    )
  )


layer4 <- model_birth_12 %>%
  group_by(Compound) %>%
  summarise(
    N = sum(is.finite(Model_12mo_vs_Birth_pct)),

    Model_12mo_Change_Median_pct =
      safe_median(Model_12mo_vs_Birth_pct),

    Model_12mo_Change_P025_pct =
      safe_quantile(Model_12mo_vs_Birth_pct, 0.025),

    Model_12mo_Change_P975_pct =
      safe_quantile(Model_12mo_vs_Birth_pct, 0.975),

    Median_Birth_C_ng_L =
      safe_median(C_Birth_ng_L),

    Median_12mo_C_ng_L =
      safe_median(C_12mo_ng_L),

    .groups = "drop"
  ) %>%
  left_join(
    wu_trend,
    by = "Compound"
  ) %>%
  mutate(
    Model_Direction =
      case_when(
        Model_12mo_Change_Median_pct > 0 ~ "Increase",
        Model_12mo_Change_Median_pct < 0 ~ "Decrease",
        TRUE ~ "No change"
      ),

    Wu_Direction =
      case_when(
        Wu_Adjusted_Annual_Change_pct > 0 ~ "Increase",
        Wu_Adjusted_Annual_Change_pct < 0 ~ "Decrease",
        TRUE ~ "No change"
      ),

    Direction_Concordant =
      Model_Direction == Wu_Direction,

    Layer = "Layer 4",
    Layer_Name =
      "Chinese childhood longitudinal plausibility benchmark",

    Interpretation = paste0(
      "Direction-only plausibility comparison. ",
      "Wu follows birth to 7 years; current PBPK stops at 12 months ",
      "and omits later diet/water/dust exposure."
    )
  )


# ============================================================
# 8. FOUR-LAYER EVIDENCE SUMMARY
# ============================================================

layer_summary <- tibble::tribble(
  ~Layer, ~Evidence, ~Target, ~Primary_Output,
  ~Interpretation_Level,

  "Layer 1",
  "Abraham et al. 2024 adult long-term kinetics",
  "Shared renal elimination structure / terminal kinetics",
  "Observed vs predicted terminal half-life",
  "Adult model qualification",

  "Layer 2",
  "Yao et al. 2023 paired neonatal cord blood + urine",
  "Net renal clearance in early infancy",
  "Observed vs PBPK-derived CLrenal",
  "Mechanistic external benchmark",

  "Layer 3",
  "Varsi et al. 2022 6-month infant serum",
  "Infant systemic concentration",
  "Predicted vs observed population medians",
  "Population-level external evaluation",

  "Layer 4",
  "Wu et al. 2026 Shanghai Birth Cohort",
  "Longitudinal direction of child PFAS burden",
  "Birth-to-12mo model direction vs birth-to-7y observed trend",
  "Longitudinal plausibility benchmark"
)


metric_summary <- bind_rows(
  adult_metrics %>%
    mutate(
      Layer = "Layer 1",
      Dataset = "Abraham 2024 — fixed compounds only",
      Metric_Unit = "days"
    ),

  layer2_metrics %>%
    mutate(
      Layer = "Layer 2",
      Dataset = "Yao 2023 — day 7 primary CLrenal comparison",
      Metric_Unit = "mL/kg/day"
    ),

  varsi_metrics %>%
    mutate(
      Layer = "Layer 3",
      Dataset = "Varsi 2022 — six compound population medians",
      Metric_Unit = "ng/mL"
    )
) %>%
  select(
    Layer, Dataset, Metric_Unit,
    everything()
  )


# ============================================================
# 9. FIGURE 1 — FOUR-LAYER VALIDATION WORKFLOW
# ============================================================

workflow_plot_df <- tibble::tibble(
  x = 1:4,
  y = 1,
  Layer = paste0("Layer ", 1:4),
  Label = c(
    "Adult renal\nqualification",
    "Neonatal renal\nclearance",
    "6-month infant\nexternal evaluation",
    "Childhood longitudinal\nplausibility"
  )
)

p_workflow <- ggplot(
  workflow_plot_df,
  aes(x = x, y = y)
) +
  annotate(
    "segment", x = 1, xend = 4, y = 1, yend = 1,
    linewidth = 0.7
  ) +
  geom_point(size = 7) +
  geom_text(
    aes(label = Layer),
    vjust = -1.2,
    fontface = "bold",
    size = 4.3
  ) +
  geom_text(
    aes(label = Label),
    vjust = 1.9,
    size = 3.8,
    lineheight = 0.95
  ) +
  scale_x_continuous(
    limits = c(0.65, 4.35),
    breaks = NULL
  ) +
  scale_y_continuous(
    limits = c(0.78, 1.22),
    breaks = NULL
  ) +
  labs(
    title = "Four-layer qualification and external-evaluation framework",
    subtitle =
      "Adult kinetics -> neonatal renal clearance -> 6-month infant concentration -> longitudinal childhood plausibility"
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5
    )
  )

save_plot_both(
  p_workflow,
  "Figure_1_four_layer_validation_workflow",
  width = 12,
  height = 4.3
)


# ============================================================
# 10. FIGURE 2 — ADULT HALF-LIFE QUALIFICATION
# ============================================================

p_adult <- ggplot(
  adult_layer,
  aes(
    x = Observed_HalfLife_d,
    y = Predicted_HalfLife_d,
    shape = Evaluation_Class,
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
    slope = 0.5,
    linetype = 2,
    linewidth = 0.6
  ) +
  geom_abline(
    intercept = 0,
    slope = 2,
    linetype = 2,
    linewidth = 0.6
  ) +
  geom_point(size = 3.5) +
  geom_text(
    nudge_x = 60,
    nudge_y = 70,
    check_overlap = TRUE,
    size = 3.8
  ) +
  coord_equal() +
  labs(
    title = "Layer 1 — Adult renal qualification",
    subtitle =
      "PFDA/PFUdA points are adult-calibrated and should not be interpreted as independent predictions",
    x = "Observed terminal half-life (days)",
    y = "PBPK terminal half-life (days)",
    shape = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

save_plot_both(
  p_adult,
  "Figure_2_layer1_adult_half_life",
  width = 8.5,
  height = 7
)


# ============================================================
# 11. FIGURE 3 — NEONATAL CLrenal
# ============================================================

# Plot only estimable predictions. Keep ALL compounds/ages (including NA)
# in the exported workbook; no missing estimate is interpreted as zero.
layer2_plot_data <- layer2_by_age %>%
  dplyr::filter(
    is.finite(Pred_CLrenal_Excretion_Median),
    is.finite(Pred_CLrenal_Excretion_P025),
    is.finite(Pred_CLrenal_Excretion_P975)
  )
renal_missing_n <- nrow(layer2_by_age) - nrow(layer2_plot_data)
if (renal_missing_n > 0L) {
  message("Layer 2 plot: ", renal_missing_n,
          " non-estimable point(s) omitted from figure only; retained in workbook.")
}

p_renal <- layer2_plot_data %>%
  ggplot(
    aes(
      x = factor(Target_Day),
      y = Pred_CLrenal_Excretion_Median,
      group = Compound
    )
  ) +
  geom_hline(
    data = yao_obs,
    aes(
      yintercept = Observed_CLrenal_mL_kg_day
    ),
    linetype = 2,
    linewidth = 0.65
  ) +
  geom_errorbar(
    aes(
      ymin = Pred_CLrenal_Excretion_P025,
      ymax = Pred_CLrenal_Excretion_P975
    ),
    width = 0.12
  ) +
  geom_point(size = 2.8) +
  geom_line(linewidth = 0.7) +
  facet_wrap(
    ~Compound,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = "Layer 2 — Neonatal renal-clearance external benchmark",
    subtitle =
      "Points: PBPK population median; bars: 2.5th–97.5th percentiles; dashed line: Yao et al. neonatal estimate",
    x = "Postnatal age used for PBPK comparison (days)",
    y = "Renal clearance (mL/kg/day)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

save_plot_both(
  p_renal,
  "Figure_3_layer2_neonatal_CLrenal",
  width = 9.5,
  height = 7.5
)


# ============================================================
# 12. FIGURE 4 — VARSI 6-MONTH EXTERNAL EVALUATION
# ============================================================

p_varsi <- ggplot(
  varsi_layer,
  aes(
    x = Observed_ng_mL,
    y = Predicted_ng_mL,
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
    slope = 0.5,
    linetype = 2,
    linewidth = 0.6
  ) +
  geom_abline(
    intercept = 0,
    slope = 2,
    linetype = 2,
    linewidth = 0.6
  ) +
  geom_point(size = 3.2) +
  geom_text(
    nudge_x = 0.07,
    nudge_y = 0.07,
    size = 3.8,
    check_overlap = TRUE
  ) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal() +
  labs(
    title = "Layer 3 — External evaluation at 6 months",
    subtitle =
      "Study-specific frozen PBPK predictions vs Varsi 2022 population medians",
    x = "Observed infant serum median (ng/mL)",
    y = "Predicted infant plasma median (ng/mL)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

save_plot_both(
  p_varsi,
  "Figure_4_layer3_Varsi_predicted_vs_observed",
  width = 7.5,
  height = 7
)


# ============================================================
# 13. FIGURE 5 — WU LONGITUDINAL PLAUSIBILITY
# ============================================================

layer4_long <- layer4 %>%
  transmute(
    Compound,
    `PBPK birth to 12 mo` =
      Model_12mo_Change_Median_pct,
    `Wu birth to 7 y (annual trend)` =
      Wu_Adjusted_Annual_Change_pct
  ) %>%
  pivot_longer(
    cols = -Compound,
    names_to = "Evidence",
    values_to = "Percent_Change"
  )

layer4_plot_data <- layer4_long %>%
  dplyr::filter(is.finite(Percent_Change))
wu_missing_n <- nrow(layer4_long) - nrow(layer4_plot_data)
if (wu_missing_n > 0L) {
  message("Layer 4 plot: ", wu_missing_n,
          " non-estimable point(s) omitted from figure only; retained in workbook.")
}
p_wu <- ggplot(
  layer4_plot_data,
  aes(
    x = Compound,
    y = Percent_Change,
    shape = Evidence,
    group = Evidence
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.6
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title = "Layer 4 — Longitudinal plausibility benchmark",
    subtitle =
      "Magnitude is not expected to match because Wu spans 0–7 years and includes post-weaning exposure; direction is the primary comparison",
    x = NULL,
    y = "Percent change",
    shape = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

save_plot_both(
  p_wu,
  "Figure_5_layer4_Wu_longitudinal_plausibility",
  width = 9,
  height = 6.2
)


# ============================================================
# 14. OPTIONAL: INSPECT WU SUPPLEMENT SHEETS
# ============================================================
#
# This section only creates an index of sheet names so the next version can
# add a direct age-2 concentration comparison without manually searching
# the supplementary workbook.
# ============================================================

wu_supp_sheet_index <- tibble::tibble()

if (file.exists(WU_SUPP_XLSX)) {
  wu_sheets <- readxl::excel_sheets(WU_SUPP_XLSX)

  wu_supp_sheet_index <- tibble::tibble(
    Sheet = wu_sheets,
    Candidate_S5 =
      str_detect(
        tolower(wu_sheets),
        "s5|table.?5"
      )
  )
}


# ============================================================
# 15. SCIENTIFIC INTERPRETATION TABLE
# ============================================================

interpretation_notes <- tibble::tribble(
  ~Item, ~Interpretation,

  "Layer 1",
  paste0(
    "Adult Abraham data qualify the shared renal architecture. ",
    "PFDA and PFUdA are calibration compounds in this layer."
  ),

  "Layer 2 primary age",
  paste0(
    "Day 7 is primary because Yao neonatal urine was collected ",
    "during the first postnatal week and the current PBPK report grid is weekly."
  ),

  "Layer 2 CLrenal definition",
  paste0(
    "Primary PBPK CLrenal = urinary excretion rate / plasma concentration, ",
    "converted to mL/kg/day. A net-filtration-minus-reabsorption diagnostic ",
    "is reported separately."
  ),

  "Layer 2 PFDA/PFUdA",
  "Unidentified KurineC makes urinary clearance NA; do not replace missing predictions with zero or compute fold errors.",

  "Layer 3 provenance",
  "Frozen Varsi-specific model medians reproduced from the preceding evaluation, not recalculated from the present infant cohort. Re-run external inputs to refresh this comparison.",

  "Layer 2 limitation",
  paste0(
    "Yao clearance estimates are based on paired neonatal blood/urine samples ",
    "and should be treated as an external mechanistic benchmark, not a parameter-fitting target."
  ),

  "Layer 3",
  paste0(
    "Varsi comparison is based on six compound-level population medians. ",
    "Metrics such as R2 therefore do not represent individual-level predictive performance."
  ),

  "Layer 4",
  paste0(
    "Wu 2026 is a plausibility benchmark only in the current version. ",
    "The PBPK model stops at 12 months, whereas Wu follows children to 7 years ",
    "and later childhood includes diet, water, dust and other exposure sources."
  ),

  "2-fold criterion",
  paste0(
    "0.5–2.0 predicted/observed is reported descriptively and is not used ",
    "as a single model-acceptance rule."
  ),

  "PFOS in Yao",
  paste0(
    "PFOS was not detected in neonatal urine in the Yao study; ",
    "therefore no PFOS infant CLrenal observation is used in Layer 2."
  )
)


# ============================================================
# 16. WRITE OUTPUT WORKBOOK
# ============================================================

out_xlsx <- file.path(
  TABLE_DIR,
  "four_layer_evaluation_summary.xlsx"
)

safe_write <- function(x, path) {
  ok <- tryCatch(
    {
      writexl::write_xlsx(x, path)
      TRUE
    },
    error = function(e) {
      message("write_xlsx failed: ", conditionMessage(e))
      FALSE
    }
  )

  if (ok) {
    return(path)
  }

  alt <- file.path(
    dirname(path),
    paste0(
      tools::file_path_sans_ext(
        basename(path)
      ),
      "_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".xlsx"
    )
  )

  writexl::write_xlsx(x, alt)

  warning(
    "Could not overwrite requested workbook. Saved instead to:\n",
    alt
  )

  alt
}


written_file <- safe_write(
  list(
    Framework = layer_summary,
    Interpretation_Notes = interpretation_notes,

    L1_Adult_Qualification = adult_layer,
    L1_Fixed_Only_Metrics = adult_metrics,

    L2_Yao_Observed = yao_obs,
    L2_CLrenal_All_Ages = layer2_by_age,
    L2_CLrenal_Primary_Day7 = layer2_primary,
    L2_Day7_Metrics = layer2_metrics,

    L3_Varsi_Comparison = varsi_layer,
    L3_Varsi_Metrics = varsi_metrics,

    L4_Wu_Trend = wu_trend,
    L4_PBPK_Subject_Changes = model_birth_12,
    L4_Wu_Plausibility = layer4,

    Overall_Metric_Summary = metric_summary,

    Wu_Supplement_Sheet_Index =
      wu_supp_sheet_index
  ),
  out_xlsx
)


# ============================================================
# 17. CONSOLE REPORT
# ============================================================

cat("\n============================================================\n")
cat("FOUR-LAYER EVALUATION COMPLETED\n")
cat("============================================================\n")

cat("\nLayer 1 — Adult fixed-compound metrics:\n")
print(adult_metrics)

cat("\nLayer 2 — Yao neonatal CLrenal, primary Day 7:\n")
print(
  layer2_primary %>%
    select(
      Compound,
      Observed_CLrenal_mL_kg_day,
      Pred_CLrenal_Excretion_Median,
      Pred_CLrenal_Excretion_P025,
      Pred_CLrenal_Excretion_P975,
      PO_Ratio_Excretion,
      Within_2Fold_Excretion
    ),
  n = Inf
)

cat("\nLayer 2 — Day-7 overall descriptive metrics:\n")
print(layer2_metrics)

cat("\nLayer 3 — Varsi six-compound external evaluation:\n")
print(varsi_layer, n = Inf)

cat("\nLayer 3 — Overall metrics:\n")
print(varsi_metrics)

cat("\nLayer 4 — Wu longitudinal plausibility:\n")
print(
  layer4 %>%
    select(
      Compound,
      Model_12mo_Change_Median_pct,
      Model_12mo_Change_P025_pct,
      Model_12mo_Change_P975_pct,
      Wu_Adjusted_Annual_Change_pct,
      Model_Direction,
      Wu_Direction,
      Direction_Concordant
    ),
  n = Inf
)

cat("\nOutput workbook:\n", written_file, "\n", sep = "")
cat("Plot directory:\n", PLOT_DIR, "\n", sep = "")

cat("\nScientific labels to use in manuscript:\n")
cat("  Layer 1: adult renal qualification\n")
cat("  Layer 2: neonatal renal-clearance external benchmark\n")
cat("  Layer 3: population-level infant external evaluation\n")
cat("  Layer 4: longitudinal childhood plausibility benchmark\n")

cat("\nDone.\n")

# Reproducibility metadata (does not alter results).
provenance <- data.frame(
  Key = c("Model_freeze_version", "Evaluation_script", "Main_result_path",
          "Main_result_md5", "Primary_scenario", "Run_timestamp", "R_version"),
  Value = c(MODEL_FREEZE_VERSION, "R/02_four_layer_evaluation.R",
            normalizePath(PBPK_RESULTS_XLSX, winslash="/"),
            unname(tools::md5sum(PBPK_RESULTS_XLSX)), PRIMARY_SCENARIO,
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), R.version.string)
)
writexl::write_xlsx(list(Provenance = provenance),
                    file.path(LOG_DIR, "evaluation_provenance.xlsx"))
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo.txt"))
