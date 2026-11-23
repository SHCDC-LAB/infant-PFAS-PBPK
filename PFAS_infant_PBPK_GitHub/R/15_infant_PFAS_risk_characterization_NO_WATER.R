#!/usr/bin/env Rscript
# =============================================================================
# R/15_infant_PFAS_risk_characterization_NO_WATER.R
# Manuscript-aligned exposure and internal benchmark characterization.
# Uses the frozen S1 PBPK output; does NOT rerun ODEs.
#
# EDI is reconstructed from the model's recorded ingestion rate:
#   EDI = Ingestion_Rate_ng_h * 24 / BW_current_kg
# This retains the modelled milk concentration, age-dependent milk intake,
# breastfeeding mixture ratio, and body weight used in the PBPK simulation.
#
# Six PFAS: EDI and organ/plasma concentration summaries.
# Four PFAS (PFOA, PFNA, PFHxS, PFOS): EFSA serum-equivalent BMDL10 comparison.
# EFSA BMDL10 = 17.5 ng/mL; this is a point of departure, not a clinical cutoff.
# No infant EDI-vs-TWI comparison; no water analyses; no organ HQ/HI.
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)
pkgs <- c("readxl", "dplyr", "tidyr", "tibble", "openxlsx")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
})

ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", unset = getwd()),
                      winslash = "/", mustWork = FALSE)
VERSION <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
MODE <- tolower(Sys.getenv("PFAS_TRAJECTORY_MODE", "proxy"))
if (MODE != "proxy") VERSION <- paste0(VERSION, "_beta_", MODE)

input <- Sys.getenv(
  "PFAS_PBPK_RESULTS",
  unset = file.path(ROOT, "results", VERSION, "tables",
                    "PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx")
)
if (!file.exists(input)) stop("Frozen PBPK workbook not found: ", input)
out_dir <- file.path(ROOT, "results", VERSION, "risk_characterization_NO_WATER")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output <- file.path(out_dir, "infant_PFAS_risk_characterization_NO_WATER.xlsx")

dat0 <- readxl::read_excel(input, sheet = "PBPK_Results")
required <- c("ID", "Compound", "Scenario", "Time_days", "BW_current_kg",
              "Milk_C", "Ingestion_Rate_ng_h", "C_Plasma", "C_Liver",
              "C_Kidney", "Mass_Balance_Relative_Error")
absent <- setdiff(required, names(dat0))
if (length(absent)) stop("PBPK_Results missing: ", paste(absent, collapse = ", "))

chem <- c("PFOA", "PFNA", "PFDA", "PFUdA", "PFHxS", "PFOS")
four <- c("PFOA", "PFNA", "PFHxS", "PFOS")
dat <- dat0 %>%
  mutate(ID = as.character(ID), Compound = as.character(Compound),
         Scenario = as.character(Scenario), Time_days = as.numeric(Time_days),
         BW_current_kg = as.numeric(BW_current_kg),
         Ingestion_Rate_ng_h = as.numeric(Ingestion_Rate_ng_h)) %>%
  filter(Scenario == "S1", Compound %in% chem) %>%
  arrange(ID, Compound, Time_days)

if (n_distinct(dat$ID) != 175L) stop("Expected 175 infants; found ", n_distinct(dat$ID))
if (!setequal(unique(dat$Compound), chem)) stop("The S1 output does not contain all six target PFAS.")
if (anyDuplicated(dat[c("ID", "Compound", "Time_days")])) stop("Duplicate ID-compound-time records.")
if (any(!is.finite(dat$BW_current_kg) | dat$BW_current_kg <= 0)) stop("Invalid body weight.")
if (any(!is.finite(dat$Ingestion_Rate_ng_h))) stop("Invalid ingestion rate.")

dat <- dat %>%
  mutate(EDI_ng_kg_day = Ingestion_Rate_ng_h * 24 / BW_current_kg)

# Trapezoidal time average with linear interpolation at exact period boundaries.
time_average <- function(t, y, start, end) {
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]; y <- y[ok]
  ord <- order(t); t <- t[ord]; y <- y[ord]
  if (length(t) < 2 || min(t) > start || max(t) < end)
    stop("Trajectory does not span requested interval ", start, "-", end)
  inside <- t > start & t < end
  tt <- c(start, t[inside], end)
  yy <- approx(t, y, xout = tt, method = "linear", ties = "ordered")$y
  sum(diff(tt) * (head(yy, -1) + tail(yy, -1)) / 2) / (end - start)
}

periods <- tibble::tribble(
  ~Period, ~Start_day, ~End_day,
  "0-180 days", 0, 180,
  "180-365 days", 180, 365,
  "0-365 days", 0, 365
)
edi_individual <- dat %>%
  group_by(ID, Compound) %>%
  group_modify(~ {
    bind_rows(lapply(seq_len(nrow(periods)), function(i) {
      tibble(Period = periods$Period[i],
             Mean_EDI_ng_kg_day = time_average(.x$Time_days, .x$EDI_ng_kg_day,
                                                periods$Start_day[i], periods$End_day[i]))
    }))
  }) %>% ungroup()

edi_summary <- edi_individual %>%
  group_by(Compound, Period) %>%
  summarise(N = n_distinct(ID),
            Median = median(Mean_EDI_ng_kg_day),
            Q25 = quantile(Mean_EDI_ng_kg_day, .25),
            Q75 = quantile(Mean_EDI_ng_kg_day, .75),
            P95 = quantile(Mean_EDI_ng_kg_day, .95), .groups = "drop")

mixture <- edi_individual %>%
  filter(Compound %in% chem) %>%
  group_by(ID, Period) %>%
  summarise(N_PFAS = n_distinct(Compound),
            Sum6_EDI = sum(Mean_EDI_ng_kg_day),
            Sum4_EDI = sum(Mean_EDI_ng_kg_day[Compound %in% four]),
            .groups = "drop")
if (any(mixture$N_PFAS != 6L)) stop("Incomplete six-compound EDI mixture.")

ratio <- tibble::tribble(
  ~Scenario, ~Compound, ~RSP,
  "Lower", "PFOA", 1.1, "Lower", "PFNA", 1.2,
  "Lower", "PFHxS", 1.0, "Lower", "PFOS", 1.1,
  "Central", "PFOA", 1.2, "Central", "PFNA", 1.3,
  "Central", "PFHxS", 1.1, "Central", "PFOS", 1.1,
  "Upper", "PFOA", 1.2, "Upper", "PFNA", 1.4,
  "Upper", "PFHxS", 1.1, "Upper", "PFOS", 1.2
)
day365 <- dat %>%
  filter(Compound %in% four) %>%
  group_by(ID, Compound) %>%
  group_modify(~ tibble(C_plasma_day365_ng_mL =
    approx(.x$Time_days, .x$C_Plasma, xout = 365, method = "linear")$y / 1000)) %>%
  ungroup()

ibq_component <- day365 %>%
  inner_join(ratio, by = "Compound") %>%
  mutate(Serum_equivalent_ng_mL = C_plasma_day365_ng_mL * RSP)
ibq_individual <- ibq_component %>%
  group_by(ID, Scenario) %>%
  summarise(Serum_equivalent_sum4_ng_mL = sum(Serum_equivalent_ng_mL),
            RBMDL10 = Serum_equivalent_sum4_ng_mL / 17.5,
            .groups = "drop")
ibq_summary <- ibq_individual %>%
  group_by(Scenario) %>%
  summarise(N = n_distinct(ID),
            Serum_sum_median = median(Serum_equivalent_sum4_ng_mL),
            Serum_sum_Q25 = quantile(Serum_equivalent_sum4_ng_mL, .25),
            Serum_sum_Q75 = quantile(Serum_equivalent_sum4_ng_mL, .75),
            RBMDL10_median = median(RBMDL10),
            RBMDL10_Q25 = quantile(RBMDL10, .25),
            RBMDL10_Q75 = quantile(RBMDL10, .75),
            N_RBMDL10_ge_1 = sum(RBMDL10 >= 1),
            Pct_RBMDL10_ge_1 = mean(RBMDL10 >= 1) * 100,
            .groups = "drop")

central_shares <- ibq_component %>%
  filter(Scenario == "Central") %>%
  group_by(ID) %>%
  mutate(Share_pct = 100 * Serum_equivalent_ng_mL / sum(Serum_equivalent_ng_mL)) %>%
  ungroup() %>%
  group_by(Compound) %>%
  summarise(N = n_distinct(ID), Median_share_pct = median(Share_pct),
            Q25 = quantile(Share_pct, .25), Q75 = quantile(Share_pct, .75),
            .groups = "drop")

mass_balance_qc <- dat %>%
  group_by(Compound) %>%
  summarise(N_ID = n_distinct(ID), N_rows = n(),
            Max_abs_relative_mass_balance_error =
              max(abs(Mass_Balance_Relative_Error), na.rm = TRUE),
            .groups = "drop")

methods_notes <- tibble(
  Item = c("Milk-derived EDI", "EFSA TWI", "Internal benchmark",
           "Serum/plasma ratios", "Organ-specific HQ/HI", "Water analyses"),
  Definition = c(
    "Ingestion_Rate_ng_h * 24 / BW_current_kg; time-averaged by trapezoidal integration.",
    "Infant milk-derived EDI is not directly compared with the maternal group TWI.",
    "Day-365 serum-equivalent sum of PFOA, PFNA, PFHxS and PFOS divided by 17.5 ng/mL.",
    "Compound-specific lower, central and upper ratios from the current manuscript.",
    "Not calculated: matching tissue-specific internal-dose toxicity benchmarks are unavailable.",
    "Excluded from this manuscript-aligned analysis."
  )
)
qc <- tibble(
  Item = c("Input workbook", "Scenario", "Unique infants", "Compounds",
           "Records", "Expected records", "Duplicate ID-compound-time",
           "PBPK ODE rerun", "Water data"),
  Result = c(input, "S1 only", n_distinct(dat$ID), paste(chem, collapse = ", "),
             nrow(dat), "175 x 6 x 54 (verify against supplied run)",
             sum(duplicated(dat[c("ID", "Compound", "Time_days")])),
             "No", "Not read or analysed")
)

wb <- openxlsx::createWorkbook()
sheets <- list(
  QC = qc, Methods_Notes = methods_notes,
  Milk_EDI_Individual = edi_individual, Milk_EDI_Summary = edi_summary,
  Milk_Mixture_Individual = mixture, IBQ_Individual = ibq_individual,
  IBQ_Summary = ibq_summary, IBQ_Component_Individual = ibq_component,
  IBQ_Central_Component_Shares = central_shares,
  Mass_Balance_QC = mass_balance_qc
)
for (nm in names(sheets)) {
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, as.data.frame(sheets[[nm]]), withFilter = TRUE)
  openxlsx::freezePane(wb, nm, firstRow = TRUE)
  if (ncol(sheets[[nm]]) > 0) openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(sheets[[nm]])), widths = "auto")
}
openxlsx::saveWorkbook(wb, output, overwrite = TRUE)
message("No-water risk characterization completed: ", output)
