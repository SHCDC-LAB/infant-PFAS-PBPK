# ============================================================================
# R/08_integrated_sensitivity_analysis.R
# Infant PFAS PBPK v2 | Integrated sensitivity and robustness analysis
# Date: 2026-11-04
# ============================================================================
# Scientific principles:
#   * Primary scenario is S1-proxy; original 01 script is NEVER sourced.
#   * Live ODE/physiology/renal functions are extracted from the actual 01 script.
#   * One-at-a-time (OAT) changes are hypothetical uncertainty scenarios, not CIs.
#   * The frozen main model parameters and its output files are NOT modified.
#   * Existing 03 / 06 results are linked, not mislabeled as new simulations.
#   * BW exclusion is a sample-influence check; BW reimputation is an input check.
#   * No imputed Ku for PFDA / PFUdA is silently presented as calibrated.
#   * No Shapley recalculation is claimed; 03 retains source attribution.
# ============================================================================

options(stringsAsFactors = FALSE)

need <- c("readxl", "writexl", "dplyr", "tibble", "deSolve", "ggplot2")
missing_pkgs <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install R packages before running 08: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2)
})

ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
MAIN <- file.path(ROOT, "R", "01_run_main_model.R")
INPUT <- Sys.getenv("PFAS_PBPK_INPUT", unset = file.path(ROOT, "data", "input-PBPK.xlsx"))
VERSION <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
BASE_DIR <- file.path(ROOT, "results", VERSION)
BASE_XLSX <- Sys.getenv("PFAS_08_BASE_RESULTS", unset = file.path(BASE_DIR, "tables", "PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx"))
OUT <- file.path(BASE_DIR, "integrated_sensitivity")
CACHE_DIR <- file.path(OUT, "cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------- USER CONFIGURATION ------------------------------
# full: all IDs; pilot: first N IDs only (for testing; NOT manuscript results)
RUN_MODE <- tolower(Sys.getenv("PFAS_08_RUN_MODE", "full"))
if (!RUN_MODE %in% c("full", "pilot")) stop("PFAS_08_RUN_MODE must be 'full' or 'pilot'.")
PILOT_N <- as.integer(Sys.getenv("PFAS_08_PILOT_N", "12"))
FORCE <- identical(Sys.getenv("PFAS_08_FORCE_RERUN", "0"), "1")
# Default focus is all 6 PFAS; e.g. PFAS_08_COMPOUNDS="PFNA,PFDA,PFUdA"
ALL_CHEM <- c("PFOA", "PFNA", "PFDA", "PFUdA", "PFHxS", "PFOS")
CHEM <- trimws(strsplit(Sys.getenv("PFAS_08_COMPOUNDS", paste(ALL_CHEM, collapse=",")), ",", fixed=TRUE)[[1]])
if (!length(CHEM) || any(!CHEM %in% ALL_CHEM) || anyDuplicated(CHEM)) stop("Invalid PFAS_08_COMPOUNDS")
# Optional modules can be limited when runtime is a concern.
MODULES <- trimws(strsplit(Sys.getenv("PFAS_08_MODULES", "BW,RENAL,TISSUE,MILK,GROWTH,URINE"), ",", fixed=TRUE)[[1]])
if (any(!MODULES %in% c("BW","RENAL","TISSUE","MILK","GROWTH","URINE"))) stop("Invalid PFAS_08_MODULES")
TARGET_DAYS <- c(0, 30, 90, 180, 270, 365.25)
MASS_TOL <- 1e-5
# Screening multipliers: deliberately symmetric exploration, not literature CIs.
PERT <- c(0.5, 1.5)
GFR_PERT <- c(0.8, 1.2)
MILK_PERT <- c(0.8, 1.2)
GROWTH_REMAINING <- c(0.05, 0.25, 1.00) # baseline = 0.10
# Conservative fail-closed baseline check; increasing requires explicit action.
BASE_QC_REL_TOL <- as.numeric(Sys.getenv("PFAS_08_BASELINE_TOL", "0.01"))
if (!is.finite(BASE_QC_REL_TOL) || BASE_QC_REL_TOL <= 0) stop("Invalid PFAS_08_BASELINE_TOL")

if (!file.exists(MAIN)) stop("Cannot find main 01 script: ", MAIN)
if (!file.exists(INPUT)) stop("Cannot find cohort input: ", INPUT)
if (!file.exists(BASE_XLSX)) stop("Cannot find frozen S1 result: ", BASE_XLSX)

# --------------------- SAFE, SELECTIVE MODEL LOADING -------------------------
# Parse the actual model. EVALUATE ONLY named functions and parameter declarations.
# Never source() 01: it begins with rm(list=ls()) and launches the main run.
exprs <- parse(file = MAIN, keep.source = FALSE)
model <- new.env(parent = globalenv())
model$SIM_YEARS <- 1
model$DAYS_PER_YEAR <- 365.25
model$HOURS_PER_YEAR <- 365.25 * 24
model$SIM_END_HOURS <- model$HOURS_PER_YEAR
model$TIME_STEP_HOURS <- 168
model$RENAL_CACHE_STEP_DAYS <- 1
model$MIN_VOLUME_L <- 1e-10
model$VFIL_L_PER_KG_BW <- 4.0e-4
model$TM_SCALING_EXPONENT <- 0.75
model$KU_SCALING_EXPONENT <- -0.25
model$KUST_SCALING_EXPONENT <- -0.25
model$KUST_MODE <- "QFIL_OVER_VFIL"
model$ALLOW_MISSING_KU <- TRUE
model$BW_INDIVIDUALIZATION_MODE <- "DECAY_TO_REFERENCE"
model$BW_BIRTH_EFFECT_REMAINING_AT_12MO <- 0.10
model$ODE_RTOL <- 1e-6
model$ODE_ATOL <- 1e-8
model$TRAJECTORY_MODE <- "proxy"
model$pfas_cols <- ALL_CHEM
model$INPUT_BW_UNIT <- "g"
model$INPUT_MILK_UNIT <- "ng/L"
model$DAYS_PER_MONTH <- 30.4375
model$MACROSOMIA_THRESHOLD_KG <- 4
model$EXTREME_HIGH_BW_THRESHOLD_KG <- 4.5
model$VERY_EXTREME_BW_THRESHOLD_KG <- 5
model$times <- sort(unique(c(seq(0, model$SIM_END_HOURS, by = model$TIME_STEP_HOURS), model$SIM_END_HOURS)))
model$renal_age_seq_yr <- sort(unique(c(seq(0, 1, by = 1 / 365.25), 1)))

fn_needed <- c(
  "get_reference_height_weight", "get_growth_adjustment", "get_infant_phys_meta",
  "get_infant_initial_burden", "get_tmkt_row", "scale_tmc_to_tm", "scale_ku",
  "safe_approxfun", "make_subject_physiology_cache", "make_renal_cache_unified",
  "pbpk_ode_unified_tmkt", "make_milk_function", "convert_milk_to_ng_L",
  "convert_bw_to_kg", "parse_numeric_censored", "monthly_to_beta"
)
param_needed <- c("TMKT_PARAMS", "transfer_data", "pfas_params", "trajectory_params")
seen_fn <- character(0)
get_lhs <- function(e) {
  if (!is.call(e) || !identical(e[[1]], as.name("<-"))) return(NA_character_)
  if (!is.symbol(e[[2]])) return(NA_character_)
  as.character(e[[2]])
}
# Functions first; values are resolved in model environment at call time.
for (e in exprs) {
  lhs <- get_lhs(e)
  if (!is.na(lhs) && lhs %in% fn_needed && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, envir = model)
    seen_fn <- c(seen_fn, lhs)
  }
}
if (length(setdiff(fn_needed, seen_fn))) {
  stop("Main-model definition(s) changed/missing: ", paste(setdiff(fn_needed, seen_fn), collapse=", "),
       ". Stop instead of using reimplemented equations.")
}
# Parameter declarations. trajectory_params is assigned twice by the main script;
# evaluate both expressions in original order to preserve the final beta fields.
seen_param <- character(0)
for (e in exprs) {
  lhs <- get_lhs(e)
  if (!is.na(lhs) && lhs %in% param_needed) {
    eval(e, envir = model)
    seen_param <- c(seen_param, lhs)
  }
}
if (length(setdiff(param_needed, seen_param))) stop("Main-model parameter definition missing: ", paste(setdiff(param_needed, seen_param), collapse=", "))
if (!all(c("Compound", "beta_mean_day", "supported_by_blomberg") %in% names(model$trajectory_params))) stop("Missing model trajectory columns")
if (abs(model$transfer_data$R_BM[match("PFNA", model$transfer_data$Compound)] - 0.0300) > 1e-9) stop("Unexpected PFNA frozen R_BM; model version mismatch")

# ------------------------------- INPUT -------------------------------------
raw <- as.data.frame(readxl::read_excel(INPUT, sheet = "Sheet1"), stringsAsFactors = FALSE)
need_input <- c("BW", "infsex", ALL_CHEM)
if (length(setdiff(need_input, names(raw)))) stop("Missing input: ", paste(setdiff(need_input, names(raw)), collapse = ", "))
for (chem in ALL_CHEM) raw[[chem]] <- model$parse_numeric_censored(raw[[chem]])
raw$infsex <- suppressWarnings(as.numeric(raw$infsex))
raw$BW <- suppressWarnings(as.numeric(raw$BW))
raw$BW_kg <- model$convert_bw_to_kg(raw$BW, model$INPUT_BW_UNIT)
# Mirror 01 sex-specific BW mean imputation, before row filtering.
sex_mean <- tapply(raw$BW_kg, raw$infsex, function(v) mean(v, na.rm=TRUE))
miss_bw <- !is.finite(raw$BW_kg)
raw$BW_kg[miss_bw] <- unname(sex_mean[as.character(raw$infsex[miss_bw])])
raw$BW_birth_kg <- raw$BW_kg
raw$ID <- seq_len(nrow(raw)) # initial ID mirrors 01 row_number()
raw$mix_ratio <- 1
if ("BFpc" %in% names(raw)) {
  bf <- suppressWarnings(as.numeric(raw$BFpc))
  raw$mix_ratio[is.finite(bf) & bf >= 0 & bf <= 1] <- bf[is.finite(bf) & bf >= 0 & bf <= 1]
}
for (chem in ALL_CHEM) raw[[chem]] <- model$convert_milk_to_ng_L(raw[[chem]], model$INPUT_MILK_UNIT)
keep <- is.finite(raw$BW_birth_kg) & raw$BW_birth_kg > 0 & is.finite(raw$infsex) &
  rowSums(sapply(ALL_CHEM, function(k) is.finite(raw[[k]]))) > 0
cohort <- raw[keep, , drop=FALSE]
# 01 assigns ID BEFORE filtering; replicate that faithfully.
cohort$Macrosomia_4kg_Flag <- cohort$BW_birth_kg >= 4
cohort$Extreme_High_BW_4p5kg_Flag <- cohort$BW_birth_kg >= 4.5
cohort$Very_Extreme_BW_5kg_Flag <- cohort$BW_birth_kg >= 5
if (RUN_MODE == "pilot") cohort <- head(cohort, PILOT_N)
if (!nrow(cohort)) stop("Empty cohort")
if (anyDuplicated(cohort$ID)) stop("Duplicate subject IDs")

# --------------------------- BASELINE RESULTS -------------------------------
cols_needed <- c("ID", "Compound", "Scenario", "Time_days", "C_Plasma", "C_Liver", "C_Kidney",
                 "A_Body_Burden", "A_Urine_Storage", "R_Urine_Excretion", "GFR_L_h", "Mass_Balance_Relative_Error")
headers <- names(readxl::read_excel(BASE_XLSX, sheet="PBPK_Results", n_max=0, .name_repair="minimal"))
if (length(setdiff(cols_needed, headers))) stop("Baseline workbook missing field(s): ", paste(setdiff(cols_needed, headers), collapse=", "))
base_types <- ifelse(headers %in% cols_needed, "guess", "skip")
base <- as.data.frame(readxl::read_excel(BASE_XLSX, sheet="PBPK_Results", col_types=base_types, .name_repair="minimal", progress=FALSE))
base <- base[, cols_needed, drop=FALSE]
base$ID <- as.integer(base$ID)
base$Compound <- as.character(base$Compound)
base$Scenario <- as.character(base$Scenario)
for (k in setdiff(cols_needed, c("ID", "Compound", "Scenario"))) base[[k]] <- suppressWarnings(as.numeric(base[[k]]))
base <- base[base$Scenario=="S1" & base$ID %in% cohort$ID & base$Compound %in% CHEM, , drop=FALSE]
if (!nrow(base)) stop("No matching S1 baseline data")
base_endpoint <- base[abs(base$Time_days - 365.25) < 0.001, , drop=FALSE]
expected_keys <- as.vector(outer(as.character(cohort$ID), CHEM, paste, sep="||"))
actual_keys <- paste(base_endpoint$ID, base_endpoint$Compound, sep="||")
if (anyDuplicated(actual_keys) || !setequal(actual_keys, expected_keys)) {
  stop("Baseline S1 endpoint IDs/compounds do not match cleaned input. Verify INPUT and BASE_XLSX refer to the same run.")
}
if (any(!is.finite(base_endpoint$C_Plasma) | base_endpoint$C_Plasma < 0)) stop("Invalid baseline plasma endpoint")

# -------------------- LIGHTWEIGHT DIRECT PBPK RUNNER -------------------------
# Uses frozen 01 physiology + renal caches + RHS, NOT its expensive extra source
# attribution counterfactuals. Cache and ODE calls must remain exact for base.
trapz <- function(x,y) {
  if (length(x) < 2 || any(!is.finite(x)) || any(!is.finite(y))) return(NA_real_)
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y,-1) + tail(y,-1))/2)
}
qsafe <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x,p, names=FALSE))
}

# Scenario schema: unique ID, module, variable, factor and option.
scenario <- function(id, module, parameter, multiplier=NA_real_, option="", interpretation="") {
  data.frame(Scenario_ID=id, Module=module, Parameter=parameter, Multiplier=multiplier,
             Option=option, Interpretation=interpretation, stringsAsFactors=FALSE)
}
sc <- list(scenario("BASE","BASE","BASE",1, "", "Frozen S1-proxy"))
if ("BW" %in% MODULES) {
  sc <- c(sc, list(scenario("BW_EXCL_4p5","BW","sample",NA,"exclude_ge_4.5kg","Cohort composition only"),
                   scenario("BW_REIMP_5","BW","birth_weight",NA,"reimpute_ge_5kg","Replace flagged BW using sex-specific nonflagged mean"),
                   scenario("BW_REIMP_6","BW","birth_weight",NA,"reimpute_ge_6kg","Replace flagged BW using sex-specific nonflagged mean")))
}
if ("RENAL" %in% MODULES) {
  for (p in c("Tmc", "Kt", "fu")) for (m in PERT)
    sc <- c(sc, list(scenario(paste0("RENAL_",p,"_",gsub("\\.","p",format(m))), "RENAL",p,m,
                              "", "Hypothetical OAT screening; not parameter CI")))
  for (m in GFR_PERT)
    sc <- c(sc, list(scenario(paste0("RENAL_GFR_",gsub("\\.","p",format(m))), "RENAL","GFR",m,
                              "", "Scales GFR filtration and physiology-linked filtrate residence")))
}
if ("URINE" %in% MODULES) for (m in PERT)
  sc <- c(sc, list(scenario(paste0("URINE_Ku_",gsub("\\.","p",format(m))),"URINE","Ku",m,
                            "", "Only where Ku is finite; missing PFDA/PFUdA remain unresolved")))
if ("TISSUE" %in% MODULES) for (p in c("P_li","P_kd","P_re")) for (m in PERT)
  sc <- c(sc, list(scenario(paste0("TISSUE_",p,"_",gsub("\\.","p",format(m))),"TISSUE",p,m,
                            "", "Total tissue / total plasma partition coefficient")))
if ("MILK" %in% MODULES) for (m in MILK_PERT)
  sc <- c(sc, list(scenario(paste0("MILK_Intake_",gsub("\\.","p",format(m))),"MILK","KMilk",m,
                            "", "Age-specific milk ingestion rate multiplier")))
if ("GROWTH" %in% MODULES) for (rem in GROWTH_REMAINING)
  sc <- c(sc, list(scenario(paste0("GROWTH_Rem_",gsub("\\.","p",format(rem))),"GROWTH","BW_remaining",rem,
                            "", "Proportional birth-weight deviation remaining at month 12")))
scenarios <- do.call(rbind, sc)
if (anyDuplicated(scenarios$Scenario_ID)) stop("Duplicate Scenario_ID in design")

# Frozen originals; changing local model env values does not change 01 script.
tmkt_original <- model$TMKT_PARAMS
pfas_original <- model$pfas_params
growth_original <- model$BW_BIRTH_EFFECT_REMAINING_AT_12MO

simulate_one <- function(sub, chem, spec) {
  id <- spec$Scenario_ID; par <- spec$Parameter; mul <- spec$Multiplier
  model$TMKT_PARAMS <- tmkt_original
  model$pfas_params <- pfas_original
  model$BW_BIRTH_EFFECT_REMAINING_AT_12MO <- growth_original
  on.exit({
    model$TMKT_PARAMS <- tmkt_original
    model$pfas_params <- pfas_original
    model$BW_BIRTH_EFFECT_REMAINING_AT_12MO <- growth_original
  }, add=TRUE)

  if (spec$Module == "BW" && spec$Option == "exclude_ge_4.5kg" && sub$BW_birth_kg >= 4.5) return(NULL)
  if (spec$Module == "BW" && spec$Option %in% c("reimpute_ge_5kg","reimpute_ge_6kg")) {
    cutoff <- if (spec$Option == "reimpute_ge_5kg") 5 else 6
    if (sub$BW_birth_kg >= cutoff) {
      peers <- cohort$BW_birth_kg[cohort$infsex == sub$infsex & cohort$BW_birth_kg < cutoff & is.finite(cohort$BW_birth_kg)]
      if (!length(peers)) stop("Cannot impute BW; no valid peers for sex=", sub$infsex)
      sub$BW_birth_kg <- mean(peers)
    }
  }
  if (par == "BW_remaining") model$BW_BIRTH_EFFECT_REMAINING_AT_12MO <- mul
  if (par %in% c("Tmc","Kt","fu","Ku")) {
    ix <- which(model$TMKT_PARAMS$Compound == chem)
    key <- switch(par, Tmc="Tmc_ug_h_bw075", Kt="Kt_ug_L", fu="Free_fraction", Ku="KurineC_h_bwminus025")
    current <- model$TMKT_PARAMS[[key]][ix]
    if (!is.finite(current)) return(NULL) # missing Ku => deliberately omitted
    val <- current * mul
    if (par == "fu") val <- min(1, val)
    model$TMKT_PARAMS[[key]][ix] <- val
  }
  if (par %in% c("P_li", "P_kd", "P_re")) model$pfas_params[[chem]][par] <- model$pfas_params[[chem]][par] * mul
  phys <- model$make_subject_physiology_cache(sub)
  if (par == "GFR") phys$phys_df$GFR_mL_min <- phys$phys_df$GFR_mL_min * mul
  renal <- model$make_renal_cache_unified(chem, phys$phys_df)
  f <- phys$fns
  # Physiology cache is from frozen 01. All tissue exchange and reabsorption
  # remain exactly in pbpk_ode_unified_tmkt; GFR affects Qfil, kust and Vfil.
  milk_fun <- model$make_milk_function(chem, as.numeric(sub[[chem]]), "S1")
  if (!is.finite(as.numeric(sub[[chem]])) || as.numeric(sub[[chem]]) <= 0) return(NULL)
  k_milk <- if (par == "KMilk") mul else 1
  pars <- list(
    P_li=as.numeric(model$pfas_params[[chem]]["P_li"]),
    P_kd=as.numeric(model$pfas_params[[chem]]["P_kd"]),
    P_re=as.numeric(model$pfas_params[[chem]]["P_re"]),
    mix_ratio=as.numeric(sub$mix_ratio),
    f_V_pl=f$f_V_pl, f_V_li=f$f_V_li, f_V_kd=f$f_V_kd, f_V_re=f$f_V_re,
    f_Q_li=f$f_Q_li, f_Q_kd=f$f_Q_kd, f_Q_re=f$f_Q_re,
    f_KMilk=function(age_yr) f$f_KMilk(age_yr) * k_milk,
    f_Free=renal$f_Free, f_Tm=renal$f_Tm, f_Kt=renal$f_Kt,
    f_Qfil=renal$f_Qfil, f_kust=renal$f_kust, f_ku=renal$f_ku,
    f_V_fil=renal$f_Vfil, f_C_milk=milk_fun
  )
  init <- model$get_infant_initial_burden(chem, as.numeric(sub[[chem]]),
      f$f_V_pl(0), f$f_V_li(0), f$f_V_kd(0), f$f_V_re(0),
      pars$P_li, pars$P_kd, pars$P_re)
  y0 <- c(Al=init[["Al"]], Ak=init[["Ak"]], Ar=init[["Ar"]], Ap=init[["Ap"]],
          Afil=0, Aust=0, Aurine=0, Afeces=0, Ainput=0)
  o <- deSolve::ode(y=y0, times=model$times, func=model$pbpk_ode_unified_tmkt,
                     parms=pars, method="lsoda", rtol=model$ODE_RTOL, atol=model$ODE_ATOL)
  o <- as.data.frame(o)
  tday <- o$time/24; tyr <- o$time/model$HOURS_PER_YEAR
  plasma <- o$Ap/pmax(f$f_V_pl(tyr), model$MIN_VOLUME_L)
  liver <- o$Al/pmax(f$f_V_li(tyr), model$MIN_VOLUME_L)
  kidney <- o$Ak/pmax(f$f_V_kd(tyr), model$MIN_VOLUME_L)
  burden <- o$Al + o$Ak + o$Ar + o$Ap
  all_amount <- burden + o$Afil + o$Aust + o$Aurine + o$Afeces
  mass_ref <- sum(y0[1:8]) + o$Ainput
  mass_err <- max(abs(all_amount-mass_ref)/pmax(abs(mass_ref),1e-18))
  if (any(!is.finite(c(plasma, liver, kidney, burden))) ||
      any(c(plasma, liver, kidney, burden) < -1e-8) ||
      !is.finite(mass_err) || mass_err > MASS_TOL) {
    stop("Numerical QC failed: ", id, ", ID=", sub$ID, ", chem=", chem,
         ", max relative mass error=", signif(mass_err,3))
  }
  idx <- which.min(abs(tday-365.25))
  if (abs(tday[idx]-365.25) > 1e-4) stop("Missing exact 365.25-day endpoint")
  d7 <- which.min(abs(tday-7))
  ku <- renal$f_ku(tyr[d7])
  urinary_rate <- if (is.finite(ku)) ku*o$Aust[d7] else NA_real_
  # Amount rate ng/h divided by plasma ng/L -> L/h, then mL/kg/day.
  cl_ur <- if (is.finite(urinary_rate) && plasma[d7] > 0) urinary_rate/plasma[d7]*1000*24/f$f_BW(tyr[d7]) else NA_real_
  out <- data.frame(
    Scenario_ID=id, ID=as.integer(sub$ID), Compound=chem,
    BW_birth_kg=as.numeric(sub$BW_birth_kg),
    C_Plasma_12m=plasma[idx], C_Liver_12m=liver[idx], C_Kidney_12m=kidney[idx],
    A_Body_Burden_12m=burden[idx], AUC_Plasma_0_12m=trapz(tday,plasma),
    CL_Urine_d7_mL_kg_day=cl_ur, A_Urine_Storage_12m=o$Aust[idx],
    Max_Mass_Error=mass_err, stringsAsFactors=FALSE
  )
  # Actual model reporting times, not invented exact day 30/90 etc.
  selected <- lapply(TARGET_DAYS, function(day) {
    i <- which.min(abs(tday-day))
    data.frame(Scenario_ID=id, ID=as.integer(sub$ID), Compound=chem,
               Target_Day=day, Actual_Day=tday[i],
               C_Plasma=plasma[i], C_Liver=liver[i], C_Kidney=kidney[i],
               stringsAsFactors=FALSE)
  })
  list(endpoint=out, selected=do.call(rbind,selected))
}

# -------------------------- BASELINE CROSS-CHECK ----------------------------
# Compare independently integrated ODE directly to published 01 baseline.
# For any mismatch > 1% stop, rather than generate a plausible but wrong Excel.
check_ids <- head(cohort$ID, 2L)
check_compounds <- CHEM
check_rows <- list()
for (pid in check_ids) for (chem in check_compounds) {
  sub <- cohort[cohort$ID==pid, , drop=FALSE]
  sim <- simulate_one(sub, chem, scenarios[1,,drop=FALSE])
  ref <- base_endpoint[base_endpoint$ID==pid & base_endpoint$Compound==chem, ,drop=FALSE]
  if (is.null(sim) || nrow(ref)!=1) stop("Baseline QC cannot locate unique ID/compound")
  for (nm in c("C_Plasma","C_Liver","C_Kidney","A_Body_Burden")) {
    pred <- sim$endpoint[[paste0(nm,"_12m")]]
    obs <- as.numeric(ref[[nm]])
    rel <- abs(pred-obs)/pmax(abs(obs),abs(pred), 1e-8)
    check_rows[[length(check_rows)+1]] <- data.frame(ID=pid,Compound=chem,Variable=nm,
                       Recomputed=pred,Frozen=obs,Relative_Difference=rel)
  }
}
check_df <- do.call(rbind,check_rows)
utils::write.csv(check_df, file.path(OUT,"baseline_reproduction_QC.csv"), row.names=FALSE)
if (any(!is.finite(check_df$Relative_Difference)) || any(check_df$Relative_Difference > BASE_QC_REL_TOL)) {
  stop("BASELINE REPRODUCTION FAILED. Compare baseline_reproduction_QC.csv. ",
       "Do not use this sensitivity result until model 01 and input versions match.")
}
message("08 baseline reproduction QC passed; max relative discrepancy = ", signif(max(check_df$Relative_Difference),4))

# ----------------------------- EXISTING RESULTS -----------------------------
# Reuse existing 01 baseline and 03 / 06 provenance; never relabel 03/06 as new.
# Summary for frozen baseline derived directly from its individual 01 records.
make_baseline_endpoints <- function() {
  groups <- split(base, interaction(base$ID, base$Compound, drop=TRUE))
  rows <- lapply(groups, function(d) {
    d <- d[order(d$Time_days),,drop=FALSE]
    last <- d[which.min(abs(d$Time_days-365.25)),,drop=FALSE]
    data.frame(Scenario_ID="BASE", ID=last$ID, Compound=last$Compound,
       BW_birth_kg=cohort$BW_birth_kg[match(last$ID,cohort$ID)],
       C_Plasma_12m=last$C_Plasma,C_Liver_12m=last$C_Liver,
       C_Kidney_12m=last$C_Kidney,A_Body_Burden_12m=last$A_Body_Burden,
       AUC_Plasma_0_12m=trapz(d$Time_days,d$C_Plasma),
       CL_Urine_d7_mL_kg_day=NA_real_,
       A_Urine_Storage_12m=last$A_Urine_Storage,
       Max_Mass_Error=max(abs(d$Mass_Balance_Relative_Error),na.rm=TRUE))
  })
  do.call(rbind,rows)
}
base_ep <- make_baseline_endpoints()
# Day-7 urinary comparison uses the model's first near-7d reported time only;
# for accurate day-7 clearances use purpose-built layer-2 diagnostics (02).

# ------------------------------- SCENARIOS ---------------------------------
results_ep <- list(BASE=base_ep)
results_day <- list()
status_rows <- list()
for (j in seq_len(nrow(scenarios))) {
  spec <- scenarios[j,,drop=FALSE]
  sid <- spec$Scenario_ID
  if (sid == "BASE") next
  cache_rds <- file.path(CACHE_DIR, paste0(sid,"_",RUN_MODE,"_",paste(CHEM,collapse="-"),".rds"))
  # Reuse only exact design/preset/subject signature; other sessions are unsafe.
  signature <- list(input_md5=unname(tools::md5sum(INPUT)),main_md5=unname(tools::md5sum(MAIN)),
                    base_md5=unname(tools::md5sum(BASE_XLSX)),
                    subjects=cohort$ID,chem=CHEM,scenario=as.list(spec),run_mode=RUN_MODE)
  cached <- NULL
  if (!FORCE && file.exists(cache_rds)) {
    obj <- tryCatch(readRDS(cache_rds),error=function(e) NULL)
    if (is.list(obj) && identical(obj$signature,signature)) cached <- obj
  }
  if (!is.null(cached)) {
    ep <- cached$endpoint; days <- cached$selected
    action <- "reused_checkpoint"
  } else if (sid == "BW_EXCL_4p5") {
    ep <- base_ep[base_ep$BW_birth_kg < 4.5, , drop=FALSE]
    ep$Scenario_ID <- sid
    days <- data.frame()
    saveRDS(list(signature=signature,endpoint=ep,selected=days),cache_rds)
    action <- "filtered_frozen_baseline_no_ODE"
  } else {
    message(sprintf("[%02d/%02d] Simulating %s; %s subjects x %s compounds",
                    j,nrow(scenarios),sid,nrow(cohort),length(CHEM)))
    ep_list <- list(); day_list <- list(); kk <- 0L
    for (i in seq_len(nrow(cohort))) {
      sub <- cohort[i,,drop=FALSE]
      for (chem in CHEM) {
        sim <- simulate_one(sub,chem,spec)
        if (is.null(sim)) next
        kk <- kk+1L
        ep_list[[kk]] <- sim$endpoint
        day_list[[kk]] <- sim$selected
      }
      if (i %% max(1,as.integer(nrow(cohort)/5)) == 0 || i == nrow(cohort))
        message("  ",sid,": ",i,"/",nrow(cohort)," subjects")
    }
    ep <- if(length(ep_list)) do.call(rbind,ep_list) else base_ep[FALSE,,drop=FALSE]
    days <- if(length(day_list)) do.call(rbind,day_list) else data.frame()
    saveRDS(list(signature=signature, endpoint=ep, selected=days),cache_rds)
    action <- "calculated"
  }
  results_ep[[sid]] <- ep
  results_day[[sid]] <- days
  status_rows[[length(status_rows)+1]] <- data.frame(Scenario_ID=sid,Module=spec$Module,
      Parameter=spec$Parameter,Multiplier=spec$Multiplier,N_Rows=nrow(ep),
      N_Subjects=length(unique(ep$ID)),Action=action,Cache=cache_rds)
  utils::write.csv(do.call(rbind,status_rows),file.path(OUT,"08_run_status.csv"),row.names=FALSE)
}

# --------------------- ROBUSTNESS / SUMMARY ---------------------------------
all_ep <- dplyr::bind_rows(results_ep)
measures <- c("C_Plasma_12m","C_Liver_12m","C_Kidney_12m",
              "A_Body_Burden_12m","AUC_Plasma_0_12m","A_Urine_Storage_12m")
summary_rows <- list()
for (sid in unique(all_ep$Scenario_ID)) for (chem in CHEM) {
  d <- all_ep[all_ep$Scenario_ID==sid & all_ep$Compound==chem, ,drop=FALSE]
  for (metric in measures) {
    v <- d[[metric]]
    summary_rows[[length(summary_rows)+1]] <- data.frame(
      Scenario_ID=sid, Compound=chem, Metric=metric, N=sum(is.finite(v)),
      Median=qsafe(v,.5),P025=qsafe(v,.025),P975=qsafe(v,.975))
  }
}
sum_df <- do.call(rbind,summary_rows)
ref_med <- sum_df[sum_df$Scenario_ID=="BASE",c("Compound","Metric","Median")]
names(ref_med)[3] <- "Baseline_Median"
comparison <- merge(sum_df,ref_med,by=c("Compound","Metric"),all.x=TRUE,sort=FALSE)
comparison$Change_pct <- ifelse(is.finite(comparison$Baseline_Median) & abs(comparison$Baseline_Median)>1e-12,
    100*(comparison$Median/comparison$Baseline_Median - 1), NA_real_)
comparison$Module <- scenarios$Module[match(comparison$Scenario_ID, scenarios$Scenario_ID)]
comparison$Parameter <- scenarios$Parameter[match(comparison$Scenario_ID, scenarios$Scenario_ID)]
comparison$Multiplier <- scenarios$Multiplier[match(comparison$Scenario_ID, scenarios$Scenario_ID)]
comparison$Interpretation <- scenarios$Interpretation[match(comparison$Scenario_ID, scenarios$Scenario_ID)]

# Paired cohort change, defined per ID; BW_EXCL is a cohort-composition check
# and intentionally must not be treated as a paired counterfactual intervention.
paired <- list()
base_for_join <- base_ep[,c("ID","Compound",measures),drop=FALSE]
names(base_for_join)[-(1:2)] <- paste0("BASE_",measures)
for (sid in setdiff(unique(all_ep$Scenario_ID),"BASE")) {
  if (sid=="BW_EXCL_4p5") next
  t <- all_ep[all_ep$Scenario_ID==sid,,drop=FALSE]
  if (!nrow(t)) next
  z <- merge(t,base_for_join,by=c("ID","Compound"))
  for (metric in measures) {
    x <- z[[metric]]; y <- z[[paste0("BASE_",metric)]]
    good <- is.finite(x)&is.finite(y)&abs(y)>1e-12
    if (!any(good)) next
    d <- 100*(x[good]/y[good]-1)
    paired[[length(paired)+1]] <- data.frame(Scenario_ID=sid,Compound=z$Compound[good],
      ID=z$ID[good],Metric=metric,Paired_Change_pct=d)
  }
}
paired_df <- if(length(paired)) do.call(rbind,paired) else data.frame()
paired_summary <- if(nrow(paired_df)) {
  paired_df %>% group_by(Scenario_ID,Compound,Metric) %>%
    summarise(N=n(),Median_Paired_Change_pct=median(Paired_Change_pct),
      P025=qsafe(Paired_Change_pct,.025),P975=qsafe(Paired_Change_pct,.975),.groups="drop")
} else data.frame()

# Selected-age compartment trajectories: aggregate actual PBPK outputs.
# Baseline is taken from frozen 01, not from new approximations.
base_selected_rows <- lapply(TARGET_DAYS, function(day) {
  all_keys <- split(base,interaction(base$ID,base$Compound,drop=TRUE))
  do.call(rbind,lapply(all_keys,function(d) {
    k <- which.min(abs(d$Time_days-day))
    data.frame(Scenario_ID="BASE", ID=d$ID[k], Compound=d$Compound[k],
      Target_Day=day,Actual_Day=d$Time_days[k],
      C_Plasma=d$C_Plasma[k],C_Liver=d$C_Liver[k],C_Kidney=d$C_Kidney[k])
  }))
})
selected_parts <- c(list(do.call(rbind,base_selected_rows)),results_day)
selected_parts <- selected_parts[vapply(selected_parts,function(z) is.data.frame(z) && nrow(z)>0,logical(1))]
selected_all <- do.call(rbind,selected_parts)
selected_metrics <- c("C_Plasma","C_Liver","C_Kidney")
selected_summary_rows <- list()
for (sid in unique(selected_all$Scenario_ID)) for (chem in CHEM) for (day in TARGET_DAYS) {
  x <- selected_all[selected_all$Scenario_ID==sid & selected_all$Compound==chem &
                     abs(selected_all$Target_Day-day)<0.001,,drop=FALSE]
  if(!nrow(x)) next
  for (metric in selected_metrics) selected_summary_rows[[length(selected_summary_rows)+1]] <-
     data.frame(Scenario_ID=sid,Compound=chem,Target_Day=day,Metric=metric,
       N=sum(is.finite(x[[metric]])),Median=qsafe(x[[metric]],.5),
       P025=qsafe(x[[metric]],.025),P975=qsafe(x[[metric]],.975),
       Actual_Day_Median=qsafe(x$Actual_Day,.5))
}
selected_summary <- do.call(rbind,selected_summary_rows)
selected_base <- selected_summary[selected_summary$Scenario_ID=="BASE",
                  c("Compound","Target_Day","Metric","Median")]
names(selected_base)[4] <- "Baseline_Median"
selected_comparison <- merge(selected_summary,selected_base,
                     by=c("Compound","Target_Day","Metric"),all.x=TRUE,sort=FALSE)
selected_comparison$Change_pct <- ifelse(is.finite(selected_comparison$Baseline_Median) &
                          abs(selected_comparison$Baseline_Median)>1e-12,
                          100*(selected_comparison$Median/selected_comparison$Baseline_Median-1),NA_real_)
utils::write.csv(selected_comparison,file.path(OUT,"08_selected_days_relative_changes.csv"),row.names=FALSE)

# Local normalized finite-difference indices, only for continuous multipliers;
# do not compute artificial elasticity for BW exclusion/reimputation or growth.
elasticity <- comparison[comparison$Module %in% c("RENAL","TISSUE","MILK","URINE") &
                           comparison$Scenario_ID!="BASE" & is.finite(comparison$Multiplier) &
                           comparison$Multiplier!=1,,drop=FALSE]
elasticity$Elasticity <- (elasticity$Change_pct/100)/(elasticity$Multiplier - 1)
# Rank by maximum magnitude across the tested side-specific perturbations.
ranking <- if(nrow(elasticity)) {
  elasticity %>% group_by(Compound,Metric,Parameter) %>%
    summarise(MaxAbsElasticity=max(abs(Elasticity), na.rm=TRUE),.groups="drop") %>%
    arrange(Compound,Metric,desc(MaxAbsElasticity))
} else data.frame()

# Scenario-level QC: no missing endpoints unless Ku was explicitly undefined.
qcs <- do.call(rbind,lapply(names(results_ep),function(sid) {
  x <- results_ep[[sid]]
  data.frame(Scenario_ID=sid, N_Rows=nrow(x), N_Infants=length(unique(x$ID)),
             Max_Mass_Error=if(nrow(x)) max(x$Max_Mass_Error,na.rm=TRUE) else NA_real_,
             Negative_Concentrations=if(nrow(x)) sum(as.matrix(x[,c("C_Plasma_12m","C_Liver_12m","C_Kidney_12m")]) < -1e-8,na.rm=TRUE) else 0)
}))

# 03 and 06 are reported as existing analysis dependencies, not re-run here.
existing <- data.frame(Module=c("03 source attribution", "06 beta trajectory"),
  File=c(file.path(BASE_DIR,"source_attribution_sensitivity","PFAS_RCM_RBM_source_attribution_sensitivity_v1_EP.xlsx"),
         file.path(ROOT,"results","beta_proxy_sensitivity_summary.xlsx")),stringsAsFactors=FALSE)
existing$Exists <- file.exists(existing$File)
existing$Note <- c("Separate nonlinear Shapley simulations; never recomputed by 08",
                   "Already-completed proxy/flat/half/strong simulations; not duplicated by 08")
# Import the existing real result sheets as addenda if their workbooks exist.
# Do not silently manufacture missing estimates or interpolate external inputs.
read_optional_sheet <- function(path, sheet) {
  if (!file.exists(path)) return(data.frame())
  available <- tryCatch(readxl::excel_sheets(path), error=function(e) character())
  if (!sheet %in% available) return(data.frame())
  x <- tryCatch(as.data.frame(readxl::read_excel(path, sheet=sheet)), error=function(e) data.frame())
  if (!nrow(x)) return(data.frame())
  x
}
existing_03 <- read_optional_sheet(existing$File[1],"Baseline_Formal_Summary")
existing_06 <- read_optional_sheet(existing$File[2],"Relative_To_Proxy")
# Explicitly indicate which actual results have been imported.
existing$Imported <- c(nrow(existing_03)>0, nrow(existing_06)>0)

# -------------------------- EXPORT RESULTS ---------------------------------
utils::write.csv(comparison,file.path(OUT,"08_sensitivity_relative_changes.csv"),row.names=FALSE)
if(nrow(paired_df)) utils::write.csv(paired_df,file.path(OUT,"08_paired_infant_changes.csv"),row.names=FALSE)
if(nrow(ranking)) utils::write.csv(ranking,file.path(OUT,"08_elasticity_ranking.csv"),row.names=FALSE)

interpretation <- data.frame(
  Topic=c("Model basis","Baseline reproduction","BW exclusion","BW reimputation","Renal perturbations",
          "Partition coefficients","Missing urinary Ku","GFR","Urinary day 7","03 and 06", "External validation",
          "Study scope"),
  Note=c(
    "Exact definitions extracted from 01 without sourcing the executable main script; S1-proxy only",
    "Sample-level model reproduction must pass before perturbation analyses are written",
    "BW_EXCL_4p5 changes sample composition, NOT the same individuals' physiology",
    "BW_REIMP_5 and BW_REIMP_6 modify BW for >= threshold using nonflagged within-sex mean; not silent winsorizing",
    "Tmc, Kt, free fraction and GFR are screened independently; multipliers are assumed values, NOT confidence limits",
    "Pli/Pkd/Pre represent TOTAL tissue / TOTAL plasma concentration ratios. Filtration alone uses free fraction.",
    "PFDA/PFUdA Ku remains NA; no quantitative urinary outlet clearance is interpreted",
    "GFR is changed downstream of organ perfusion and BSA; this is partial uncertainty, not whole-body physiology rescaling",
    "Day 7 of weekly grid may be approximate; prefer dedicated model 02 renal diagnostics for Yao comparisons",
    "03 Shapley and 06 beta ensemble are independently existing modules, not newly rerun here",
    "Parameter sensitivity alone cannot validate infant organ or plasma concentrations",
    if(RUN_MODE=="pilot") "PILOT ONLY; cannot be used as manuscript full-cohort results" else "FULL cohort, provided input IDs match frozen baseline"
  ),stringsAsFactors=FALSE)

wb <- list(Design=scenarios, Run_Status=if(length(status_rows))do.call(rbind,status_rows) else data.frame(),
           Baseline_Repro_QC=check_df, Scenario_QC=qcs,
           Endpoint_Summary=sum_df, Relative_To_Baseline=comparison,
           Selected_Days=selected_comparison, Paired_Change_Summary=paired_summary,
           Elasticity_Ranking=ranking,
           Existing_03_06=existing, Attribution_03=existing_03,
           Milk_Beta_06=existing_06, Interpretation=interpretation)
# writexl sheet names never exceed 31 chars.
wb <- lapply(wb,function(d) if(ncol(d)==0) data.frame(Information="No applicable rows") else d)
wb <- lapply(wb,function(d) { d[] <- lapply(d,function(x) if(is.factor(x)) as.character(x) else x); d })
workbook <- file.path(OUT,"integrated_sensitivity_analysis.xlsx")
writexl::write_xlsx(wb,workbook)

# Journal-supporting figure: normalized plasma endpoint across modules.
plot_df <- comparison[comparison$Scenario_ID!="BASE" & comparison$Metric=="C_Plasma_12m" &
                        is.finite(comparison$Change_pct),,drop=FALSE]
if (nrow(plot_df)) {
  p <- ggplot(plot_df, aes(x=reorder(Scenario_ID,Change_pct),y=Change_pct)) +
    geom_col(fill="#426E9D",width=.7) + coord_flip() + facet_wrap(~Compound,scales="free_y",ncol=3) +
    geom_hline(yintercept=0,linewidth=.25) + labs(x=NULL,y="Change vs frozen S1-proxy (%)",
      title="Sensitivity of 12-month plasma PFAS concentrations") +
    theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())
  ggsave(file.path(OUT,"Figure_08_Plasma_Sensitivity.png"),p,width=13,height=10,dpi=350)
  ggsave(file.path(OUT,"Figure_08_Plasma_Sensitivity.pdf"),p,width=13,height=10)
}

message("\n08 INTEGRATED SENSITIVITY COMPLETE")
message("  Run mode: ",RUN_MODE)
message("  Subjects: ",nrow(cohort),"; compounds: ",paste(CHEM,collapse=", "))
message("  Workbook: ",workbook)
message("  Checkpoints: ",CACHE_DIR)
message("  Separate 03/06 outputs were not overwritten.")
