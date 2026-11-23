# ============================================================================
# R/10_global_sensitivity_uncertainty_FAST175_PARALLEL.R
# Full-cohort (175 infants) adaptive batched LHS / PRCC / conditional
# parameter-range uncertainty and PFUdA Tmc x Kt interaction.
# IMPORTANT: uses the frozen 01 ODE functions; does NOT source the run script.
# Authoring version: 2026-11-04 | parameter ranges HYPOTHETICAL, not CIs.
# PRCC = partial rank correlation, NOT contribution to explained variance.
# ============================================================================
options(stringsAsFactors = FALSE)
needed <- c("readxl","writexl","dplyr","tibble","deSolve","ggplot2")
missing_packages <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)) stop("Install R packages: ",paste(missing_packages,collapse=", "))
suppressPackageStartupMessages({library(dplyr);library(tibble);library(ggplot2)})
ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT",unset=getwd()),winslash="/",mustWork=TRUE)
MAIN <- file.path(ROOT,"R","01_run_main_model.R")
INPUT <- Sys.getenv("PFAS_PBPK_INPUT",unset=file.path(ROOT,"data","input-PBPK.xlsx"))
VERSION <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
BASE_XLSX <- Sys.getenv("PFAS_10_BASE_RESULTS",unset=file.path(ROOT,"results",VERSION,"tables","PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx"))
# Reuse original results/cache location so previously completed PFOA draws are retained.
OUT <- file.path(ROOT,"results",VERSION,"global_sensitivity_FAST175")
WORKERS <- as.integer(Sys.getenv("PFAS_10_WORKERS", "4"))
if(is.na(WORKERS) || WORKERS < 1L || WORKERS > 24L) stop("PFAS_10_WORKERS must be 1..24")
CACHE_EVERY <- as.integer(Sys.getenv("PFAS_10_CACHE_EVERY", "5"))
if(is.na(CACHE_EVERY)||CACHE_EVERY < 1L) stop("PFAS_10_CACHE_EVERY must be >=1")
ADOPT_CACHE <- identical(Sys.getenv("PFAS_10_ADOPT_EXISTING_CACHE", "0"), "1")
RUN_MODE <- tolower(Sys.getenv("PFAS_10_RUN_MODE","pilot"))
if(!RUN_MODE %in% c("pilot","full")) stop("PFAS_10_RUN_MODE must be pilot or full")
PILOT_N <- as.integer(Sys.getenv("PFAS_10_PILOT_N","12"))
ALL_CHEM <- c("PFOA","PFNA","PFDA","PFUdA","PFHxS","PFOS")
CHEM <- trimws(strsplit(Sys.getenv("PFAS_10_COMPOUNDS",paste(ALL_CHEM,collapse=",")),",",fixed=TRUE)[[1]])
if(!length(CHEM)||any(!CHEM %in% ALL_CHEM)||anyDuplicated(CHEM)) stop("Invalid compounds")
SEED <- as.integer(Sys.getenv("PFAS_10_SEED","20261104"))
FORCE <- identical(Sys.getenv("PFAS_10_FORCE_RERUN","0"),"1")
# Each stage adds NEW independent stratified LHS blocks; accumulated design is
# blockwise LHS (not mathematically a single nested LHS of final size).
STAGES <- as.integer(trimws(strsplit(Sys.getenv("PFAS_10_STAGES",if(RUN_MODE=="pilot") "20,40" else "100,150,200,300,400,600"),",",fixed=TRUE)[[1]]))
if(anyNA(STAGES)||any(diff(STAGES)<=0)||STAGES[1]<15L) stop("Stages must be increasing integers >=15")
MAX_DRAWS <- tail(STAGES,1)
BATCH_SIZE <- as.integer(Sys.getenv("PFAS_10_BATCH_SIZE",if(RUN_MODE=="pilot") "10" else "25"))
STABLE_PRCC_TOL <- as.numeric(Sys.getenv("PFAS_10_PRCC_TOL","0.05"))
STABLE_Q_TOL <- as.numeric(Sys.getenv("PFAS_10_QUANTILE_TOL","0.05"))
MIN_CONVERGENCE_N <- as.integer(Sys.getenv("PFAS_10_MIN_CONVERGENCE_N",if(RUN_MODE=="pilot") "40" else "300"))
# 2 consecutive eligible comparisons: does not stop at 300 just because one comparison agrees.
N_STABLE_CHECKS <- as.integer(Sys.getenv("PFAS_10_STABLE_CHECKS","2"))
TOP_K <- as.integer(Sys.getenv("PFAS_10_TOP_K","5"))
BASE_QC_REL_TOL <- as.numeric(Sys.getenv("PFAS_10_BASELINE_TOL","0.01"))
DO_INTERACTION <- Sys.getenv("PFAS_10_INTERACTION",if(RUN_MODE=="pilot") "0" else "1") == "1"
if(!is.finite(SEED)||!is.finite(BATCH_SIZE)||BATCH_SIZE<1||!is.finite(MIN_CONVERGENCE_N)||
   !is.finite(N_STABLE_CHECKS)||N_STABLE_CHECKS<1||!is.finite(TOP_K)||TOP_K<1||
   !is.finite(BASE_QC_REL_TOL)||BASE_QC_REL_TOL<=0||!is.finite(STABLE_Q_TOL)||STABLE_Q_TOL<=0||
   !is.finite(STABLE_PRCC_TOL)||STABLE_PRCC_TOL<=0) stop("Invalid FAST175 options")
if(!file.exists(MAIN)||!file.exists(INPUT)||!file.exists(BASE_XLSX)) stop(
  "Missing MAIN/INPUT/BASE_XLSX; set PFAS_PBPK_ROOT or PFAS_10_BASE_RESULTS")
# Distinct pilot/full dirs. Never overwrite 01-09 or previous 10 outputs.
OUT <- file.path(OUT,RUN_MODE)
CACHE_DIR <- file.path(OUT,"cache")
dir.create(CACHE_DIR,recursive=TRUE,showWarnings=FALSE)

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


# ---------------------------- FREEZE CHECKS --------------------------------
# Fail closed: an old 01 can silently produce plausible but wrong kidney curves.
main_text <- paste(readLines(MAIN,warn=FALSE,encoding="UTF-8"),collapse="\n")
if(!grepl("dAk\\s*<-\\s*Q_kd\\s*\\*\\s*\\(\\s*Cp\\s*-\\s*Ck\\s*/\\s*P_kd\\s*\\)",main_text,perl=TRUE))
  stop("Frozen 01 kidney correction missing: dAk must use Cp - Ck/P_kd (no fu in tissue exchange).")
if(any(!is.finite(model$trajectory_params$beta_mean_day))) stop("Missing beta_mean_day")
if(any(abs(model$trajectory_params$beta_mean_day[match(c("PFNA","PFDA","PFUdA"),model$trajectory_params$Compound)])<1e-12))
  stop("The 01 script is not the S1-proxy version: PFNA/PFDA/PFUdA beta must be nonzero")
if(!identical(as.character(model$TRAJECTORY_MODE),"proxy")) stop("S1 proxy trajectory mode required")
if(RUN_MODE=="full" && nrow(cohort)!=175L)
  stop("FAST175 full mode requires exactly 175 retained infants; got ",nrow(cohort))
if(RUN_MODE=="pilot") message("PILOT ONLY: ",nrow(cohort)," infants; NOT publishable full-cohort results")

# ------------------------- PARAMETER SCREENING -----------------------------
# Log-uniform is a deliberately bounded screening distribution; it is NOT
# an evidence-calibrated probabilistic uncertainty distribution.
PARAMS <- data.frame(
  Parameter=c("Milk_Intake","R_CM","R_BM","Tmc","Kt","fu","GFR","P_li","P_kd","P_re","Beta_Proxy"),
  Low=c(.8,.7,.7,.5,.5,.7,.8,.5,.5,.5,0),
  High=c(1.2,1.3,1.3,1.5,1.5,1.3,1.2,1.5,1.5,1.5,1.5),
  Scale=c(rep("log_uniform",10),"uniform"),stringsAsFactors=FALSE)
stopifnot(nrow(PARAMS) == 11L, !anyNA(PARAMS$Low), !anyNA(PARAMS$High),
          all(PARAMS$Low < PARAMS$High),
          all(PARAMS$Low[PARAMS$Scale == "log_uniform"] > 0))
PARAMS$Evidence <- "Hypothetical parameter screening range (NOT empirical CI)"
PARAM_NAMES <- PARAMS$Parameter
TARGET_DAYS <- c(0,7,30,180,365.25)
TIME_VECTOR <- model$times
TARGET_ACTUAL_DAY <- vapply(TARGET_DAYS,function(d){
  x<-TIME_VECTOR/24; x[which.min(abs(x-d))]
},numeric(1))
METRICS <- c("C_Plasma_12m","C_Plasma_6m","AUC_Plasma_0_12m","C_Liver_12m","C_Kidney_12m")

qfast <- function(x,p) {
  x<-x[is.finite(x)]
  if(!length(x)) return(NA_real_)
  unname(stats::quantile(x,p,type=7))
}
trapfast <- function(x,y){
  if(any(!is.finite(y))||length(x)!=length(y)) return(NA_real_)
  sum(diff(x)*(head(y,-1L)+tail(y,-1L))*.5)
}
atomic_rds <- function(value,path){
  tmp<-paste0(path,".writing")
  saveRDS(value,tmp,compress=TRUE)
  if(!file.rename(tmp,path)){
    if(!file.copy(tmp,path,overwrite=TRUE)) stop("Cannot write checkpoint: ",path)
    unlink(tmp)
  }
}
# Use separate stable LHS blocks for every interval between requested stages.
# Identical design used across all six PFAS; Beta_Proxy inactive for 3 with literature slopes.
make_block_design <- function(n,seed){
  set.seed(seed)
  out<-matrix(NA_real_,n,length(PARAM_NAMES),dimnames=list(NULL,PARAM_NAMES))
  for(j in seq_along(PARAM_NAMES)){
    u<-(sample.int(n)-1+runif(n))/n
    lo<-PARAMS$Low[j];hi<-PARAMS$High[j]
    out[,j]<-if(PARAMS$Scale[j]=="log_uniform")
      exp(log(lo)+u*(log(hi)-log(lo))) else lo+u*(hi-lo)
  }
  as.data.frame(out,check.names=FALSE)
}
blocks<-diff(c(0L,STAGES))
design_list<-lapply(seq_along(blocks),function(b)make_block_design(blocks[b],SEED+b*10007L))
DESIGN<-do.call(rbind,design_list)
DESIGN$Draw<-seq_len(nrow(DESIGN))
DESIGN<-DESIGN[,c("Draw",PARAM_NAMES),drop=FALSE]
stopifnot(nrow(DESIGN)==MAX_DRAWS)
utils::write.csv(DESIGN,file.path(OUT,"FAST175_LHS_design.csv"),row.names=FALSE)

# Prevent unsafe reuse when code/input/freeze/seed/selection/metrics changes.
# Readable text signature includes complete numeric design and MD5s.
script_path <- file.path(ROOT,"R","10_global_sensitivity_uncertainty_FAST175_PARALLEL.R")
script_md5 <- if(file.exists(script_path)) unname(tools::md5sum(script_path)) else "external_source_file"
FINGERPRINT <- list(
  engine="FAST175_v2",main_md5=unname(tools::md5sum(MAIN)),
  input_md5=unname(tools::md5sum(INPUT)),
  baseline_md5=unname(tools::md5sum(BASE_XLSX)),
  design=DESIGN,cohort_ids=cohort$ID,compounds=CHEM,
  model_rtol=model$ODE_RTOL,model_atol=model$ODE_ATOL,
  script_md5=script_md5
)

# Pre-compute expensive age-dependent physiology ONCE for each infant.
# This is the principal exact (non-emulator) speed improvement relative to 10.
message("Precomputing infant physiology for ",nrow(cohort)," subjects...")
PHYS<-setNames(lapply(seq_len(nrow(cohort)),function(i)
  model$make_subject_physiology_cache(cohort[i,,drop=FALSE])),as.character(cohort$ID))

# Store baseline params and restore them on every simulation, for safety.
TMKT0<-model$TMKT_PARAMS
PFAS0<-model$pfas_params
TRANSFER0<-model$transfer_data
TRAJ0<-model$trajectory_params
unit_multipliers <- function() as.list(setNames(rep(1,length(PARAM_NAMES)),PARAM_NAMES))

sim_fast175 <- function(sub,chem,mult,with_curve=FALSE){
  # Model functions execute in the isolated 'model' environment.
  model$TMKT_PARAMS<-TMKT0
  model$pfas_params<-PFAS0
  model$transfer_data<-TRANSFER0
  model$trajectory_params<-TRAJ0
  on.exit({model$TMKT_PARAMS<-TMKT0;model$pfas_params<-PFAS0;
  model$transfer_data<-TRANSFER0;model$trajectory_params<-TRAJ0},add=TRUE)
  ix<-match(chem,model$TMKT_PARAMS$Compound)
  for(pair in list(c("Tmc","Tmc_ug_h_bw075"),c("Kt","Kt_ug_L"),c("fu","Free_fraction")))
    model$TMKT_PARAMS[[pair[2]]][ix]<-model$TMKT_PARAMS[[pair[2]]][ix]*mult[[pair[1]]]
  model$TMKT_PARAMS$Free_fraction[ix]<-min(1,max(1e-6,model$TMKT_PARAMS$Free_fraction[ix]))
  for(p in c("P_li","P_kd","P_re")) model$pfas_params[[chem]][p]<-model$pfas_params[[chem]][p]*mult[[p]]
  tx<-match(chem,model$transfer_data$Compound)
  model$transfer_data$R_CM[tx]<-model$transfer_data$R_CM[tx]*mult$R_CM
  model$transfer_data$R_BM[tx]<-model$transfer_data$R_BM[tx]*mult$R_BM
  if(!isTRUE(TRAJ0$supported_by_blomberg[match(chem,TRAJ0$Compound)])){
    k<-match(chem,model$trajectory_params$Compound)
    model$trajectory_params$beta_mean_day[k]<-model$trajectory_params$beta_mean_day[k]*mult$Beta_Proxy
  }
  phy<-PHYS[[as.character(sub$ID)]]
  p_df<-phy$phys_df
  if(mult$GFR!=1){p_df$GFR_mL_min<-p_df$GFR_mL_min*mult$GFR}
  renal<-model$make_renal_cache_unified(chem,p_df)
  f<-phy$fns
  c0<-as.numeric(sub[[chem]])
  if(!is.finite(c0)||c0<=0) stop("Nonpositive milk C0: ",chem," ID ",sub$ID)
  milk<-model$make_milk_function(chem,c0,"S1")
  pars<-list(
    P_li=as.numeric(model$pfas_params[[chem]]["P_li"]),
    P_kd=as.numeric(model$pfas_params[[chem]]["P_kd"]),
    P_re=as.numeric(model$pfas_params[[chem]]["P_re"]),
    mix_ratio=as.numeric(sub$mix_ratio),
    f_V_pl=f$f_V_pl,f_V_li=f$f_V_li,f_V_kd=f$f_V_kd,f_V_re=f$f_V_re,
    f_Q_li=f$f_Q_li,f_Q_kd=f$f_Q_kd,f_Q_re=f$f_Q_re,
    f_KMilk=function(age_yr) f$f_KMilk(age_yr)*mult$Milk_Intake,
    f_Free=renal$f_Free,f_Tm=renal$f_Tm,f_Kt=renal$f_Kt,
    f_Qfil=renal$f_Qfil,f_kust=renal$f_kust,f_ku=renal$f_ku,
    f_V_fil=renal$f_Vfil,f_C_milk=milk)
  init<-model$get_infant_initial_burden(chem,c0,f$f_V_pl(0),f$f_V_li(0),f$f_V_kd(0),f$f_V_re(0),
                                        pars$P_li,pars$P_kd,pars$P_re)
  y0<-c(Al=init[["Al"]],Ak=init[["Ak"]],Ar=init[["Ar"]],Ap=init[["Ap"]],
        Afil=0,Aust=0,Aurine=0,Afeces=0,Ainput=0)
  o<-as.data.frame(deSolve::ode(y=y0,times=TIME_VECTOR,
                                func=model$pbpk_ode_unified_tmkt,parms=pars,
                                method="lsoda",rtol=model$ODE_RTOL,atol=model$ODE_ATOL))
  days<-o$time/24;yrs<-o$time/model$HOURS_PER_YEAR
  plasma<-o$Ap/pmax(f$f_V_pl(yrs),model$MIN_VOLUME_L)
  liver<-o$Al/pmax(f$f_V_li(yrs),model$MIN_VOLUME_L)
  kidney<-o$Ak/pmax(f$f_V_kd(yrs),model$MIN_VOLUME_L)
  burden<-o$Ap+o$Al+o$Ak+o$Ar
  mass<-burden+o$Afil+o$Aust+o$Aurine+o$Afeces
  expected<-sum(y0[1:8])+o$Ainput
  err<-max(abs(mass-expected)/pmax(abs(expected),1e-18))
  if(!is.finite(err)||err>1e-5||any(!is.finite(c(plasma,liver,kidney)))||
     any(c(plasma,liver,kidney) < -1e-8)) stop("ODE QC failed: ",chem," ID ",sub$ID," err=",err)
  i_end<-which.min(abs(days-365.25))
  if(abs(days[i_end]-365.25)>0.001) stop("Missing exact endpoint")
  i_180<-which.min(abs(days-180))
  ep<-c(C_Plasma_12m=plasma[i_end],C_Plasma_6m=plasma[i_180],
        AUC_Plasma_0_12m=trapfast(days,plasma),
        C_Liver_12m=liver[i_end],C_Kidney_12m=kidney[i_end],
        Max_Mass_Error=err)
  if(!with_curve)return(list(endpoint=ep,curve=NULL))
  indices<-vapply(TARGET_DAYS,function(z)which.min(abs(days-z)),integer(1))
  cr<-data.frame(Target_Day=TARGET_DAYS,Actual_Day=days[indices],
                 C_Plasma=plasma[indices],C_Liver=liver[indices],C_Kidney=kidney[indices])
  list(endpoint=ep,curve=cr)
}

# --------- PARALLEL ENGINE (Windows PSOCK, persistent workers) ---------------
# One shared parameter draw -> 175 independent infant ODE solves, dispatched to
# persistent R workers. R/model state is isolated per worker. The full cohort
# is always used; the numerical ODE and all parameters are UNCHANGED.
# No worker writes the cache or Excel: writes happen on the master only.
cl10 <- NULL
if (WORKERS > 1L) {
  message("Starting ",WORKERS," persistent PSOCK workers; exporting physiology once...")
  cl10 <- parallel::makePSOCKcluster(WORKERS, outfile = "")
  worker_subject10 <- function(i, chem, mult, curve=TRUE) {
    sim_fast175(cohort[i, , drop=FALSE], chem, mult, with_curve=curve)
  }
  parallel::clusterExport(
    cl10,
    varlist=c("model","PHYS","TMKT0","PFAS0","TRANSFER0","TRAJ0",
              "TIME_VECTOR","TARGET_DAYS","TARGET_ACTUAL_DAY","cohort",
              "METRICS","trapfast","sim_fast175","worker_subject10"),
    envir=.GlobalEnv
  )
  parallel::clusterEvalQ(cl10, {
    if(!requireNamespace("deSolve",quietly=TRUE)) stop("deSolve missing on worker")
    if(!requireNamespace("dplyr",quietly=TRUE)) stop("dplyr missing on worker")
    NULL
  })
}
run_cohort10 <- function(chem, mult, curve=TRUE) {
  if(is.null(cl10)) {
    return(lapply(seq_len(nrow(cohort)),function(i)
      sim_fast175(cohort[i,,drop=FALSE],chem,mult,with_curve=curve)))
  }
  # parLapplyLB distributes uneven infant runtimes across workers.
  parallel::parLapplyLB(cl10,seq_len(nrow(cohort)),worker_subject10,
                        chem=chem,mult=mult,curve=curve)
}

# ---------------------------- BASELINE REPRO QC -----------------------------
# Test 2 complete cohort subjects per chemical before expensive GSA.
BASEQC<-list()
for(chem in CHEM){
  valid<-cohort[is.finite(cohort[[chem]]) & cohort[[chem]]>0,,drop=FALSE]
  if(nrow(valid)!=nrow(cohort)) stop("Full cohort has missing/nonpositive ",chem,
                                     "; do not silently vary analysis N; inspect the input before GSA")
  idx<-unique(c(1L,nrow(valid)))
  for(i in idx){
    sub<-valid[i,,drop=FALSE]
    sim<-sim_fast175(sub,chem,unit_multipliers())$endpoint
    ref<-base_endpoint[base_endpoint$ID==sub$ID & base_endpoint$Compound==chem,,drop=FALSE]
    if(nrow(ref)!=1L)stop("Missing baseline reference for ",chem," ID=",sub$ID)
    for(m in c("C_Plasma","C_Liver","C_Kidney")){
      pred<-unname(sim[paste0(m,"_12m")]);obs<-as.numeric(ref[[m]])
      difference<-abs(pred-obs)/max(abs(pred),abs(obs),1e-8)
      BASEQC[[length(BASEQC)+1L]]<-data.frame(Compound=chem,ID=sub$ID,
                                              Metric=m,Recalculated=pred,Frozen=obs,Relative_Difference=difference)
    }
  }
}
BASEQC<-do.call(rbind,BASEQC)
utils::write.csv(BASEQC,file.path(OUT,"FAST175_baseline_QC.csv"),row.names=FALSE)
if(any(!is.finite(BASEQC$Relative_Difference))||any(BASEQC$Relative_Difference>BASE_QC_REL_TOL))
  stop("FAST175 baseline QC FAILED. Check 01 version and input. See FAST175_baseline_QC.csv")
message("FAST175 baseline QC OK. Max relative diff: ",signif(max(BASEQC$Relative_Difference),3))

# --------------------------- RANK-BASED PRCC -------------------------------
prcc_fast <- function(X,y){
  X<-as.matrix(X);y<-as.numeric(y)
  if(any(!is.finite(X))||any(!is.finite(y)))stop("Nonfinite PRCC inputs")
  n<-nrow(X);k<-ncol(X)
  if(n<=k+8)stop("Not enough parameter draws for PRCC")
  if(!is.matrix(X)||k<1)stop("PRCC needs input matrix")
  ranks<-apply(X,2,rank,ties.method="average")
  if(k==1) ranks<-matrix(ranks,ncol=1)
  yr<-rank(y,ties.method="average")
  out<-vector("list",k)
  for(j in seq_len(k)){
    rest<-setdiff(seq_len(k),j)
    if(length(rest)){
      z<-cbind(Intercept=1,ranks[,rest,drop=FALSE])
      rx<-stats::lm.fit(z,ranks[,j])$residuals
      ry<-stats::lm.fit(z,yr)$residuals
    } else {rx<-ranks[,j]-mean(ranks[,j]);ry<-yr-mean(yr)}
    r<-if(stats::sd(rx)>0&&stats::sd(ry)>0)stats::cor(rx,ry) else NA_real_
    df<-n-k-1L
    tv<-if(is.finite(r)&&abs(r)<1)r*sqrt(df/max(1e-15,1-r*r)) else NA_real_
    out[[j]]<-data.frame(Parameter=colnames(X)[j],PRCC=r,
                         Approx_P=if(is.finite(tv))2*stats::pt(-abs(tv),df) else NA_real_)
  }
  do.call(rbind,out)
}
active_params <- function(chem){
  if(isTRUE(TRAJ0$supported_by_blomberg[match(chem,TRAJ0$Compound)]))
    setdiff(PARAM_NAMES,"Beta_Proxy") else PARAM_NAMES
}

# Convergence checks for each of 5 endpoints, NOT just plasma at 12 months.
# Relative stability of 2.5% and 97.5% across stages is checked separately.
get_stage_stats <- function(chem,dat,n){
  d<-dat[seq_len(n),,drop=FALSE]
  ps<-active_params(chem)
  X<-DESIGN[seq_len(n),ps,drop=FALSE]
  rows<-list(); quants<-list()
  for(m in METRICS){
    p<-prcc_fast(X,d[[m]])
    p$Compound<-chem;p$Metric<-m;p$N_Draw<-n
    rows[[m]]<-p
    quants[[m]]<-data.frame(Compound=chem,Metric=m,N_Draw=n,
                            P2_5=qfast(d[[m]],.025),Median=qfast(d[[m]],.5),P97_5=qfast(d[[m]],.975))
  }
  list(prcc=do.call(rbind,rows),quantiles=do.call(rbind,quants))
}
compare_stages <- function(a,b,prcc_tol,q_tol,top_k){
  pp<-merge(a$prcc,b$prcc,by=c("Compound","Metric","Parameter"),suffixes=c("_old","_new"))
  qq<-merge(a$quantiles,b$quantiles,by=c("Compound","Metric"),suffixes=c("_old","_new"))
  detail<-list()
  for(metric in unique(pp$Metric)){
    p<-pp[pp$Metric==metric,,drop=FALSE]
    p<-p[is.finite(p$PRCC_old)&is.finite(p$PRCC_new),,drop=FALSE]
    rank_old<-head(p$Parameter[order(-abs(p$PRCC_old))],top_k)
    rank_new<-head(p$Parameter[order(-abs(p$PRCC_new))],top_k)
    key<-union(rank_old,rank_new)
    pkey<-p[p$Parameter %in% key,,drop=FALSE]
    top_overlap<-setequal(rank_old,rank_new)
    stable_signs<-all(sign(pkey$PRCC_old)==sign(pkey$PRCC_new))
    max_prcc_change<-if(nrow(pkey))max(abs(pkey$PRCC_new-pkey$PRCC_old)) else Inf
    q<-qq[qq$Metric==metric,,drop=FALSE]
    q_relative<-vapply(c("P2_5","Median","P97_5"),function(col){
      old<-q[[paste0(col,"_old")]];new<-q[[paste0(col,"_new")]]
      max(abs(new-old)/pmax(abs(old),abs(new),1e-12))
    },numeric(1))
    detail[[metric]]<-data.frame(Metric=metric,
                                 Top_K_Same=top_overlap,Signs_Same=stable_signs,
                                 Max_Top_PRCC_Delta=max_prcc_change,
                                 Max_Quantile_Rel_Delta=max(q_relative),
                                 PRCC_Stable=top_overlap && stable_signs && max_prcc_change<=prcc_tol,
                                 Quantile_Stable=all(q_relative<=q_tol))
  }
  do.call(rbind,detail)
}

# -------------------- GSA: 175 INFANTS PER PARAMETER DRAW --------------------
STAGE_RESULT<-list();FINAL_BY_CHEM<-list();PROGRESS<-list();ALL_ENDPOINTS<-list()
for(chem in CHEM){
  message("FAST175 global sensitivity: ",chem," | n=",nrow(cohort)," | max draws=",MAX_DRAWS)
  # Read/write exactly one cache per chem; completed entries reused if fingerprint matches.
  cp<-file.path(CACHE_DIR,paste0("FAST175_",chem,".rds"))
  cache<-NULL
  if(file.exists(cp)&&!FORCE){
    old<-tryCatch(readRDS(cp),error=function(e)NULL)
    if(!is.null(old)&&identical(old$compound,chem)){
      if(identical(old$fingerprint,FINGERPRINT)) {
        cache<-old
      } else if(is.list(old$fingerprint)) {
        # Only allow reuse when numerical model/input/design/cohort/ODE tolerances
        # match exactly. Script hash differs solely because of parallel scheduling.
        f_old<-old$fingerprint; f_new<-FINGERPRINT
        prior_hash<-as.character(f_old$script_md5)
        f_old$script_md5<-NULL; f_new$script_md5<-NULL
        trusted_prior<-prior_hash %in% c(
          "b8103458602417377eecf0a638dfba23", # issued FAST175_FIXED ZIP
          "28a6133cbc5afbd1584bb7e0f89427cb"  # prior FAST175 source
        )
        if(identical(f_old,f_new) && (trusted_prior || ADOPT_CACHE)) {
          if(!trusted_prior) warning("User-authorized cache migration: unknown previous script hash. Ensure ODE logic was unchanged.")
          old$fingerprint<-FINGERPRINT
          cache<-old
          atomic_rds(cache,cp)
          message("Reusing verified original cache for ",chem,": ",cache$completed," draws")
        } else {
          stop("Found cache for ",chem," with incompatible fingerprint; will NOT overwrite previous work. ",
               "Check frozen model/input/LHS design, or back up cache before starting fresh.")
        }
      }
    }
  }
  if(is.null(cache))cache<-list(fingerprint=FINGERPRINT,compound=chem,completed=0L,
                                endpoint=list(),curve=list(),stages=list())
  # If existing completed stage meets convergence rules, reuse stopping decision.
  last_stage<-0L;stable_hits<-0L;finished<-FALSE
  previous_stats<-NULL
  for(s in seq_along(STAGES)){
    n<-STAGES[s]
    if(cache$completed<n){
      for(draw in seq.int(cache$completed+1L,n)){
        m<-as.list(DESIGN[draw,PARAM_NAMES,drop=FALSE])
        # For ALL 175 infants, same global parameter vector.
        endpoint_mat<-matrix(NA_real_,nrow(cohort),length(METRICS)+1L,
                             dimnames=list(NULL,c(METRICS,"Max_Mass_Error")))
        # No full 175 x 54 x (many columns) R data frame is retained per draw.
        all_curves<-array(NA_real_,dim=c(length(TARGET_DAYS),3L,nrow(cohort)))
        sims <- run_cohort10(chem,m,curve=TRUE)
        for(i in seq_along(sims)){
          z<-sims[[i]]
          endpoint_mat[i,]<-z$endpoint[colnames(endpoint_mat)]
          all_curves[,,i]<-as.matrix(z$curve[,c("C_Plasma","C_Liver","C_Kidney")])
        }
        rm(sims)
        medians<-apply(endpoint_mat[,METRICS,drop=FALSE],2,median)
        # Population P95 reflects variability between 175 infants at fixed parameters.
        high<-apply(endpoint_mat[,METRICS,drop=FALSE],2,qfast,p=.95)
        ep<-as.data.frame(as.list(c(Compound=chem,Draw=as.character(draw),
                                    N_Infants=as.character(nrow(cohort)),medians)))
        # enforce numerics (data.frame constructed from all-character c())
        ep$Draw<-as.integer(ep$Draw);ep$N_Infants<-as.integer(ep$N_Infants)
        for(met in METRICS){ep[[met]]<-as.numeric(ep[[met]])
        ep[[paste0(met,"_Population_P95")]]<-as.numeric(high[[met]])}
        ep$Max_Mass_Error<-max(endpoint_mat[,"Max_Mass_Error"])
        curve<-data.frame(Compound=chem,Draw=draw,Target_Day=TARGET_DAYS,
                          Actual_Day=TARGET_ACTUAL_DAY,
                          C_Plasma=apply(all_curves[,1,,drop=FALSE],1,median),
                          C_Liver=apply(all_curves[,2,,drop=FALSE],1,median),
                          C_Kidney=apply(all_curves[,3,,drop=FALSE],1,median))
        cache$endpoint[[draw]]<-ep
        cache$curve[[draw]]<-curve
        cache$completed<-draw
        if(draw%%CACHE_EVERY==0L||draw==n){
          atomic_rds(cache,cp)
          message(sprintf("  %s draw %d/%d | %d infants | %s",
                          chem,draw,MAX_DRAWS,nrow(cohort),format(Sys.time(),"%H:%M:%S")))
        }
      }
    }
    data_now<-do.call(rbind,cache$endpoint[seq_len(n)])
    stats_now<-get_stage_stats(chem,data_now,n)
    conv<-if(!is.null(previous_stats)) compare_stages(previous_stats,stats_now,
                                                      STABLE_PRCC_TOL,STABLE_Q_TOL,TOP_K) else data.frame(
                                                        Metric=METRICS,Top_K_Same=NA,Signs_Same=NA,
                                                        Max_Top_PRCC_Delta=NA_real_,Max_Quantile_Rel_Delta=NA_real_,
                                                        PRCC_Stable=FALSE,Quantile_Stable=FALSE)
    eligible<-n>=MIN_CONVERGENCE_N
    passed<-eligible && all(conv$PRCC_Stable) && all(conv$Quantile_Stable)
    stable_hits<-if(passed)stable_hits+1L else 0L
    conv$Compound<-chem;conv$N_Draw<-n;conv$Stage_Eligible<-eligible
    conv$Consecutive_Stable<-stable_hits;conv$All_Endpoints_Passed<-passed
    STAGE_RESULT[[length(STAGE_RESULT)+1L]]<-conv
    previous_stats<-stats_now
    last_stage<-n
    # Persist stage checks; stops automatically only when BOTH rank and extreme
    # conditional-parameter quantile stability hold for consecutive stages.
    cache$stages[[as.character(n)]]<-conv
    atomic_rds(cache,cp)
    if(stable_hits>=N_STABLE_CHECKS){finished<-TRUE;break}
  }
  final<-do.call(rbind,cache$endpoint[seq_len(last_stage)])
  FINAL_BY_CHEM[[chem]]<-list(stats=previous_stats,data=final,
                              curves=do.call(rbind,cache$curve[seq_len(last_stage)]))
  ALL_ENDPOINTS[[chem]]<-final
  PROGRESS[[chem]]<-data.frame(Compound=chem,Completed_Draws=last_stage,
                               Max_Draws=MAX_DRAWS,N_Infants=nrow(cohort),
                               PRCC_And_Quantiles_Stable=finished,
                               Stopped_At_Cap=!finished)
  message("Completed ",chem,": ",last_stage," draws. All stability criteria met: ",finished)
}

# ------------------------ COMPILE PRCC AND UNCERTAINTY ------------------------
prcc_all<-do.call(rbind,lapply(FINAL_BY_CHEM,function(z)z$stats$prcc))
prcc_all$Abs_PRCC<-abs(prcc_all$PRCC)
prcc_all$BH_FDR<-NA_real_
for(ch in unique(prcc_all$Compound))for(m in METRICS){
  ix<-which(prcc_all$Compound==ch&prcc_all$Metric==m)
  prcc_all$BH_FDR[ix]<-p.adjust(prcc_all$Approx_P[ix],"BH")
}
unc_all<-do.call(rbind,lapply(FINAL_BY_CHEM,function(z)z$stats$quantiles))
unc_all$Interpretation<-"Sensitivity range under assumed joint screening distributions; not calibrated 95% CI"
status<-do.call(rbind,PROGRESS)
stability<-do.call(rbind,STAGE_RESULT)
ep_all<-do.call(rbind,ALL_ENDPOINTS)
selected_bands<-do.call(rbind,lapply(names(FINAL_BY_CHEM),function(chem){
  curves<-FINAL_BY_CHEM[[chem]]$curves
  rows<-list()
  for(day in TARGET_DAYS)for(met in c("C_Plasma","C_Liver","C_Kidney")){
    vals<-curves[curves$Target_Day==day,met]
    rows[[length(rows)+1L]]<-data.frame(Compound=chem,Target_Day=day,
                                        Actual_Day=TARGET_ACTUAL_DAY[which(TARGET_DAYS==day)],Metric=met,
                                        N_Draw=length(vals),P2_5=qfast(vals,.025),Median=qfast(vals,.5),P97_5=qfast(vals,.975))
  }
  do.call(rbind,rows)
}))

# ----------------------- PFUdA INTERACTION | ALL 175 -------------------------
interaction<-data.frame()
if(DO_INTERACTION&&"PFUdA" %in% CHEM){
  vals<-c(.5,.75,1,1.25,1.5)
  grid<-expand.grid(Tmc=vals,Kt=vals)
  icp<-file.path(CACHE_DIR,"PFUdA_Tmc_Kt_grid175.rds")
  signature<-list(fingerprint=FINGERPRINT,grid=grid)
  ci<-NULL
  if(file.exists(icp)&&!FORCE){
    old<-tryCatch(readRDS(icp),error=function(e)NULL)
    if(!is.null(old)&&identical(old$signature,signature))ci<-old
  }
  if(is.null(ci))ci<-list(signature=signature,rows=list())
  for(j in seq_len(nrow(grid))){
    if(length(ci$rows)>=j && !is.null(ci$rows[[j]]))next
    mult<-unit_multipliers();mult$Tmc<-grid$Tmc[j];mult$Kt<-grid$Kt[j]
    sims <- run_cohort10("PFUdA",mult,curve=FALSE)
    data <- vapply(sims,function(z){
      e<-z$endpoint
      e[c("C_Plasma_12m","C_Kidney_12m","AUC_Plasma_0_12m","Max_Mass_Error")]
    },numeric(4))
    rm(sims)
    ci$rows[[j]]<-data.frame(Compound="PFUdA",Tmc=grid$Tmc[j],Kt=grid$Kt[j],
                             N_Infants=nrow(cohort),Plasma_Median=median(data[1,]),
                             Kidney_Median=median(data[2,]),AUC_Median=median(data[3,]),
                             Max_Mass_Error=max(data[4,]))
    atomic_rds(ci,icp)
    message("PFUdA renal grid: ",j,"/",nrow(grid))
  }
  interaction<-do.call(rbind,ci$rows)
}

if(!is.null(cl10)){
  parallel::stopCluster(cl10)
  cl10 <- NULL
  message("Parallel workers stopped cleanly.")
}

# --------------------------- PUBLICATION OUTPUT ------------------------------
notes<-data.frame(
  Item=c("Population","Method","LHS blocks","Adaptation","PRCC convergence",
         "Uncertainty convergence","Sampling distributions","Conditional uncertainty",
         "Population variation","No pseudoreplication","Frozen renal model",
         "Input slope assumptions","Parametric dependence","Ku missing","Time-point matching","Status"),
  Interpretation=c(
    paste0("Full cohort (",nrow(cohort)," matched infants) per parameter draw"),
    "Common parameter vectors applied to all infants; PRCC computed on per-draw population median",
    "Independent stratified LHS blocks; combined stages are batched LHS, not exact nested LHS",
    paste(STAGES,collapse=" -> "),
    paste0("Top ",TOP_K," same, sign stable, max absolute top PRCC difference <=",STABLE_PRCC_TOL,
           "; ",N_STABLE_CHECKS," consecutive eligible comparisons"),
    paste0("For every output: P2.5, median, P97.5 relative change <=",STABLE_Q_TOL,
           "; ",N_STABLE_CHECKS," consecutive eligible comparisons"),
    "Hypothetical bounded independent parameter screening distributions, not literature-based probabilities",
    "Across draws: parameter-screening sensitivity of cohort median, NOT calibrated 95% prediction CI",
    "Each draw also stores across-infant 95th percentile; differs from cross-draw percentiles",
    "N_PRCC equals parameter draws, NOT 175 times number of draws",
    "01 model ODE functions extracted; baseline and mass-balance QC fail closed",
    "PFNA/PFDA/PFUdA beta multiplicative 0 to 1.5; published-slope compounds beta fixed",
    "Independent sampled parameters are screening approximation; physiological/kinetic correlations not established",
    "PFDA/PFUdA Ku remain NA; no arbitrary empirical imputation",
    "Selected days use nearest actual weekly output; exact matched Actual_Day is exported",
    "Execution is complete only when user successfully runs R script; check status and QC"
  ),stringsAsFactors=FALSE)
utils::write.csv(prcc_all,file.path(OUT,"FAST175_PRCC.csv"),row.names=FALSE)
utils::write.csv(unc_all,file.path(OUT,"FAST175_parameter_screening_bounds.csv"),row.names=FALSE)
utils::write.csv(stability,file.path(OUT,"FAST175_convergence.csv"),row.names=FALSE)
utils::write.csv(status,file.path(OUT,"FAST175_run_status.csv"),row.names=FALSE)
utils::write.csv(ep_all,file.path(OUT,"FAST175_population_metrics_by_draw.csv"),row.names=FALSE)
if(nrow(interaction))utils::write.csv(interaction,file.path(OUT,"FAST175_PFUdA_interaction.csv"),row.names=FALSE)

pr_plot<-prcc_all[prcc_all$Metric=="C_Plasma_12m",,drop=FALSE]
pr_plot$Parameter<-factor(pr_plot$Parameter,levels=rev(PARAM_NAMES))
p<-ggplot2::ggplot(pr_plot,ggplot2::aes(PRCC,Parameter,fill=PRCC>=0))+
  ggplot2::geom_col(width=.75)+ggplot2::facet_wrap(~Compound,ncol=3)+
  ggplot2::geom_vline(xintercept=0,linewidth=.25)+
  ggplot2::coord_cartesian(xlim=c(-1,1))+
  ggplot2::scale_fill_manual(values=c("TRUE"="#BA5B4B","FALSE"="#4579A4"),guide="none")+
  ggplot2::theme_bw(base_size=11)+
  ggplot2::labs(x="Partial rank correlation coefficient",y=NULL,
                title="FAST175: sensitivity of 12-month plasma concentration")
ggplot2::ggsave(file.path(OUT,"FAST175_PRCC_plasma.png"),p,width=12,height=8,dpi=300)
ggplot2::ggsave(file.path(OUT,"FAST175_PRCC_plasma.pdf"),p,width=12,height=8)
if(nrow(interaction)){
  g<-ggplot2::ggplot(interaction,ggplot2::aes(Tmc,Kt,fill=Plasma_Median))+
    ggplot2::geom_tile()+ggplot2::geom_contour(
      ggplot2::aes(z=Plasma_Median),colour="white",linewidth=.35)+
    ggplot2::theme_bw(base_size=11)+
    ggplot2::labs(title="PFUdA: full-cohort renal parameter interaction",
                  x="Tmc multiplier",y="Kt multiplier",fill="Plasma ng/L")
  ggplot2::ggsave(file.path(OUT,"FAST175_PFUdA_Tmc_Kt.png"),g,width=7,height=5,dpi=300)
  ggplot2::ggsave(file.path(OUT,"FAST175_PFUdA_Tmc_Kt.pdf"),g,width=7,height=5)
}
wb<-list(Run_Status=status,Baseline_QC=BASEQC,Parameter_Ranges=PARAMS,
         PRCC_All=prcc_all,Endpoint_Ranges=unc_all,Selected_Day_Bands=selected_bands,
         Convergence=stability,PFUdA_Tmc_Kt=interaction,Interpretation=notes)
wb<-wb[vapply(wb,function(d)is.data.frame(d)&&nrow(d)>0,logical(1))]
stopifnot(all(nchar(names(wb))<=31))
excel<-file.path(OUT,"global_sensitivity_FAST175.xlsx")
writexl::write_xlsx(wb,excel)
message("FAST175 completed: ",excel)
if(any(!status$PRCC_And_Quantiles_Stable))
  warning("Some chemicals reached max draws without meeting ALL convergence criteria. See Convergence sheet; do not claim convergence.")
