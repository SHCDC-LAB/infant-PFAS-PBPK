# Revised evidence evaluation for S1 kidney-corrected model (2026-11-03)
# No stale frozen Varsi predictions may be reused as current predictions.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(readxl);library(dplyr);library(ggplot2);library(tidyr);library(writexl)})
ROOT <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT",getwd()),mustWork=FALSE,winslash="/")
VERSION <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
MODE <- tolower(Sys.getenv("PFAS_TRAJECTORY_MODE","proxy"))
if(MODE != "proxy") VERSION <- paste0(VERSION,"_beta_",MODE)
OUT <- file.path(ROOT,"results",VERSION)
TBL <- file.path(OUT,"tables")
dir.create(TBL,recursive=TRUE,showWarnings=FALSE)
PATH <- file.path(TBL,"PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx")
if(!file.exists(PATH)) stop("Run R/01_run_main_model.R first. Missing: ",PATH)
d <- readxl::read_excel(PATH,sheet="PBPK_Results") %>% filter(Scenario=="S1")
need <- c("ID","Compound","Time_days","C_Plasma","C_Liver","C_Kidney","BW_current_kg","R_Urine_Excretion","R_Filtration","R_Reabsorption","Mass_Balance_Relative_Error")
missing <- setdiff(need,names(d));if(length(missing)) stop("Missing columns: ",paste(missing,collapse=", "))
quant <- function(x,p) {v<-x[is.finite(x)];if(!length(v))return(NA_real_);as.numeric(quantile(v,p))}
median_ok <- function(x) {v<-x[is.finite(x)];if(!length(v))return(NA_real_);median(v)}
# Numerical QC / full individual-level predicted plasma and organ concentrations.
qc <- d %>% group_by(Compound) %>% summarise(N_ID=n_distinct(ID),N_Rows=n(),
  Max_Abs_MassBalance=if(all(is.na(Mass_Balance_Relative_Error))) NA_real_ else max(abs(Mass_Balance_Relative_Error),na.rm=TRUE),
  Min_Plasma=min(C_Plasma,na.rm=TRUE),Min_Liver=min(C_Liver,na.rm=TRUE),Min_Kidney=min(C_Kidney,na.rm=TRUE),.groups="drop")
organs <- d %>% select(ID,Compound,Time_days,C_Plasma,C_Liver,C_Kidney) %>%
 pivot_longer(cols=c(C_Plasma,C_Liver,C_Kidney),names_to="Compartment",values_to="Concentration_ng_L") %>%
 group_by(Compound,Time_days,Compartment) %>% summarise(N=sum(is.finite(Concentration_ng_L)),Median=median_ok(Concentration_ng_L),
 P025=quant(Concentration_ng_L,.025),P975=quant(Concentration_ng_L,.975),.groups="drop")
# Paired ratios (not ratio of two population medians).
paired <- d %>% filter(C_Plasma>0) %>% mutate(Liver_Plasma=C_Liver/C_Plasma,Kidney_Plasma=C_Kidney/C_Plasma) %>%
 group_by(Compound,Time_days) %>% summarise(N=n_distinct(ID),Liver_Median=median_ok(Liver_Plasma),Kidney_Median=median_ok(Kidney_Plasma),
 Liver_P025=quant(Liver_Plasma,.025),Liver_P975=quant(Liver_Plasma,.975),Kidney_P025=quant(Kidney_Plasma,.025),Kidney_P975=quant(Kidney_Plasma,.975),.groups="drop")
# Neonatal external benchmark, cited study-level observations, NOT same-subject pairing.
yao <- data.frame(Compound=c("PFOA","PFNA","PFDA","PFHxS"),Observed_CLrenal_mL_kg_day=c(.049,.047,.022,.270),
 Source="Yao 2023 neonatal cohort summary; historical benchmark; verify sampling window")
cl <- d %>% filter(Time_days>0) %>% mutate(
 CL_Urine=ifelse(C_Plasma>0 & BW_current_kg>0,R_Urine_Excretion/C_Plasma*1000*24/BW_current_kg,NA_real_),
 CL_Net=ifelse(C_Plasma>0 & BW_current_kg>0,pmax(R_Filtration-R_Reabsorption,0)/C_Plasma*1000*24/BW_current_kg,NA_real_)) %>%
 group_by(ID,Compound) %>% slice_min(abs(Time_days-7),with_ties=FALSE) %>% ungroup() %>%
 group_by(Compound) %>% summarise(Day=median(Time_days),N=n_distinct(ID),Pred_Urine=median_ok(CL_Urine),Pred_Net=median_ok(CL_Net),.groups="drop") %>%
 left_join(yao,by="Compound") %>% mutate(PO_Urine=Pred_Urine/Observed_CLrenal_mL_kg_day,PO_Net=Pred_Net/Observed_CLrenal_mL_kg_day,
 Evidence="Neonatal renal clearance external benchmark; urine and net renal entry separately")
# Adult half-lives from old benchmark are not shown as CURRENT predictions; modified ODE invalidates them.
adult_old <- data.frame(Compound=c("PFOA","PFNA","PFDA","PFUdA","PFHxS","PFOS"),
 Old_Frozen_Adult_Predicted_days=c(1392.9,1670.6,1122.6,568.2,3756.5,2412.8),
 Adult_Observed_days=c(2348.1,1581.3,1122.6,568.2,1885.2,1420.8),
 Current_Status="ARCHIVE ONLY: rerun adult model under corrected kidney exchange for current qualification")
# Study-specific Varsi external evaluation MUST be regenerated from corrected physiology.
# If independently produced revised workbook exists, load it; otherwise do not claim any new plasma agreement.
vpath <- file.path(OUT,"Varsi_revised","tables","Varsi2022_final_freeze_validation_summary_v1.xlsx")
varsi_status <- data.frame(Comparison="Study-specific Varsi 6-month infant plasma",Status=if(file.exists(vpath))"Revised workbook present; see R/07 output, still verify trajectory/input provenance" else
 "NOT RECOMPUTED for corrected model; historical Varsi predictions excluded",Current_File=vpath)
varsi_current <- data.frame(Status="Not rerun", Note="Revised Varsi study-specific simulation required")
if(file.exists(vpath)) {
  varsi_current <- tryCatch(readxl::read_excel(vpath,sheet="Final_Primary_Validation"),error=function(e){warning("Varsi workbook unreadable: ",conditionMessage(e));data.frame()})
  if(nrow(varsi_current)) {
    if(!"Model_Version" %in% names(varsi_current) || any(varsi_current$Model_Version != VERSION)) {
      warning("Varsi workbook model version mismatch: DO NOT interpret as current model validation")
      varsi_status$Status <- "VERSION MISMATCH: results not eligible for current-model validation"
      varsi_current <- data.frame(Status="Not rerun", Note="Revised Varsi study-specific simulation required")
    }
  }
}

# Longitudinal benchmark: different age horizon, direction only.
wu <- data.frame(Compound=c("PFOA","PFNA","PFDA","PFUdA","PFHxS","PFOS"),Wu_Annual_Change_pct=c(9.6,11.4,10.8,6.7,31.1,9.6))
firstlast <- d %>% group_by(ID,Compound) %>% arrange(Time_days,.by_group=TRUE) %>% summarise(C0=first(C_Plasma),C12=last(C_Plasma),.groups="drop") %>%
 mutate(Change_pct=ifelse(C0>0,100*(C12/C0-1),NA_real_)) %>% group_by(Compound) %>% summarise(Model_Median_Change_pct=median_ok(Change_pct),.groups="drop") %>%
 left_join(wu,by="Compound") %>% mutate(Level="Direction-only plausibility; 0-12 months vs 0-7 years is not direct absolute validation")
notes <- data.frame(Check=c("Scope","Varsi","Organ","Neonatal CL","Adult","Longitudinal"),
Meaning=c("S1 primary; PFNA/PFDA/PFUdA milk trends are PFOA-homolog proxies, not observations",
 "Do not reuse frozen Varsi prediction numbers; revised study-specific simulation required",
 "Actual individual predictions exported; independent paired infant liver/kidney observations are not available",
 "Predicted urine outflow is not identical to net-filtrate entry; missing Ku can make output zero",
 "Historic adult predictions are archived, not valid under new exchange equation",
 "Wu is beyond infant model age range; direction only"))
writexl::write_xlsx(list(Readme=notes,Numerical_QC=qc,Organ_Timecourse=organs,Paired_Partition_Ratios=paired,
 Yao_Day7=cl,Varsi_Status=varsi_status,Varsi_Recomputed_If_Available=varsi_current,Adult_Archive_NotCurrent=adult_old,Wu_Direction=firstlast),file.path(TBL,"S1_external_and_numerical_evaluation.xlsx"))
if(any(qc$Max_Abs_MassBalance>1e-5,na.rm=TRUE)) warning("Mass balance QC exceeded threshold")
message("S1 evaluation complete. Varsi has NOT been claimed validated unless rerun separately.")
