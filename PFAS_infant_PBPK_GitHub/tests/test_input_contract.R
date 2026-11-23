# Basic input contract checks for the S1 model result workbook.
# Run after the main model: Rscript tests/test_input_contract.R
root <- normalizePath(Sys.getenv("PFAS_PBPK_ROOT", getwd()), winslash = "/", mustWork = TRUE)
version <- "PFAS_infant_PBPK_S1_renal_totalplasma_v2"
p <- Sys.getenv("PFAS_PBPK_RESULTS", unset = file.path(
  root, "results", version, "tables",
  "PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx"
))
stopifnot(file.exists(p))
stopifnot(requireNamespace("readxl", quietly = TRUE))
x <- readxl::read_excel(p, sheet = "PBPK_Results")
required <- c("ID","Compound","Scenario","Time_days","BW_current_kg",
              "Ingestion_Rate_ng_h","C_Plasma","C_Liver","C_Kidney",
              "Mass_Balance_Relative_Error")
missing <- setdiff(required, names(x))
if (length(missing)) stop("Missing PBPK result columns: ", paste(missing, collapse = ", "))
x <- x[x$Scenario == "S1", ]
stopifnot(length(unique(x$ID)) == 175L)
stopifnot(setequal(unique(x$Compound), c("PFOA","PFNA","PFDA","PFUdA","PFHxS","PFOS")))
stopifnot(!anyDuplicated(x[c("ID","Compound","Time_days")]))
message("Input contract passed. This check does not establish scientific validation.")
