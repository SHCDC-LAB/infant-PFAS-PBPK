# Infant PFAS PBPK model and screening-level risk characterization

This repository contains the R code selected for the current manuscript:

**Breast-milk PFAS exposure and PBPK-predicted internal concentrations in infants: a six-compound exposure and risk characterization in Shanghai**

The analysis covers 175 mother–infant pairs and six PFAS (PFOA, PFNA, PFDA, PFUdA/PFUnDA, PFHxS, PFOS) during the first postnatal year. The internal benchmark comparison uses the four-PFAS mixture PFOA + PFNA + PFHxS + PFOS.

## Scope

The public analysis workflow is deliberately **NO-WATER**. Drinking-water monitoring, water–milk regional comparisons, water–PBPK ecological linkage, and water ND sensitivity analyses are not part of the current manuscript and are not included in this repository snapshot.

The main analysis sequence is:

1. S1 infant PBPK simulation with the unified renal Tm/Kt and urinary transit/storage structure.
2. Milk-derived EDI reconstruction from the model's recorded ingestion rate and age-specific body weight.
3. Plasma, liver, and kidney concentration-time summaries and within-organ six-PFAS CTWA composition.
4. Day-365 serum-equivalent four-PFAS sum and exploratory comparison with the EFSA 1-year serum BMDL10 point of departure (17.5 ng/mL), using the manuscript's compound-specific serum-to-plasma ratio scenarios.
5. Manuscript figures 1–4 and Supplementary Figure S1.

No infant EDI is compared directly with the EFSA maternal group TWI. No organ-specific HQ/HI is calculated because matching tissue-specific internal-dose toxicity benchmarks are not available. The BMDL10 ratio is a screening-level benchmark comparison, not a clinical classification or probability of disease.

## Repository layout

```text
R/
  00_run_pipeline.R
  01_run_main_model.R
  02_four_layer_evaluation.R
  06_trajectory_proxy_sensitivity.R
  07_Varsi_study_specific_revalidation.R
  08_integrated_sensitivity_analysis.R
  10_global_sensitivity_uncertainty_FAST175.R
  10b_FAST175_posthoc_convergence_noODE.R
  15_infant_PFAS_risk_characterization_NO_WATER.R
  16_export_main_figures_1_2_NO_WATER.R
  17_export_IBQ_organ_contribution_figures_NO_WATER.R
  18_export_Supplementary_Figure_S1_NO_WATER.R
data/
  README.md
results/
  .gitkeep
docs/
  MANUSCRIPT_CODE_MAP.md
  VERSION_AUDIT.md
tests/
  test_input_contract.R
```

## Data privacy

The cohort workbook is **not included**. Place the authorized `input-PBPK.xlsx` file in the local `data/` directory. Do not commit participant-level data, identifiers, raw laboratory files, or non-public results to a public GitHub repository.

The analysis scripts also consume generated results under `results/`; this directory is ignored by Git except for `.gitkeep`.

## Requirements

Use R (recommended R 4.3 or later). Install packages used by the scripts:

```r
install.packages(c(
  "readxl", "writexl", "openxlsx", "dplyr", "tidyr", "tibble",
  "ggplot2", "scales", "stringr", "deSolve", "future",
  "future.apply", "minpack.lm", "purrr"
))
```

Some optional evaluation and sensitivity modules may require additional packages. Each script checks its own dependencies.

## Run the primary manuscript pipeline

Open the repository root as the working directory, then run:

```r
source("R/00_run_pipeline.R")
```

Or from a shell:

```bash
Rscript R/00_run_pipeline.R
```

The primary pipeline runs the S1 model and then the NO-WATER risk and figure scripts. The main PBPK run is computationally intensive.

### Optional supporting analyses

Set an environment variable to `1` before running the pipeline to activate a supporting module:

```r
Sys.setenv(PFAS_RUN_FOUR_LAYER = "1")
Sys.setenv(PFAS_RUN_VARSI = "1")
Sys.setenv(PFAS_RUN_TRAJECTORY_SENSITIVITY = "1")
Sys.setenv(PFAS_RUN_INTEGRATED_SENSITIVITY = "1")
Sys.setenv(PFAS_RUN_FAST175 = "1")
Sys.setenv(PFAS_RUN_FAST175_POSTHOC = "1")
source("R/00_run_pipeline.R")
```

FAST175 and integrated sensitivity analyses can require substantial computation. The post-hoc FAST175 diagnostic uses existing FAST175 CSV files and does not run PBPK ODEs.

## Important model-version note

Before using this repository to claim exact numerical reproduction of the submitted manuscript, resolve the version discrepancy documented in [`docs/VERSION_AUDIT.md`](docs/VERSION_AUDIT.md). The supplied main-model script identifies itself as a revised S1 candidate/v2 implementation, while the submitted Supplementary Methods names `PFAS_infant_PBPK_final_freeze_v1_1_endpoint365p25.R`. The exact source script and result workbook used to generate the manuscript tables must be confirmed. Do not assume that a filename rename establishes equivalence.

## Outputs

The main model writes its S1 time series to:

`results/PFAS_infant_PBPK_S1_renal_totalplasma_v2/tables/PBPK_all_S1_unified_TmKt_UrinaryStorage_results.xlsx`

The no-water risk script writes:

`results/PFAS_infant_PBPK_S1_renal_totalplasma_v2/risk_characterization_NO_WATER/infant_PFAS_risk_characterization_NO_WATER.xlsx`

Figure scripts write journal-ready PDF and 900-dpi TIFF files beneath `results/`. Inspect all QC sheets and compare generated summaries against the manuscript before submission.

## Reproducibility and interpretation

- The model's predicted liver and kidney concentrations are model outputs, not infant tissue measurements.
- The human organ literature is contextual evidence and does not constitute infant absolute-concentration validation.
- The six-PFAS within-organ percentages are concentration-composition metrics, not mass-based organ burden fractions or toxicity-weighted contributions.
- The serum-to-plasma ratios are literature-informed scenarios, not infant-specific measured conversion factors.
- The day-365 BMDL10 comparison is exploratory and does not estimate disease probability.
- The PFNA/PFDA/PFUdA milk trajectory uses the PFOA monthly decline as a structural-homolog proxy, as described in the manuscript; this assumption is not independently validated.

## Citation

Please cite the associated manuscript. A repository-specific citation and license should be added after the authors confirm the public release version and licensing terms.
