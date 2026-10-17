# Infant PFAS PBPK model

This repository contains the R code used for the infant PFAS PBPK analysis of PFOA, PFNA, PFDA, PFUdA (PFUnDA), PFHxS, and PFOS from birth to 12 months.

## Analysis scripts

- `R/01_run_main_model.R` — primary six-PFAS PBPK simulation, S0/S1 exposure scenarios, QC, summaries, and plots.
- `R/02_four_layer_evaluation.R` — adult, neonatal renal, 6-month infant, and childhood longitudinal evaluation/qualification framework.
- `R/03_source_attribution_sensitivity.R` — prenatal/postnatal source-attribution and transfer-ratio sensitivity analysis.

## Primary model features

- age-dependent infant physiology and growth;
- measured colostrum PFAS concentrations as exposure anchors;
- S0: constant milk concentration;
- S1: literature-informed longitudinal milk trajectories;
- plasma, liver, kidney tissue, rest-of-body, glomerular filtrate, urinary transit/storage, cumulative urine, and fecal states;
- glomerular filtration and saturable Tm/Kt renal reabsorption;
- strict mass-balance checking;
- daily physiology/renal cache with adaptive LSODA integration;
- 7-day reporting plus the exact 365.25-day endpoint.

The 2.5th–97.5th percentile bands are between-infant prediction distributions, not confidence intervals.

## Required cohort input

Place the cohort workbook at:

`data/input-PBPK.xlsx`

The workbook must contain a worksheet named `Sheet1`.

Required columns:

`BW`, `infsex`, `PFOA`, `PFNA`, `PFDA`, `PFUdA`, `PFHxS`, `PFOS`

Optional:

`BFpc`

Frozen input conventions:

- birth weight (`BW`) is in grams;
- milk PFAS concentrations are in ng/L;
- `infsex`: 1 = male, 2 = female;
- `BFpc`, when present, is constrained to 0–1; missing/invalid values default to 1.

Censored values written as `<x` are replaced by `x/sqrt(2)`. Values written as `>x` are parsed as `x`.

## Running

From the repository root:

```r
source("R/01_run_main_model.R")
```

Optional follow-up analyses:

```r
source("R/02_four_layer_evaluation.R")
source("R/03_source_attribution_sensitivity.R")
```

You can override the project root with `PFAS_PBPK_ROOT` and the cohort workbook with `PFAS_PBPK_INPUT`.

## Required R packages

`readxl`, `dplyr`, `tidyr`, `ggplot2`, `deSolve`, `writexl`, `future`, `future.apply`

The scripts install missing packages automatically.

## External evidence

External literature files are intentionally not bundled unless redistribution is permitted. Optional files can be placed under `data/external/`.

The evaluation framework uses published evidence from Abraham et al. (2024), Yao et al. (2023), Varsi et al. (2022), and Wu et al. (2026). Complete parameter provenance and references should be reported in the manuscript and Supplementary Material.

## Interpretation notes

- PFDA and PFUdA effective Tmc values are current-study calibrated parameters; corresponding adult comparisons are calibration checks, not independent validation.
- PFDA/PFUdA downstream Ku is unresolved in the frozen model; plasma kinetics remain usable, but urinary-output interpretation is limited.
- The downstream urinary transit/storage state is an effective model compartment, not a validated anatomical bladder compartment.
- S1 is flat for PFNA, PFDA, and PFUdA because comparable longitudinal milk data were unavailable in the frozen analysis.
- Childhood comparisons extending beyond 12 months are plausibility checks rather than formal validation.

## Data privacy

Do not commit participant-level cohort data to a public GitHub repository unless the ethics approval and participant consent explicitly allow public release. The supplied `.gitignore` excludes the cohort workbook and external literature files by default.

## Reproducibility

Each analysis writes `sessionInfo.txt` to its output directory.
