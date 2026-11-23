# Model and analysis version audit

## Blocking discrepancy

The supplied main-model script `01_run_main_model(1).R` begins by describing a revised S1 candidate/v2 model and writes to the result version directory `PFAS_infant_PBPK_S1_renal_totalplasma_v2`. Its changelog describes the corrected pediatric GFR equation, revised transfer parameters, urinary transit/storage structure, mass-balance bookkeeping, and growth individualization.

The supplied current Supplementary Methods, however, names the final model source as:

`PFAS_infant_PBPK_final_freeze_v1_1_endpoint365p25.R`

These identifiers are not interchangeable. The code file supplied for this repository is renamed to `01_run_main_model.R` only for repository organization; this does not prove it is byte-for-byte or numerically identical to the source used for the submitted manuscript.

## Required author-side confirmation

Before describing the repository as an exact reproduction package, verify:

1. Which exact R source file generated the PBPK workbook used for the final manuscript.
2. Whether that workbook contains 175 infants × 6 compounds × 54 time points, and whether its checksum/date matches the final analytical workbook.
3. Whether the manuscript's reported EDI values were calculated from the same `Ingestion_Rate_ng_h`, body-weight trajectory, milk trajectory and breastfeeding mixture ratio as the supplied script.
4. Whether the reported CTWA values use AUC over days 0–365 divided by 365, with linear interpolation at day 365.
5. Whether the compound-specific serum-to-plasma ratios in the manuscript are exactly those implemented in the NO-WATER risk and figure scripts.

## Other version conflicts resolved in this repository snapshot

- The earlier R15 scripts include water analysis and/or generic serum-to-plasma scenarios. They are not used in the primary manuscript pipeline.
- The manuscript's internal benchmark uses compound-specific lower/central/upper serum-to-plasma ratios (PFOA/PFNA/PFHxS/PFOS: 1.1/1.2/1.0/1.1; 1.2/1.3/1.1/1.1; 1.2/1.4/1.1/1.2). The new NO-WATER script implements those values.
- The manuscript does not compare infant milk-derived EDI directly with EFSA's maternal group TWI of 4.4 ng/kg bw/week.
- The manuscript does not calculate liver/kidney HQ or HI.
- Figure 3 is the within-organ six-PFAS CTWA composition; Figure 4 is the individual four-PFAS RBMDL10 distribution.

## Release rule

Do not remove this audit note until the exact model source and output workbook used for the manuscript have been confirmed by the authors.
