# MASLD CMRF research code

Maintainer and research-code author: mjwang.

This repository provides the statistical analysis code for a multicentre cross-sectional study of cardiometabolic risk-factor (CMRF) count, components and glycaemic categories in relation to liver histology. It includes the frozen complete-component analysis, multiple-imputation sensitivity analyses, exploratory comparisons and plotting source used during manuscript preparation.

## Data and scope

No individual-level data, imputed datasets, source workbooks, manuscript files or credentials are included. Reproducing study results requires authorised access to the original source workbook. The header-only templates describe the interface; they contain no observations and are not a demonstration cohort. The scripts assert the original study dimensions and are not a general clinical calculator.

The executed definitions are preserved: Asian adiposity thresholds, the combined glycaemic criterion, recorded hypertension, measured triglycerides and sex-specific strict-less-than HDL-C thresholds. The source extract does not contain separate medication or blood-pressure fields. The code must not be interpreted as having reconstructed unavailable branches of the consensus definition. ZS source prefixes are assigned to RJ in the derived centre variable; they are not analysed as a fourth centre.

## Analysis sequence

Run `python3 run_analysis.py --help` for the interface. The default mode prints the plan without copying data or running models. Use an empty work-directory path outside this repository and an authorised workbook to execute:

```sh
python3 run_analysis.py --input-workbook /path/to/authorised_workbook.xlsx --work-dir /path/to/new_private_run --execute
```

The runner copies the restricted workbook into the private work directory without modifying the source. Generated participant-level files and imputations remain there; do not upload that directory.

1. `20260809_152651/01_reanalysis.R`: adult eligibility, conservative five-domain reconstruction, complete-component models, fibrosis probabilities and raw-variable multiple imputation (50 datasets, 20 iterations).
2. `20260809_152651/03_qc_exploratory.R`: diagnostics, within-count contrasts, BMI interactions and imputation diagnostics.
3. `20261003_165814/04_correct_mi_pooling.R`: final MI pooling using the fitted model rank for the complete-data degrees of freedom.
4. `20261018_173329/center_comparison.R`: centre-specific descriptive comparisons.
5. `20261023_001156/01_steatosis_components.R`: complete-component and imputed steatosis component models.
6. `20261016_225206/01_pattern_vs_rest.R`: 11 common profiles versus all other participants; unadjusted, exploratory comparisons.
7. `20261027_230711/01_steatosis_probabilities.R`: marginal-standardized S≥2 and S3 probabilities with centre-stratified bootstrap intervals.

The original R environment was 4.5.3. The recorded core versions are listed in `requirements-analysis.csv`; no dependencies are installed automatically. The validation command reports the locally available versions. Seeds and model formulas are retained in the scripts. The historical plotting scripts additionally use grid, ggplot2, patchwork, dplyr, tidyr, svglite, ragg, ggsci and magick.

## Plotting source

The `source_scripts` tree preserves the original project-relative paths. Historical plotting scripts read version-specific aggregate source-data tables and sometimes previous exports for parity checks. They are included for transparency; the runner executes the seven statistical stages above, not manuscript/PDF assembly or every historical rendering step. The source summaries and prior exports needed by those plotting scripts are not bundled. A byte-identical figure regeneration is therefore not asserted for this code-only release.

## Reproducibility notes

Reproduction requires authorised source data and the stated software environment. The study models were not rerun while preparing this repository, and the historical figure dependencies described above are not bundled. In the frozen MI preprocessing, a nonbinary HCV value is converted to missing by binary-factor coding. `source_manifest.json` records source hashes and the limited packaging changes: mjwang author headers and an optional historical display-parity check; the statistical implementation is otherwise retained.
