# MASLD CMRF research

Author and maintainer: mjwang.

This repository contains one analysis pipeline for a multicentre cross-sectional study of cardiometabolic risk-factor (CMRF) count, individual components and glycaemic categories in relation to liver histology. The numbered files are successive analysis steps, not competing versions of the study. Earlier code history is retained by Git rather than duplicated in the current directory tree.

## Layout

- `analysis/`: seven statistical steps in execution order.
- `figures/`: plotting scripts used for the manuscript figures and tables.
- `templates/`: header-only data interface and variable definitions, with no study observations.
- `run_analysis.py`: the single analysis entry point.
- `requirements-analysis.csv`: recorded R and package versions.
- `source_manifest.json`: provenance and hashes of the unchanged statistical scripts.

## Run the analysis

An authorised copy of the original workbook is required. The repository does not contain individual-level data, imputed datasets or manuscript files. The scripts assert the study's original dimensions and are not a general clinical risk calculator.

```sh
python3 run_analysis.py --check-environment
python3 run_analysis.py --work-dir /path/to/new_private_run
python3 run_analysis.py --input-workbook /path/to/authorised_workbook.xlsx --work-dir /path/to/new_private_run --execute
```

The default command prints the seven-step plan without writing data or running models. Execution requires a new private work directory outside the code repository. The runner copies the workbook into that directory, leaves the source unchanged, and reconstructs the internal paths expected by the original scripts. Participant-level outputs and imputations stay in that private workspace and must not be uploaded. Dependencies are not installed automatically.

## Analysis steps

1. Reconstruct eligibility and CMRFs; fit complete-component models and generate the 50 raw-variable imputations.
2. Run diagnostics and exploratory within-count comparisons.
3. Pool imputed models using the fitted model rank for finite-sample degrees of freedom.
4. Describe differences across centres.
5. Fit steatosis component models in complete-component and imputed data.
6. Describe common profiles versus the remaining participants.
7. Estimate standardized steatosis probabilities with centre-stratified bootstrap intervals.

The executed definitions and statistical implementation are unchanged: Asian adiposity thresholds, the combined glycaemic criterion, recorded hypertension, measured triglycerides and sex-specific strict-less-than HDL-C thresholds. Separate medication and blood-pressure fields were unavailable for reconstruction. ZS source prefixes map to RJ in the derived centre variable. The nonbinary HCV entry is converted to missing by the original MI factor coding.

## Reproduction scope

Reproducing the study estimates requires authorised data access. The plotting scripts also read aggregate intermediate tables and, in some cases, earlier exports; the statistical runner executes the seven analysis steps rather than every manuscript-assembly or historical rendering step. Exact regeneration of all figure files is therefore not asserted by this code-only distribution. Consolidating the public directory layout does not refit the models or alter their estimates.
