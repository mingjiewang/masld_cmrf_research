# Restricted workbook interface

The source workbook uses a sheet named `Sheet1`. `source_workbook_columns.csv` contains only the 51 expected column headings; `variable_mapping.csv` maps those headings to the internal analysis variables. They do not contain study observations. Columns that are not used in a particular model remain part of the frozen workbook interface.

The core coding is: sex `0 = male`, `1 = female`; recorded hypertension and T2DM `0 = no`, `1 = yes`; steatosis grades `0–3`; fibrosis stages `0–4`. Source pseudonymous IDs identify the original centre prefix (`HZ`, `JL`, `RJ`, `ZS`); the derived analysis maps `ZS` to `RJ`. No participant IDs are supplied in this repository.

BMI is in kg/m², waist circumference in cm, glucose/triglycerides/HDL-C in mmol/L, HbA1c in percent, and age in years. The exact preprocessing, invalid-zero handling, conservative domain status reconstruction and binary-factor coercion are retained in the source scripts. Clinical eligibility must be interpreted in the context of the study protocol; this interface does not independently certify protocol execution.

The source workbook and all generated participant-level CSV/RDS files require the relevant institutional permission. The runner creates them only in the explicitly selected private work directory outside the repository. Do not commit them or replace these empty templates with study observations.
