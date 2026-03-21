# Chordoma Survival Calculator

A Cox proportional hazards–based survival prediction project for intracranial chordoma using the National Cancer Database (NCDB) Bone Participant Use File.

## Overview

This project develops and reports an at-diagnosis prognostic model for overall survival in adults with intracranial chordoma. The final model uses prespecified baseline predictors and models age and tumor size nonlinearly with restricted cubic splines. Outputs include a complete-case analytic cohort, fitted model objects, manuscript-ready tables and figures, and a Shiny app interface for individual prediction.

## Data source

Analyses use a curated NCDB bone tumor extract stored locally as:

`data/raw/NCDB Bone.csv`

The NCDB is a hospital-based oncology registry jointly sponsored by the American College of Surgeons Commission on Cancer and the American Cancer Society.

## Cohort definition

Eligible patients meet all of the following:
- ICD-O-3 histology/behavior codes: `9370/3`, `9371/3`, or `9372/3`
- Primary site standardized to `C41.0` (bones of skull and face)
- Age at diagnosis `>= 18` years
- Valid survival time and vital status

Tumor size is harmonized across legacy and 2016+ NCDB fields using year-aware selection. Exact millimeter values from `1` to `988` are retained. Tumor size is capped at `70 mm` as a prespecified plausibility constraint.

## Final model

The final model is a multivariable Cox proportional hazards model for overall survival with the following predictors:
- Age at diagnosis
- Sex
- Race
- Hispanic ethnicity
- Insurance category
- Median household income quartile
- Charlson-Deyo comorbidity score
- Tumor size

Age and tumor size are modeled with restricted cubic splines using 4 knots. No automated variable selection is used.

## Validation

Internal validation includes:
- Global apparent and optimism-corrected Harrell’s C-index
- Global apparent and optimism-corrected calibration slope
- Horizon-specific IPCW AUC at 5 and 10 years
- Horizon-specific IPCW Brier score at 5 and 10 years
- Bootstrap optimism correction using 300 resamples

## Project structure

```text
chordoma/
├── R/
│   ├── utils_chordoma.R
│   ├── 01_screen.R
│   ├── 02_clean.R
│   ├── 03_tables_figures.R
│   └── app.R
├── data/
│   ├── raw/
│   │   └── NCDB Bone.csv
│   └── processed/
│       ├── chordoma_screened.csv
│       ├── chordoma_cohort_completecase.csv
│       ├── chordoma_model_objects.rds
│       ├── tables/
│       │   └── manuscript/
│       └── figures/
│           └── manuscript/
├── www/
│   └── ohsu_logo.png
├── README.md
└── .gitignore
