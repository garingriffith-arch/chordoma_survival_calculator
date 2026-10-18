# Chordoma mortality prediction from the National Cancer Database

Analysis code and an investigational Shiny calculator for estimating 5- and 10-year **all-cause mortality** among adults with chordoma recorded at NCDB primary site **C41.0 (bones of the skull and face)**. The analysis uses the **2024 NCDB Bone/Joint Participant User File**, covering diagnoses in **2004–2023**.

**Development and internal validation only.** The model has not undergone independent external validation or a prospective clinical-impact evaluation. It is not intended to direct treatment decisions. The C41.0 registry site does not independently establish clival, intracranial, or skull-base involvement.

## Study and model

The complete-case clinical-model cohort comprised **1,243 patients and 310 deaths**. The Cox proportional hazards model uses five degrees of freedom: age at diagnosis, recorded tumor diameter (per 10 mm), sex, and Charlson–Deyo comorbidity score (0, 1, or ≥2). Recorded tumor size was restricted to **1–100 mm** for the primary analysis, with larger-size-range sensitivity analyses. Values above 100 mm were excluded, not truncated.

| Internal-validation measure | Optimism-corrected estimate (95% CI) |
|---|---|
| Harrell C-index | 0.743 (0.717–0.775) |
| 5-year time-dependent AUC | 0.760 (0.726–0.808) |
| 10-year time-dependent AUC | 0.760 (0.731–0.806) |
| 5-year IPCW Brier score | 0.128 (0.114–0.139) |
| 10-year IPCW Brier score | 0.177 (0.155–0.191) |

A larger comparator model also included registry-coded race, Hispanic ethnicity, insurance, and area-level income quartile; the comparison used the same 946-patient complete-case cohort. This research does not establish causal effects or treatment benefit.

## Repository contents

| Path | Contents |
|---|---|
| `R/01_cohort.R`–`R/06_tables_figures.R` | Cohort construction, models, internal validation, uncertainty and sensitivity analyses |
| `R/07_publication_exhibits.R` | Main and supplementary tables and figures |
| `R/08_histology_metastasis.R` | Exploratory coded-histology and baseline-metastasis sensitivity analyses |
| `R/09_histology_reporting.R` | Source-derived, disclosure-controlled Supplementary Table S9 |
| `R/run_all.R` | Ordered local analysis pipeline |
| `model/*.csv` | Locked aggregate Cox coefficients and reference survival |
| `app.R`, `shiny/app.R` | Calculator entry point and user interface |
| `tests/test_equation.R` | Data-free prediction-equation tests |
| `docs/REPRODUCIBILITY.md` | Reproduction details and interpretive limits |
| `docs/TRIPOD_AI_TRACEABILITY.md` | Code-to-report map |

## Reproducing the analyses

The **restricted NCDB PUF is not included**. Authorized researchers must independently obtain the same PUF release and comply with its data-use terms. From the repository root, with R installed:

```r
install.packages(c(
  "readr", "survival", "prodlim", "riskRegression",
  "ggplot2", "mice", "pmsampsize", "patchwork",
  "officer", "flextable"
))
```

```powershell
$env:NCDB_BONE_PUF = "C:\path\to\NCDBPUF_BoneJont.0.2024.0.csv"
Rscript R/run_all.R
```

The full pipeline includes computationally intensive bootstrap analyses. To regenerate the final reviewer sensitivity table after running the prior stages:

```powershell
Rscript R/08_histology_metastasis.R
Rscript R/09_histology_reporting.R
```

Intermediate fitted objects, derived patient-level datasets and publication exhibits are written under the local **Git-ignored** `outputs/` directory. **Do not commit outputs without NCDB disclosure review**: even aggregate counts may disclose small groups through subtraction across nested cohorts.

## Calculator

The Shiny app reads only aggregate model coefficients and 60-/120-month reference survival. It requires complete patient inputs and will not estimate risk for recorded sizes outside 1–100 mm. It does not retrain the model or need access to NCDB records.

```r
install.packages(c("shiny", "bslib"))
shiny::runApp(".")
```

The locked equation can be checked independently of NCDB data or Shiny:

```sh
Rscript tests/test_equation.R
```

See [DEPLOYMENT.md](DEPLOYMENT.md) for a minimal deployment bundle and [CODE_AVAILABILITY.md](CODE_AVAILABILITY.md) for data-access conditions.

## Interpretation

The endpoint is death **from any cause**, not chordoma-specific mortality, recurrence, or treatment response. Registry-recorded tumor diameter may not be a pretreatment measurement. Missing or miscoded predictor data, selection of complete cases, temporal changes in treatment, and the absence of external validation limit transportability. This research calculator is intended for investigation and further validation, not individual clinical decision-making.
