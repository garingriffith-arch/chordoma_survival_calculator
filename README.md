# Chordoma survival prediction: 2024 NCDB analysis

Code supporting development and internal validation of a 5- and 10-year all-cause mortality prediction model in adults with chordoma recorded at NCDB primary site C41.0 (bones of skull and face).

**Status:** Manuscript revision; investigational model. Internal validation does not establish accuracy for individual clinical counseling, treatment selection, or use outside the development data. C41.0 is not a verified intracranial or clival location.

## Model

The primary model is a Cox proportional hazards model developed among patients with recorded tumor diameter of **1–100 mm**, without capping or winsorization. Prespecified predictors are age, sex, tumor size per 10 mm, and Charlson–Deyo comorbidity (0, 1, ≥2). A secondary model adds race, Hispanic ethnicity, insurance, and area-level income quartile. Both are compared on a common complete-case cohort. Outcome is all-cause death, with risk estimated at 60 and 120 months.

The 100-mm cutoff follows the SEER Tumor Size Over Time Recode's site-specific C41.0 harmonization range. Recorded sizes above 100 mm are not necessarily invalid; 150-mm and all-exact-size (up to 988-mm) models are included as sensitivity analyses. Tumor-size fields are registry-derived and not assured to be available before surgery.

## Run

Requires R and access to the restricted **2024 NCDB Bone/Joint Participant User File**. No NCDB records are included here.

From the repository root:

~~~r
install.packages(c("readr","survival","prodlim","riskRegression","ggplot2","mice","pmsampsize","patchwork"))
~~~

Set the input file path, then run:

~~~powershell
$env:NCDB_BONE_PUF = "C:\path\to\NCDBPUF_BoneJont.0.2024.0.csv"
Rscript R/run_all.R
~~~

On RStudio, set the working directory to the repository root, set `Sys.setenv(NCDB_BONE_PUF = "...")`, and `source("R/run_all.R")`.

The code **stops** if cohort counts or coefficients fail the supplied locked-results parity checks; a new PUF revision, variable recode, or other changed population requires investigation, not silent acceptance. To run the smaller smoke test use `Sys.setenv(CHORDOMA_QUICK = "true")`; smoke-test output must not be cited in a manuscript.

## Repository layout

| Path | Purpose |
|---|---|
| `R/01_cohort.R` | Import, registry coding, cohort screening and missingness |
| `R/02_models.R` | Primary and secondary Cox models, coefficients and prediction equation |
| `R/03_validation.R` | Bootstrap optimism correction and held-out risk assessment |
| `R/04_uncertainty.R` | Nested-bootstrap uncertainty and paired comparison |\n| `R/05_sensitivity.R` | Tumor-size and missing-data sensitivity analyses |
| `R/06_tables_figures.R` | Reproducible aggregate reporting and figures |
| `R/07_publication_exhibits.R` | Numbered manuscript exhibits (Tables 1–4, S1–S9; Figures 1–2, S1–S2) |
| `R/run_all.R` | Single ordered entry point |
| `model/` | Aggregate coefficients and reference survival; no fitted patient-level RDS |
| `shiny/app.R` | Research-only app using published equation rather than private records |
| `CODE_AVAILABILITY.md` | Data-access, reproducibility, and verification limitations |

Outputs, patient-level derived records, imputation checkpoints, and fitted models containing NCDB data are written **locally** under `outputs/` and are ignored by Git. The code uses no machine-specific Windows paths.


### Generate exhibits without repeating bootstraps

After a completed run of stages 01–06, run from the repository root:

~~~powershell
& "C:\\Program Files\\R\\R-4.6.0\\bin\\Rscript.exe" R/07_publication_exhibits.R
~~~

The tables (CSV, optionally DOCX) and figures (PNG/PDF) are saved under `outputs/private/manuscript_exhibits_DISCLOSURE_REVIEW/`. **Do not publish that directory without reviewing NCDB disclosure restrictions**, especially small-cell counts and differences between size-threshold cohorts. If `pmsampsize` or `patchwork` is missing, install it before running stage 07. The ordered `R/run_all.R` also includes this stage; using it starts the computational pipeline again.

## Investigational Shiny calculator

Deploy the repository root using `app.R`, which loads the styled application from `shiny/app.R`. The app requires only the two published `model/*.csv` files and R packages `shiny` and `bslib`; no patient-level NCDB files or model RDS are used. See [DEPLOYMENT.md](DEPLOYMENT.md) for packaging, local testing and Posit Connect/shinyapps.io steps. Run the data-free equation test with `Rscript tests/test_equation.R`.

## Histology and baseline-stage reviewer analyses

The extended ordered pipeline now includes `R/08_histology_metastasis.R` (era-aware clinical M adjudication, adjusted two-code histology, and primary-model restriction sensitivities) and `R/09_histology_reporting.R` (private, source-driven Supplementary Table S9 panels). Run only these two scripts after the completed 2024 model pipeline:

```powershell
& $Rscript R/08_histology_metastasis.R
& $Rscript R/09_histology_reporting.R
```

The original locked equation is unchanged. These additions have passed static syntax checks only; the public-branch Stage 08/09 implementation still requires a local numerical comparison with the previously completed Stage 08 results. Stage 09 writes only into the Git-ignored disclosure-review folder. Avoid uploading exact rare subtype/stage counts.

## Reproducibility boundary

This repository includes code and a reference prediction equation, not the proprietary NCDB PUF. A local rerun against the designated 2024 PUF completed on 18 October 2026. Its fitted coefficients, reference survival, multiple-imputation coefficients, cohort counts, and apparent performance agree with the locked analysis; corrected resampling estimates differ slightly across random-number streams. See [REPRODUCIBILITY_AUDIT.md](REPRODUCIBILITY_AUDIT.md). No independent external validation is available. Stage 07 includes code for the full numbered exhibit set; the new code requires an independent local execution and disclosure audit before the exhibits can be treated as verified manuscript artifacts. Its S9 reports disclosure and field-adjudication limitations rather than inventing rare histology or metastatic-status counts.

## Reporting standard

Reporting is organized for **TRIPOD+AI (2024)**. The checklist remains a separate submission document that must cite final manuscript pages and line numbers. See [TRIPOD_AI_TRACEABILITY.md](TRIPOD_AI_TRACEABILITY.md) for the code-to-report map.

## Historical files

The earlier repository version used an obsolete 70-mm cap, socioeconomic predictors in the primary equation, and a fitted-object Shiny app. Those artifacts are not part of the revised primary analysis. Historical Git commits can still contain removed files; do not redistribute fitted objects from old revisions.
