# Reporting traceability: TRIPOD+AI

This is a map from the analysis code to reporting evidence. It is **not** a replacement for the completed TRIPOD+AI checklist submitted with a manuscript, and it does not imply all validation criteria were satisfied.

| Reporting area | Analysis source |
|---|---|
| NCDB source, diagnosis years, eligibility, cohort flow | `R/01_cohort.R`; `R/07_publication_exhibits.R` |
| Outcome, censoring and follow-up | `R/01_cohort.R`; `R/03_validation.R` |
| Predictor coding and missingness | `R/01_cohort.R`; `R/07_publication_exhibits.R` |
| Model formulation and reference survival | `R/02_models.R`; `model/coefficients.csv`; `model/reference_survival.csv` |
| Optimism-corrected discrimination and error | `R/03_validation.R`; `R/04_uncertainty.R` |
| Calibration and subgroup assessment | `R/03_validation.R`; `R/06_tables_figures.R`; `R/07_publication_exhibits.R` |
| Comparator models and missing-data sensitivity | `R/02_models.R`; `R/05_sensitivity.R` |
| Histology and metastatic-status sensitivity | `R/08_histology_metastasis.R`; `R/09_histology_reporting.R` |
| Manuscript tables and figures | `R/07_publication_exhibits.R`; `R/09_histology_reporting.R` |
| Public code, restricted data and deployment | `CODE_AVAILABILITY.md`; `DEPLOYMENT.md` |

The study was not prospectively registered and had no prior formal written protocol. No formal IRB determination was sought because IRB review was not required for the research. Patient/public representatives were not involved in the research. There was no project-specific funding or declared conflict of interest. Predictor values were obtained from registry records, and sociodemographic subgroup evaluation was exploratory.

Not performed: independent external validation, center-level performance heterogeneity analysis, or prospective clinical-impact assessment. The intended clinical audience is conditional on future validation and appropriate clinical-utility studies. The publicly available calculator refuses incomplete predictor inputs rather than imputing missing values.

For item-by-item checklist status and final manuscript page locations, refer to the submission's separately maintained **52-item TRIPOD+AI checklist**.
