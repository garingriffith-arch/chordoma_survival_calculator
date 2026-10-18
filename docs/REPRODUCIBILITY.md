# Reproducibility and reporting notes

## Data and cohort

The development dataset was the licensed **2024 NCDB Bone/Joint Participant User File** (diagnoses 2004–2023). Eligible adults had chordoma histology codes 9370/3, 9371/3, or 9372/3 at recorded primary site C41.0, with nonmissing survival information. The primary complete-case cohort comprised **1,243 individuals (310 deaths)** with recorded tumor size 1–100 mm. The common-cohort clinical-versus-expanded comparison used **946 individuals (242 deaths)**.

C41.0 is a registry site grouping, not a verified intracranial or clival classification. All-cause mortality—not disease-specific survival—is the model endpoint.

## Development and internal assessment

A five-degree-of-freedom Cox model used age, tumor diameter per 10 mm, sex, and categorical Charlson–Deyo comorbidity. The expanded comparator added sociodemographic variables. Internally validated performance was estimated using 300 bootstrap replications for optimism correction and a 120-outer-iteration, 30-inner-iteration nested bootstrap for uncertainty, with additional repeated cross-fitting for calibration. A 20-imputation analysis evaluated coefficient sensitivity to missing data.

| Measure | Primary-model optimism-corrected estimate | Nested-bootstrap 95% interval |
|---|---:|---:|
| Harrell C-index | 0.742899 | 0.716571–0.775452 |
| AUC, 60 months | 0.760143 | 0.726152–0.808474 |
| AUC, 120 months | 0.759789 | 0.730599–0.805902 |
| IPCW Brier, 60 months | 0.128481 | 0.113613–0.139181 |
| IPCW Brier, 120 months | 0.176675 | 0.154565–0.191470 |

These estimates arose from the designated local 2024 PUF analysis. Resampled estimates from different random-number streams need not be bitwise identical; reproducibility depends on matching the exact PUF release and code. None of these measures constitutes external validation.

## Histology and metastatic-status sensitivity

An exploratory adjusted analysis compared 9371/3 with 9370/3 after excluding sparse 9372/3: **HR 1.055 (95% CI 0.694–1.602), likelihood-ratio p=0.805**. Histology 9370/3 does not separately identify poorly differentiated disease. Sparse metastatic deaths did not support a stable adjusted M1-versus-M0 estimate. Sensitivity refits excluding 9372/3, excluding documented M1/conflicting stage records, and restricting to documented M0 produced similar clinical-model coefficients. The restriction to documented M0 changes the analyzed population; these are not new validated prediction models.

The year-specific clinical M variables and available metastatic-site flags are reconciled in `R/08_histology_metastasis.R`, with missing status retained as unknown. `R/09_histology_reporting.R` generates publication-oriented summaries **without rare-cell denominators**.

## Privacy and verification

The scripts require the licensed PUF and produce intermediate records locally under `outputs/private/`. Do not publish patient-level files, fitted RDS objects, or disclosure-limited aggregates. In particular, nested tumor-size thresholds can reveal small counts through subtraction. Publication tables should undergo a separate small-cell and complementary-suppression review.

Automated checks parse all public R files, test the locked equation numerically, and load the Shiny application without NCDB data. The end-to-end model pipeline and Stage 08 histology estimate have previously been run against the licensed PUF; these automated checks **cannot** repeat the PUF analysis.

Current limitations include no external validation, no independent clinical utility evaluation, incomplete verification of pre-treatment tumor-size measurement, and no facility-level heterogeneity assessment.
