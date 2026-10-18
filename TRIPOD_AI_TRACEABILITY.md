# TRIPOD+AI code-to-report map

This map identifies code and outputs, not final checklist page references. The official TRIPOD+AI checklist should be completed against the paginated accepted manuscript.

| Reporting domain | Source |
|---|---|
| Source population, setting, eligibility, diagnosis years, flow | \`R/01_cohort.R\`; \`outputs/aggregate/cohort_flow.csv\` |
| Outcome definition and censoring | \`R/01_cohort.R\`; \`R/03_validation.R\` |
| Predictor definitions, categorical encodings, 100-mm policy | \`R/01_cohort.R\` and \`R/02_models.R\` |
| Missingness and included-versus-excluded comparisons | \`R/01_cohort.R\`; \`R/07_publication_exhibits.R\` (Tables S1 and S8) |
| Sample size and events, degrees of freedom | \`R/02_models.R\`; \`R/07_publication_exhibits.R\` (Table S7; assumed-R2 scenarios only) |
| Full coefficients, reference groups, prediction equation | \`model/\` and \`R/02_models.R\` |
| Discrimination, prediction error, internal bootstrap | \`R/03_validation.R\` |
| Calibration and subgroup performance | \`R/03_validation.R\`; \`R/06_tables_figures.R\`; \`R/07_publication_exhibits.R\` (Figure 1, Table S4) |
| Socioeconomic predictor comparison | \`R/02_models.R\`; \`R/03_validation.R\` |
| Size-threshold and MI sensitivity | \`R/05_sensitivity.R\`; \`R/07_publication_exhibits.R\` (Tables S3, S5; Figure 2) |
| Clinical implementation limits | \`README.md\` and \`shiny/app.R\` |
| Public availability of code and restricted data | \`CODE_AVAILABILITY.md\` |

Not covered by code: external validation, impact analysis, any independent evidence of ethical or social fairness, journal page/line references, funding and conflict declarations, and author-level research governance statements. These require accurate reporting in the manuscript and submission forms.

## Exhibit coverage and exceptions

The final reporting generator is `R/07_publication_exhibits.R`: main Tables 1–4, supplemental Tables S1–S9, Figures 1–2 and S1–S2. Table S9 is an explicit status disclosure, not a histology/metastasis effect estimate; field-specific review is still needed. Figure S1 is illustrative development-data Kaplan–Meier stratification, not external validation. The generator must be run and reviewed locally before any statement that every publication exhibit is reproduced. All new exhibits remain under Git-ignored `outputs/private/manuscript_exhibits_DISCLOSURE_REVIEW` until NCDB disclosure review.

The 18 October 2026 primary pipeline rerun is described in [REPRODUCIBILITY_AUDIT.md](REPRODUCIBILITY_AUDIT.md). Its bootstrap performance estimates are the designated canonical series for a future manuscript update; do not mix these with earlier random-seed CIs.

## Reviewer closeout: histology and clinical metastatic status

`R/08_histology_metastasis.R` reproduces year-aware clinical-M adjudication, exploratory 9371/3 versus 9370/3 adjustment, and three restricted-cohort clinical-Cox refits. A three-category histology model and an adjusted M1/M0 effect estimate were not supported by rare events. `R/09_histology_reporting.R` generates two Supplementary Table S9 panels from restricted Stage 08 numerical CSVs, withholding rare subgroup counts. The published calculator's five-variable equation and validation are unchanged. These additional GitHub scripts require local numerical verification before release.
