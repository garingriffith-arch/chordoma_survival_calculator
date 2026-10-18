# Code and data availability

The analysis code and aggregate equation required to reconstruct risk calculations are provided in this repository. The 2024 ACS NCDB Bone/Joint Participant User File is governed by its data-use agreement and is **not distributed**. Authorized investigators must obtain the PUF independently and configure the input path using \`NCDB_BONE_PUF\`.

No raw NCDB data, patient-level analytic cohorts, fitted model objects embedding patient information, or resampling checkpoints are intended for this public repository. The \`outputs/\` directory and all \`.rds\`, \`.RData\`, and similar files are excluded from version control.

Run \`R/run_all.R\` from the repository root. Stages 01–06 produce cohort QC, coefficients, equation checks, validation, and sensitivity analyses. Stage 07 contains the generator for the numbered manuscript tables and figures but still needs local execution and audit; its histology/metastasis table is an explicit status disclosure, not completed subgroup estimation. The 18 October 2026 local rerun and its quantitative comparison with locked results are documented in [REPRODUCIBILITY_AUDIT.md](REPRODUCIBILITY_AUDIT.md). Reproduction requires the same PUF release and analysis definitions; locked-equation parity checks fail on incompatible inputs.

The published risk equation is included as aggregate coefficients and baseline survival; its use is investigational. The model has only internal validation. No claims of transportability, clinical utility, or chordoma-specific mortality prediction are made.

**Repository history:** Earlier commits contain legacy fitted-object \`.rds\` files. Removing them from the current tree does not remove them from Git history. Repository history should be cleaned separately if those legacy objects include redistributable restrictions.

**Disclosure controls:** Aggregate CSV output is not automatically safe for redistribution. Sensitivity summaries can reveal very small categories by comparing nested thresholds, and subgroup output is explicitly marked for disclosure review. Do not commit locally generated aggregate output without a data-use and small-cell review.
