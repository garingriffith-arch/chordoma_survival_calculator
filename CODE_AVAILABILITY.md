# Code and data availability

Source code for cohort construction, Cox model development, internal validation, uncertainty estimation, sensitivity analyses, and manuscript exhibit generation is provided in this repository. The frozen aggregate prediction coefficients and reference survival values required by the calculator are located under `model/`.

The **2024 National Cancer Database (NCDB) Bone/Joint Participant User File** is subject to an individual data-use agreement. The participant-level PUF cannot be redistributed through this repository. Authorized investigators must obtain access separately and provide the local PUF path through the `NCDB_BONE_PUF` environment variable.

The public repository intentionally contains **no** PUF extract, analytic patient-level dataset, fitted Cox object with embedded records, bootstrap checkpoint, or generated `outputs/` directory. Outputs from the R pipeline are excluded by `.gitignore` and require separate NCDB disclosure review before dissemination. Apparent aggregate output is not automatically disclosure-safe, particularly when comparing nested, small cohort differences.

The locked model equation can be tested without access to the PUF using `Rscript tests/test_equation.R`. Reproducing the development and validation results requires the specified PUF release, the code in `R/`, and compatible R package versions. The Shiny calculator uses only the aggregate model files and has **not** been externally validated or evaluated for clinical impact.
