# Investigational Shiny calculator: deployment

The calculator loads the locked Cox equation from aggregate CSV files; it does **not** load NCDB participant records or fitted model RDS files. It is for research and demonstration pending independent external validation.

## Required files

Deploy from the repository root with this **minimal bundle**:

```text
app.R
shiny/app.R
model/coefficients.csv
model/reference_survival.csv
```

`app.R` is the deployment entry point and sources the implementation in `shiny/app.R`. All four files must retain their relative paths.

## Run locally

```r
install.packages(c("shiny", "bslib", "rsconnect"))
shiny::runApp(".")
```

Verify the frozen equation separately:

```sh
Rscript tests/test_equation.R
```

## Deploy

After configuring an authorized shinyapps.io or Posit Connect account, from the repository root:

```r
rsconnect::deployApp(
  appDir = getwd(),
  appPrimaryDoc = "app.R",
  appFiles = c(
    "app.R", "shiny/app.R",
    "model/coefficients.csv", "model/reference_survival.csv"
  )
)
```

Do not deploy the full `R/` development pipeline, `outputs/`, private NCDB files, or a local R package library. The application needs only `shiny` and `bslib` at runtime.

## Deploying from GitHub to Posit Connect Cloud

Select the `main` branch and `app.R` as the primary file. The committed `manifest.json` contains the calculator's runtime dependencies and its deployment files.

Shiny normally loads every top-level R script inside `R/`. Here, that folder contains the licensed-data **analysis pipeline**, which should not run on a hosted calculator. The standard Shiny `R/_disable_autoload.R` file prevents automatic sourcing. There is no need to configure `NCDB_BONE_PUF` or upload the licensed dataset to Connect Cloud. The GitHub Actions smoke test checks the opt-out with that environment variable absent.

## Use restrictions

Inputs must be supplied explicitly; missing values are not imputed. Recorded tumor diameter must be between 1 and 100 mm. Estimates represent **all-cause mortality** at 60 and 120 months in an internally validated, registry-defined cohort with site C41.0. The app does not establish chordoma-specific mortality risk, clinical utility, treatment benefit, or performance in a new hospital or patient population.
