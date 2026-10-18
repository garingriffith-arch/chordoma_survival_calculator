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

## Use restrictions

Inputs must be supplied explicitly; missing values are not imputed. Recorded tumor diameter must be between 1 and 100 mm. Estimates represent **all-cause mortality** at 60 and 120 months in an internally validated, registry-defined cohort with site C41.0. The app does not establish chordoma-specific mortality risk, clinical utility, treatment benefit, or performance in a new hospital or patient population.
