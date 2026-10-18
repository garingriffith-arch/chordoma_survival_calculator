# Deploying the investigational chordoma calculator

Deploy from the **repository root**, where `app.R` and the `model/` directory are siblings. The root `app.R` loads the actual implementation in `shiny/app.R`. Do not deploy `shiny/` alone: the two model CSV files are required.

## Files included in the deployment bundle

```text
app.R
shiny/app.R
model/coefficients.csv
model/reference_survival.csv
```

No NCDB Participant User File, analytic cohort, fitted patient-level `.rds` object, bootstrap checkpoint, training script, manuscript, or `outputs/` directory is needed to run the calculator.

## Packages

The app requires **shiny** and **bslib**. The deployment workstation needs **rsconnect** to publish to a Shiny hosting service supported by rsconnect. Package dependencies are installed on the server from the deployment bundle; do not upload a copy of a Windows R library.

Install locally:

```r
install.packages(c("shiny", "bslib", "rsconnect"),
                 repos = "https://cloud.r-project.org")
```

Check the equation without launching Shiny:

```powershell
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" tests/test_equation.R
```

Launch locally from the repository root:

```r
shiny::runApp(".")
```

Deploy from an R session in the repository root after configuring your Posit Connect or shinyapps.io account:

```r
rsconnect::deployApp(
  appDir = getwd(),
  appPrimaryDoc = "app.R",
  appFiles = c("app.R", "shiny/app.R",
               "model/coefficients.csv", "model/reference_survival.csv")
)
```

Check the returned bundle manifest before confirming deployment. If your hosting provider requires additional metadata, use the provider-specific deployment procedure.

## Model and limitations

Only five locked Cox coefficients and reference survival at 60 and 120 months are loaded. The app does not retrain or refit and uses no NCDB participant-level records. Estimates are all-cause mortality for a registry-defined C41.0 cohort with recorded tumor size 1–100 mm. C41.0 does not confirm clival/intracranial anatomy, size may not be measured before treatment, and external validation remains unavailable. Use is investigational, not for clinical decisions.

## Publication and repository hygiene

This draft PR remains unmerged. Historical public commits of the older repository contain fitted model objects that may embed restricted NCDB data; see `HISTORY_REMEDIATION.md`. Resolve historical exposure and small-cell disclosure before citing the default branch as a fully clean published release. Changing the working tree does not clear repository history.
