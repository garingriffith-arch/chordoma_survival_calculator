# Deploy from the repository root so both model CSV files are included.
# This entry point keeps the implementation in shiny/app.R.
if (!file.exists(file.path("model","coefficients.csv")) ||
    !file.exists(file.path("model","reference_survival.csv")))
  stop("Deploy or launch the entire repository root, not the shiny/ directory alone.")
source(file.path("shiny","app.R"),local=TRUE)$value
