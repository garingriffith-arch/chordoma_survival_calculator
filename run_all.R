# Execute from the repository root using Rscript R/run_all.R.
# Data access: environment variable NCDB_BONE_PUF.
source(file.path("R","00_config.R"))
run_stage<-function(script){
  cat("\n== ",script," ==\n",sep="")
  source(file.path("R",script),local=globalenv(),chdir=FALSE)
}
run_stage("01_cohort.R")
run_stage("02_models.R")
run_stage("03_validation.R")
run_stage("04_uncertainty.R")
run_stage("05_sensitivity.R")
run_stage("06_tables_figures.R")
run_stage("07_publication_exhibits.R")
run_stage("08_histology_metastasis.R")
run_stage("09_histology_reporting.R")
cat("\nAll ordered stages finished. Audit outputs before publication.\n")
