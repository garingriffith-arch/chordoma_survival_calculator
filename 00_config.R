# Run scripts from the repository root.
ROOT <- normalizePath(Sys.getenv("CHORDOMA_ROOT", unset=getwd()),
                      winslash="/", mustWork=TRUE)
if (!file.exists(file.path(ROOT,"R","run_all.R"))) {
  stop("Set the working directory to the chordoma repository root.")
}
PUF_FILE <- Sys.getenv("NCDB_BONE_PUF", unset="")
if (!nzchar(PUF_FILE)) {
  stop("Set NCDB_BONE_PUF to the locally licensed 2024 Bone/Joint PUF CSV.")
}
if (!file.exists(PUF_FILE)) stop("NCDB_BONE_PUF file does not exist.")
OUTPUT_ROOT <- file.path(ROOT,"outputs")
PRIVATE_DIR <- file.path(OUTPUT_ROOT,"private")
AGG_DIR <- file.path(OUTPUT_ROOT,"aggregate")
FIG_DIR <- file.path(OUTPUT_ROOT,"figures")
for (d in c(PRIVATE_DIR,AGG_DIR,FIG_DIR)) dir.create(d, recursive=TRUE,showWarnings=FALSE)
MAX_MM <- 100L
HORIZONS <- c(60,120)
SEED <- 20261017L
QUICK <- identical(tolower(Sys.getenv("CHORDOMA_QUICK","false")),"true")
B_BOOT <- if (QUICK) 12L else 300L
CV_REPEATS <- if (QUICK) 2L else 10L
CV_FOLDS <- 5L
MI_M <- if (QUICK) 3L else 20L
PKGS <- c("readr","survival","prodlim","riskRegression","ggplot2","mice")
missing_pkgs <- PKGS[!vapply(PKGS,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing_pkgs)) stop("Install R packages: ",paste(missing_pkgs,collapse=", "))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(prodlim))
options(stringsAsFactors=FALSE)
write_aggregate <- function(x, filename) {
  utils::write.csv(as.data.frame(x),file.path(AGG_DIR,filename),row.names=FALSE,na="")
}
write_private <- function(x, filename) {
  saveRDS(x,file.path(PRIVATE_DIR,filename))
}
