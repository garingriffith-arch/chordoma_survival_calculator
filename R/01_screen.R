source("R/utils_chordoma.R")

ensure_project_dirs()

df0 <- read_ncdb_bone()
df1 <- derive_screen_variables(df0)
df2 <- apply_cohort_filters(df1)

write_csv(df2, file.path(PROCESSED_DIR, "chordoma_screened.csv"))

cat("Rows in raw file:", nrow(df0), "\n")
cat("Rows after screening:", nrow(df2), "\n")
cat("Tumor size capped at", TUMOR_SIZE_CAP_MM, "mm:", sum(df2$tumor_size_was_capped == 1L, na.rm = TRUE), "\n")