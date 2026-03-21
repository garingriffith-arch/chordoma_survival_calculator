source("R/utils_chordoma.R")

ensure_project_dirs()

set.seed(20260316)
B_BOOT <- 300

df0 <- read_ncdb_bone()
df1 <- derive_screen_variables(df0)
df2_base <- apply_cohort_filters(df1)
df2 <- make_complete_case_model_df(df2_base)

cat("N after complete-case restriction:", nrow(df2), "\n")
cat("Events after complete-case restriction:", sum(df2$event == 1), "\n")

set_datadist(df2)
form <- model_formula()

cox_fit <- coxph(form, data = df2, x = TRUE, y = TRUE, ties = "efron", model = TRUE)
cph_fit <- rms::cph(form, data = df2, x = TRUE, y = TRUE, surv = TRUE)

lp <- predict(cox_fit, type = "lp")
app_c <- survival::concordance(Surv(time_months, event) ~ I(-lp), data = df2)$concordance

val <- rms::validate(cph_fit, method = "boot", B = B_BOOT, dxy = TRUE)
c_corr <- as.numeric(val["Dxy", "index.corrected"] / 2 + 0.5)

cat("Apparent C-index:", round(app_c, 3), "\n")
cat("Optimism-corrected C-index:", round(c_corr, 3), "\n")

write_csv(df2, file.path(PROCESSED_DIR, "chordoma_cohort_completecase.csv"))

saveRDS(
  list(
    df2 = df2,
    cox_fit = cox_fit,
    cph_fit = cph_fit,
    validate = val,
    apparent_cindex = app_c,
    optimism_corrected_cindex = c_corr
  ),
  file.path(PROCESSED_DIR, "chordoma_model_objects.rds")
)