# R/03_tables_figures.R
# Manuscript tables and figures for the chordoma Cox model

source("R/utils_chordoma.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(survival)
  library(rms)
  library(ggplot2)
  library(riskRegression)
  library(prodlim)
  library(scales)
})

ensure_project_dirs()

SEED <- 20260316
B_CINDEX <- 300
B_STRICT <- 300
HORIZONS_M <- c(`5-year` = 60, `10-year` = 120)
times_vec <- as.numeric(HORIZONS_M)

tables_dir <- file.path(PROCESSED_DIR, "tables", "manuscript")
figures_dir <- file.path(PROCESSED_DIR, "figures", "manuscript")
base_rds <- file.path(PROCESSED_DIR, "chordoma_model_objects.rds")

stopifnot(file.exists(base_rds))

# ----------------------------
# Helpers
# ----------------------------
coerce_event01 <- function(x) {
  x_chr <- trimws(as.character(x))
  out <- dplyr::case_when(
    x_chr %in% c("0", "Alive", "alive", "FALSE", "False", "false", "No", "no") ~ 0L,
    x_chr %in% c("1", "Dead", "dead", "TRUE", "True", "true", "Yes", "yes") ~ 1L,
    TRUE ~ suppressWarnings(as.integer(x_chr))
  )
  out
}

coerce_numeric_safely <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# ----------------------------
# Table 1 and Table 1b: cohort flow and missingness
# ----------------------------
df0 <- read_ncdb_bone()
df1 <- derive_screen_variables(df0)

step_tbl <- function(name, before, after) {
  tibble(
    step = name,
    N_before = nrow(before),
    N_after = nrow(after),
    excluded = nrow(before) - nrow(after)
  )
}

flow <- list()
flow[[1]] <- tibble(
  step = "Start: raw NCDB file",
  N_before = nrow(df0),
  N_after = nrow(df0),
  excluded = 0
)

a1 <- df1 %>% filter(PRIMARY_SITE_STD == SITE_KEEP)
flow[[2]] <- step_tbl("Restrict to primary site C41.0", df1, a1)

a2 <- a1 %>% filter(HIST_BEH %in% HIST_KEEP)
flow[[3]] <- step_tbl("Restrict to chordoma histology/behavior (9370/3, 9371/3, 9372/3)", a1, a2)

a3 <- a2 %>% filter(!is.na(AGE_NUM) & AGE_NUM >= 18)
flow[[4]] <- step_tbl("Exclude age <18 or missing age", a2, a3)

a4 <- a3 %>% filter(
  !is.na(time_months) & time_months >= 0,
  vital_status %in% c(0L, 1L),
  !is.na(event)
)
flow[[5]] <- step_tbl("Exclude missing/invalid survival time or vital status", a3, a4)

a5 <- recode_model_variables(a4)
a6 <- a5 %>% filter(complete.cases(across(all_of(MODEL_VARS))))
flow[[6]] <- step_tbl("Complete-case restriction on all model variables", a5, a6)

write_csv(bind_rows(flow), file.path(tables_dir, "Table1_CohortFlow.csv"))

missing_tbl <- a5 %>%
  summarise(across(all_of(MODEL_VARS), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(100 * n_missing / nrow(a5), 1))

write_csv(missing_tbl, file.path(tables_dir, "Table1b_Missingness.csv"))

# ----------------------------
# Load analysis cohort and refit models in this session
# ----------------------------
obj <- readRDS(base_rds)

df <- obj$df2 %>%
  as.data.frame() %>%
  dplyr::select(
    time_months, event, AGE_NUM, SEX3, RACE3, HISP2,
    INS3, MED_INC_QUAR, cdcc, tumor_size_mm
  ) %>%
  dplyr::mutate(
    time_months = coerce_numeric_safely(time_months),
    event = coerce_event01(event),
    AGE_NUM = coerce_numeric_safely(AGE_NUM),
    tumor_size_mm = coerce_numeric_safely(tumor_size_mm),
    cdcc = coerce_numeric_safely(cdcc),
    SEX3 = factor(as.character(SEX3), levels = c("Male", "Female")),
    RACE3 = factor(as.character(RACE3), levels = c("White", "Black", "Other")),
    HISP2 = factor(as.character(HISP2), levels = c("Non-Hispanic", "Hispanic")),
    INS3 = factor(as.character(INS3), levels = c("Private", "Medicare", "Medicaid/OtherGov", "Uninsured")),
    MED_INC_QUAR = factor(as.character(MED_INC_QUAR), levels = c("1", "2", "3", "4"))
  )

stopifnot(is.numeric(df$time_months))
stopifnot(all(stats::na.omit(unique(df$event)) %in% c(0L, 1L)))

cox_form <- survival::Surv(time_months, event) ~
  rms::rcs(AGE_NUM, 4) +
  rms::rcs(tumor_size_mm, 4) +
  SEX3 + RACE3 + HISP2 + INS3 + MED_INC_QUAR + cdcc

dd <- rms::datadist(df[, c(
  "time_months", "event", "AGE_NUM", "SEX3", "RACE3", "HISP2",
  "INS3", "MED_INC_QUAR", "cdcc", "tumor_size_mm"
), drop = FALSE])
options(datadist = "dd")

cox_fit <- survival::coxph(
  cox_form,
  data = df,
  ties = "efron",
  x = TRUE,
  y = TRUE,
  model = TRUE,
  singular.ok = TRUE
)

cph_fit <- rms::cph(
  cox_form,
  data = df,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  model = TRUE,
  time.inc = max(times_vec),
  singular.ok = TRUE
)

# ----------------------------
# Table 2: baseline characteristics
# ----------------------------
tab2_cont <- tibble(
  Characteristic = c("Sample size", "Deaths", "Follow-up, months", "Age, years", "Tumor size, mm"),
  Summary = c(
    as.character(nrow(df)),
    as.character(sum(df$event == 1, na.rm = TRUE)),
    fmt_med_iqr(df$time_months),
    paste0(fmt_med_iqr(df$AGE_NUM), " ; mean (SD): ", fmt_mean_sd(df$AGE_NUM)),
    paste0(fmt_med_iqr(df$tumor_size_mm), " ; mean (SD): ", fmt_mean_sd(df$tumor_size_mm))
  )
)

tab2_cat <- bind_rows(
  tibble(
    Characteristic = paste0("Sex: ", levels(df$SEX3)),
    Summary = map_chr(levels(df$SEX3), ~ fmt_n_pct(df$SEX3 == .x))
  ),
  tibble(
    Characteristic = paste0("Race: ", levels(df$RACE3)),
    Summary = map_chr(levels(df$RACE3), ~ fmt_n_pct(df$RACE3 == .x))
  ),
  tibble(
    Characteristic = paste0("Ethnicity: ", levels(df$HISP2)),
    Summary = map_chr(levels(df$HISP2), ~ fmt_n_pct(df$HISP2 == .x))
  ),
  tibble(
    Characteristic = paste0("Insurance: ", levels(df$INS3)),
    Summary = map_chr(levels(df$INS3), ~ fmt_n_pct(df$INS3 == .x))
  ),
  tibble(
    Characteristic = paste0("Income quartile: ", levels(df$MED_INC_QUAR)),
    Summary = map_chr(levels(df$MED_INC_QUAR), ~ fmt_n_pct(df$MED_INC_QUAR == .x))
  ),
  tibble(
    Characteristic = c("CDCC: 0", "CDCC: 1", "CDCC: 2+"),
    Summary = c(
      fmt_n_pct(df$cdcc == 0),
      fmt_n_pct(df$cdcc == 1),
      fmt_n_pct(df$cdcc >= 2)
    )
  )
)

write_csv(bind_rows(tab2_cont, tab2_cat), file.path(tables_dir, "Table2_BaselineCharacteristics.csv"))

# ----------------------------
# Table 3: final model
# ----------------------------
s <- summary(cox_fit)

coef_tab <- as.data.frame(s$coefficients)
coef_tab$term <- rownames(coef_tab)

ci_tab <- as.data.frame(suppressWarnings(confint(cox_fit)))
ci_tab$term <- rownames(ci_tab)

if (ncol(ci_tab) >= 2) {
  names(ci_tab)[1:2] <- c("ci_low", "ci_high")
} else {
  stop("confint(cox_fit) did not return the expected two columns.")
}

coef_tab <- dplyr::left_join(coef_tab, ci_tab, by = "term")
rn <- coef_tab$term

is_spline_basis <- grepl(
  "rcs\\(AGE_NUM, 4\\)|rcs\\(tumor_size_mm, 4\\)|AGE_NUM'|tumor_size_mm'",
  rn
)
hr_rows <- which(!is_spline_basis)

hr_part <- tibble(
  Predictor = map_chr(rn[hr_rows], clean_term_label, df = df),
  HR = round(exp(coef_tab$coef[hr_rows]), 2),
  CI_95 = paste0(
    sprintf("%.2f", exp(coef_tab$ci_low[hr_rows])),
    "–",
    sprintf("%.2f", exp(coef_tab$ci_high[hr_rows]))
  ),
  P_value = fmt_p(coef_tab$`Pr(>|z|)`[hr_rows])
)

ref_nd <- make_reference_row(df)

age50 <- ref_nd
age50$AGE_NUM <- 50

age70 <- ref_nd
age70$AGE_NUM <- 70

sz30 <- ref_nd
sz30$tumor_size_mm <- 30

sz60 <- ref_nd
sz60$tumor_size_mm <- 60

get_contrast_est <- function(cph_fit, new_hi, new_lo) {
  con <- rms::contrast(cph_fit, new_hi, new_lo)
  
  if (!is.null(con$Contrast) && !is.null(con$SE)) {
    return(list(est = as.numeric(con$Contrast[1]), se = as.numeric(con$SE[1])))
  }
  
  con_u <- unclass(con)
  est_name <- names(con_u)[grepl("contrast|effect", names(con_u), ignore.case = TRUE)][1]
  se_name  <- names(con_u)[grepl("^se$|se|std", names(con_u), ignore.case = TRUE)][1]
  
  if (!is.na(est_name) && !is.na(se_name)) {
    return(list(
      est = as.numeric(con_u[[est_name]][1]),
      se = as.numeric(con_u[[se_name]][1])
    ))
  }
  
  stop("Could not extract estimate/SE from rms::contrast() output.")
}

age_con <- get_contrast_est(cph_fit, age70, age50)
sz_con  <- get_contrast_est(cph_fit, sz60, sz30)

contrast_part <- tibble(
  Predictor = c("Age: 70 vs 50 years", "Tumor size: 60 vs 30 mm"),
  HR = round(exp(c(age_con$est, sz_con$est)), 2),
  CI_95 = c(
    paste0(
      sprintf("%.2f", exp(age_con$est - 1.96 * age_con$se)),
      "–",
      sprintf("%.2f", exp(age_con$est + 1.96 * age_con$se))
    ),
    paste0(
      sprintf("%.2f", exp(sz_con$est - 1.96 * sz_con$se)),
      "–",
      sprintf("%.2f", exp(sz_con$est + 1.96 * sz_con$se))
    )
  ),
  P_value = c(NA_character_, NA_character_)
)

write_csv(bind_rows(contrast_part, hr_part), file.path(tables_dir, "Table3_FinalModel.csv"))

# ----------------------------
# Table 4: global internal validation
# ----------------------------
lp <- predict(cox_fit, type = "lp")
c_app <- survival::concordance(
  survival::Surv(time_months, event) ~ I(-lp),
  data = df
)$concordance

set.seed(SEED)
val_global <- rms::validate(cph_fit, method = "boot", B = B_CINDEX, dxy = TRUE)

global_c_corr <- as.numeric(val_global["Dxy", "index.corrected"] / 2 + 0.5)
global_c_opt <- c_app - global_c_corr

global_slope_app <- as.numeric(val_global["Slope", "index.orig"])
global_slope_corr <- as.numeric(val_global["Slope", "index.corrected"])
global_slope_opt <- global_slope_app - global_slope_corr

tab4 <- tibble(
  Metric = c("Global Harrell C-index", "Global calibration slope"),
  Apparent = c(round(c_app, 3), round(global_slope_app, 3)),
  Optimism = c(round(global_c_opt, 3), round(global_slope_opt, 3)),
  OptimismCorrected = c(round(global_c_corr, 3), round(global_slope_corr, 3)),
  BootstrapResamples = B_CINDEX
)

write_csv(tab4, file.path(tables_dir, "Table4_GlobalValidation.csv"))

# ----------------------------
# Table 5: strict horizon-specific validation
# ----------------------------
factor_vars <- names(df)[vapply(df, is.factor, logical(1))]
factor_levels <- lapply(df[factor_vars], levels)

df_scored <- apply_factor_levels(df, factor_vars, factor_levels) %>%
  as.data.frame()

df_scored$time_months <- coerce_numeric_safely(df_scored$time_months)
df_scored$event <- coerce_event01(df_scored$event)

stopifnot(all(stats::na.omit(unique(df_scored$event)) %in% c(0L, 1L)))

get_score_table <- function(score_obj, slot_name) {
  slot <- score_obj[[slot_name]]
  if (is.null(slot)) stop("Score() output missing slot: ", slot_name)
  if (is.list(slot) && !is.null(slot$score)) {
    out <- as.data.frame(slot$score)
  } else {
    out <- as.data.frame(slot)
  }
  tibble::as_tibble(out)
}

find_metric_col <- function(tbl, candidates) {
  nm <- names(tbl)
  hit <- nm[tolower(nm) %in% tolower(candidates)]
  if (length(hit) == 0) {
    hit <- nm[grepl(paste(candidates, collapse = "|"), nm, ignore.case = TRUE)]
  }
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

extract_auc_brier <- function(score_obj, times_vec, model_name = "CoxPH") {
  auc_tbl <- get_score_table(score_obj, "AUC")
  brier_tbl <- get_score_table(score_obj, "Brier")
  
  time_col_auc <- find_metric_col(auc_tbl, c("times", "time"))
  time_col_brier <- find_metric_col(brier_tbl, c("times", "time"))
  model_col_auc <- find_metric_col(auc_tbl, c("model"))
  model_col_brier <- find_metric_col(brier_tbl, c("model"))
  auc_col <- find_metric_col(auc_tbl, c("auc", "score"))
  brier_col <- find_metric_col(brier_tbl, c("brier", "score"))
  
  if (any(is.na(c(time_col_auc, time_col_brier, model_col_auc, model_col_brier, auc_col, brier_col)))) {
    stop("Could not parse riskRegression::Score() output columns.")
  }
  
  auc_out <- auc_tbl %>%
    filter(.data[[model_col_auc]] == model_name, .data[[time_col_auc]] %in% times_vec) %>%
    transmute(
      horizon_m = as.numeric(.data[[time_col_auc]]),
      auc = as.numeric(.data[[auc_col]])
    ) %>%
    arrange(horizon_m)
  
  brier_out <- brier_tbl %>%
    filter(.data[[model_col_brier]] == model_name, .data[[time_col_brier]] %in% times_vec) %>%
    transmute(
      horizon_m = as.numeric(.data[[time_col_brier]]),
      brier = as.numeric(.data[[brier_col]])
    ) %>%
    arrange(horizon_m)
  
  full_join(auc_out, brier_out, by = "horizon_m") %>%
    arrange(horizon_m)
}

strict_horizon_metrics <- function(fit, data, times_vec) {
  data <- as.data.frame(data)
  
  sc <- riskRegression::Score(
    object = list(CoxPH = fit),
    formula = Hist(time_months, event) ~ 1,
    data = data,
    cause = 1,
    metrics = c("auc", "brier"),
    times = times_vec,
    cens.model = "km",
    conf.int = FALSE,
    null.model = FALSE,
    summary = "risk"
  )
  
  extract_auc_brier(sc, times_vec = times_vec, model_name = "CoxPH") %>%
    mutate(
      Total_analytic_N = nrow(data),
      Deaths_by_horizon = map_dbl(
        horizon_m,
        ~ sum(data$event == 1 & data$time_months <= .x, na.rm = TRUE)
      )
    ) %>%
    select(horizon_m, Total_analytic_N, Deaths_by_horizon, auc, brier)
}

fit_coxph_boot <- function(data_boot) {
  data_boot <- apply_factor_levels(data_boot, factor_vars, factor_levels)
  data_boot <- as.data.frame(data_boot)
  data_boot$time_months <- coerce_numeric_safely(data_boot$time_months)
  data_boot$event <- coerce_event01(data_boot$event)
  
  tryCatch(
    survival::coxph(
      cox_form,
      data = data_boot,
      ties = "efron",
      x = TRUE,
      y = TRUE,
      model = TRUE,
      singular.ok = TRUE
    ),
    error = function(e) NULL
  )
}

bootstrap_strict_horizon_validation <- function(data, times_vec, B = 300, seed = 1) {
  set.seed(seed)
  
  data <- as.data.frame(data)
  data$time_months <- coerce_numeric_safely(data$time_months)
  data$event <- coerce_event01(data$event)
  
  n <- nrow(data)
  
  fit_orig <- survival::coxph(
    cox_form,
    data = data,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE,
    singular.ok = TRUE
  )
  
  apparent <- strict_horizon_metrics(fit_orig, data, times_vec)
  
  boot_list <- vector("list", B)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot_dat <- apply_factor_levels(data[idx, , drop = FALSE], factor_vars, factor_levels)
    boot_dat <- as.data.frame(boot_dat)
    boot_dat$time_months <- coerce_numeric_safely(boot_dat$time_months)
    boot_dat$event <- coerce_event01(boot_dat$event)
    
    fit_boot <- fit_coxph_boot(boot_dat)
    if (is.null(fit_boot)) next
    
    m_in <- tryCatch(strict_horizon_metrics(fit_boot, boot_dat, times_vec), error = function(e) NULL)
    m_out <- tryCatch(strict_horizon_metrics(fit_boot, data, times_vec), error = function(e) NULL)
    if (is.null(m_in) || is.null(m_out)) next
    
    boot_list[[b]] <- m_in %>%
      select(horizon_m, auc_in = auc, brier_in = brier) %>%
      inner_join(
        m_out %>% select(horizon_m, auc_out = auc, brier_out = brier),
        by = "horizon_m"
      ) %>%
      mutate(
        auc_optimism = auc_in - auc_out,
        brier_optimism = brier_in - brier_out,
        boot_id = b
      ) %>%
      select(boot_id, horizon_m, auc_optimism, brier_optimism)
  }
  
  boot_res <- bind_rows(boot_list)
  
  if (nrow(boot_res) == 0) {
    stop("All bootstrap strict-horizon validation iterations failed.")
  }
  
  optimism_tbl <- boot_res %>%
    group_by(horizon_m) %>%
    summarise(
      auc_optimism = mean(auc_optimism, na.rm = TRUE),
      brier_optimism = mean(brier_optimism, na.rm = TRUE),
      n_successful_boot = n_distinct(boot_id),
      .groups = "drop"
    )
  
  apparent %>%
    left_join(optimism_tbl, by = "horizon_m") %>%
    mutate(
      Horizon = paste0(horizon_m / 12, "-year"),
      AUC_apparent = auc,
      AUC_optimism = auc_optimism,
      AUC_corrected = auc - auc_optimism,
      Brier_apparent = brier,
      Brier_optimism = brier_optimism,
      Brier_corrected = brier - brier_optimism
    ) %>%
    select(
      Horizon,
      Total_analytic_N,
      Deaths_by_horizon,
      AUC_apparent,
      AUC_optimism,
      AUC_corrected,
      Brier_apparent,
      Brier_optimism,
      Brier_corrected,
      n_successful_boot
    ) %>%
    mutate(across(
      c(AUC_apparent, AUC_optimism, AUC_corrected,
        Brier_apparent, Brier_optimism, Brier_corrected),
      ~ round(.x, 3)
    ))
}

tab5 <- bootstrap_strict_horizon_validation(
  data = df_scored,
  times_vec = times_vec,
  B = B_STRICT,
  seed = SEED + 1
)

write_csv(tab5, file.path(tables_dir, "Table5_HorizonValidation_Strict.csv"))

# ----------------------------
# Figures
# ----------------------------
pred_risk <- setNames(
  lapply(times_vec, function(tt) predict_risk_horizon_coxph(cox_fit, df, tt)),
  as.character(times_vec)
)

# Figure 1: calibration
bins_all <- map_dfr(times_vec, function(tt) {
  calibration_bins(df$time_months, df$event, pred_risk[[as.character(tt)]], horizon = tt, groups = 5)
}) %>%
  mutate(
    horizon = factor(horizon_m, levels = times_vec, labels = names(HORIZONS_M))
  )

p_cal <- ggplot(bins_all, aes(pred_mean, obs_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n), alpha = 0.95) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  facet_wrap(~ horizon, nrow = 1) +
  labs(
    title = "Calibration of predicted mortality risk",
    x = "Predicted mortality risk",
    y = "Observed mortality risk (Kaplan-Meier)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(figures_dir, "Figure1_Calibration_5y_10y.png"), p_cal, width = 9, height = 4.8, dpi = 300)
ggsave(file.path(figures_dir, "Figure1_Calibration_5y_10y.pdf"), p_cal, width = 9, height = 4.8)

# Figure 2: age spline
p_age <- as.data.frame(rms::Predict(cph_fit, AGE_NUM, fun = exp, ref.zero = TRUE))
age_ref <- median(df$AGE_NUM, na.rm = TRUE)

g_age <- ggplot(p_age, aes(x = AGE_NUM, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.7) +
  labs(
    title = "Adjusted hazard ratio vs age",
    subtitle = paste0("Reference age = median age (", round(age_ref, 1), " years)"),
    x = "Age (years)",
    y = "Hazard ratio"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(figures_dir, "Figure2_Spline_Age.png"), g_age, width = 7, height = 4.8, dpi = 300)
ggsave(file.path(figures_dir, "Figure2_Spline_Age.pdf"), g_age, width = 7, height = 4.8)

# Figure 3: tumor size spline
p_size <- as.data.frame(rms::Predict(cph_fit, tumor_size_mm, fun = exp, ref.zero = TRUE))
size_ref <- median(df$tumor_size_mm, na.rm = TRUE)

g_size <- ggplot(p_size, aes(x = tumor_size_mm, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.7) +
  labs(
    title = "Adjusted hazard ratio vs tumor size",
    subtitle = paste0("Reference tumor size = median tumor size (", round(size_ref, 1), " mm)"),
    x = "Tumor size (mm)",
    y = "Hazard ratio"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(figures_dir, "Figure3_Spline_TumorSize.png"), g_size, width = 7, height = 4.8, dpi = 300)
ggsave(file.path(figures_dir, "Figure3_Spline_TumorSize.pdf"), g_size, width = 7, height = 4.8)

# Figure 4: KM by predicted 10-year risk quintile
df_risk <- df %>%
  mutate(
    risk_10y = as.numeric(pred_risk[["120"]]),
    risk_group = dplyr::ntile(risk_10y, 5),
    risk_group = factor(
      risk_group,
      levels = 1:5,
      labels = c("Q1 (lowest risk)", "Q2", "Q3", "Q4", "Q5 (highest risk)")
    )
  )

km_fit <- survival::survfit(survival::Surv(time_months, event) ~ risk_group, data = df_risk)

km_df <- tibble(
  time = km_fit$time,
  surv = km_fit$surv,
  strata = rep(names(km_fit$strata), km_fit$strata)
) %>%
  mutate(risk_group = stringr::str_replace(strata, "^risk_group=", ""))

g_km <- ggplot(km_df, aes(x = time, y = surv, color = risk_group)) +
  geom_step(linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Observed overall survival by CoxPH-predicted risk quintile",
    subtitle = "Risk quintiles defined by predicted mortality risk at 10 years (120 months)",
    x = "Months since diagnosis",
    y = "Overall survival",
    color = "Risk group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(figures_dir, "Figure4_KM_RiskQuintiles_CoxPH.png"), g_km, width = 8.4, height = 5.6, dpi = 300)
ggsave(file.path(figures_dir, "Figure4_KM_RiskQuintiles_CoxPH.pdf"), g_km, width = 8.4, height = 5.6)

message("Done.")
message("Tables written to: ", tables_dir)
message("Figures written to: ", figures_dir)
message("Final analytic cohort: N = ", nrow(df), "; events = ", sum(df$event == 1, na.rm = TRUE))
message("Global apparent C-index: ", round(c_app, 3))
message("Global optimism-corrected C-index: ", round(global_c_corr, 3))