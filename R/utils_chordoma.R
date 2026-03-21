suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(survival)
  library(rms)
})

RAW_PATH <- "C:/OneDrive/garingriffith/OneDrive/OHSU/Research/NRSG/predictive calculator/chordoma/data/raw/NCDB Bone.csv"
PROCESSED_DIR <- "C:/OneDrive/garingriffith/OneDrive/OHSU/Research/NRSG/predictive calculator/chordoma/data/processed"

TUMOR_SIZE_CAP_MM <- 70
HIST_KEEP <- c("9370/3", "9371/3", "9372/3")
SITE_KEEP <- "C41.0"
MODEL_VARS <- c(
  "time_months", "event", "AGE_NUM", "SEX3", "RACE3", "HISP2",
  "INS3", "MED_INC_QUAR", "cdcc", "tumor_size_mm"
)

ensure_project_dirs <- function() {
  dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PROCESSED_DIR, "tables", "manuscript"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PROCESSED_DIR, "figures", "manuscript"), recursive = TRUE, showWarnings = FALSE)
}

normalize_primary_site <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x[x == ""] <- NA_character_
  x <- str_replace_all(x, "[^A-Z0-9\\.]", "")
  x <- str_replace_all(x, "\\.", "")
  x <- ifelse(!is.na(x) & !str_starts(x, "C") & str_detect(x, "^[0-9]+$"), paste0("C", x), x)
  x <- ifelse(!is.na(x) & nchar(x) == 3 & str_detect(x, "^C\\d{2}$"), paste0(x, "0"), x)
  x <- ifelse(
    !is.na(x) & nchar(x) == 4 & str_detect(x, "^C\\d{3}$"),
    paste0(substr(x, 1, 3), ".", substr(x, 4, 4)),
    x
  )
  x
}

parse_mm_exact <- function(x) {
  x <- str_trim(as.character(x))
  v <- suppressWarnings(as.integer(x))
  ifelse(!is.na(v) & v >= 1L & v <= 988L, as.double(v), NA_real_)
}

read_ncdb_bone <- function(path = RAW_PATH) {
  read_csv(
    path,
    col_types = cols(.default = col_character()),
    na = c("", "NA", "N/A")
  )
}

derive_screen_variables <- function(df) {
  stopifnot("TUMOR_SIZE_SUMMARY_2016" %in% names(df))
  
  df %>%
    mutate(
      HIST_BEH = paste0(str_trim(HISTOLOGY), "/", str_trim(BEHAVIOR)),
      PRIMARY_SITE_STD = normalize_primary_site(PRIMARY_SITE),
      
      YEAR_NUM = suppressWarnings(as.integer(YEAR_OF_DIAGNOSIS)),
      AGE_NUM = suppressWarnings(as.numeric(AGE)),
      
      time_months = suppressWarnings(as.numeric(DX_LASTCONTACT_DEATH_MONTHS)),
      vital_status = suppressWarnings(as.integer(PUF_VITAL_STATUS)),
      event = case_when(
        vital_status == 0L ~ 1L,
        vital_status == 1L ~ 0L,
        TRUE ~ NA_integer_
      ),
      
      ts_old_exact = parse_mm_exact(TUMOR_SIZE),
      ts_2016_exact = parse_mm_exact(TUMOR_SIZE_SUMMARY_2016),
      
      tumor_size_mm_raw = case_when(
        !is.na(YEAR_NUM) & YEAR_NUM >= 2016L ~ coalesce(ts_2016_exact, ts_old_exact),
        !is.na(YEAR_NUM) & YEAR_NUM <= 2015L ~ coalesce(ts_old_exact, ts_2016_exact),
        TRUE ~ coalesce(ts_2016_exact, ts_old_exact)
      ),
      
      tumor_size_mm = if_else(
        !is.na(tumor_size_mm_raw) & tumor_size_mm_raw > TUMOR_SIZE_CAP_MM,
        as.double(TUMOR_SIZE_CAP_MM),
        tumor_size_mm_raw
      ),
      
      tumor_size_was_capped = if_else(
        !is.na(tumor_size_mm_raw) & tumor_size_mm_raw > TUMOR_SIZE_CAP_MM,
        1L, 0L, missing = 0L
      )
    )
}

apply_cohort_filters <- function(df) {
  df %>%
    filter(
      HIST_BEH %in% HIST_KEEP,
      PRIMARY_SITE_STD == SITE_KEEP,
      !is.na(AGE_NUM) & AGE_NUM >= 18,
      !is.na(time_months) & time_months >= 0,
      vital_status %in% c(0L, 1L),
      !is.na(event)
    )
}

recode_model_variables <- function(df) {
  df %>%
    mutate(
      MED_INC_QUAR = coalesce(MED_INC_QUAR_2020, MED_INC_QUAR_2016, MED_INC_QUAR_12, MED_INC_QUAR_00),
      MED_INC_QUAR = ifelse(MED_INC_QUAR %in% c("9", "", NA), NA_character_, MED_INC_QUAR),
      
      cdcc = suppressWarnings(as.numeric(CDCC_TOTAL_BEST)),
      
      SEX3 = case_when(
        SEX == "1" ~ "Male",
        SEX == "2" ~ "Female",
        TRUE ~ NA_character_
      ),
      RACE3 = case_when(
        RACE == "01" ~ "White",
        RACE == "02" ~ "Black",
        RACE == "99" ~ NA_character_,
        TRUE ~ "Other"
      ),
      HISP2 = case_when(
        SPANISH_HISPANIC_ORIGIN == "0" ~ "Non-Hispanic",
        SPANISH_HISPANIC_ORIGIN %in% as.character(1:8) ~ "Hispanic",
        SPANISH_HISPANIC_ORIGIN == "9" ~ NA_character_,
        TRUE ~ NA_character_
      ),
      INS3 = case_when(
        INSURANCE_STATUS == "1" ~ "Private",
        INSURANCE_STATUS == "3" ~ "Medicare",
        INSURANCE_STATUS %in% c("2", "4") ~ "Medicaid/OtherGov",
        INSURANCE_STATUS == "0" ~ "Uninsured",
        INSURANCE_STATUS == "9" ~ NA_character_,
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(
      SEX3 = factor(SEX3, levels = c("Male", "Female")),
      RACE3 = factor(RACE3, levels = c("White", "Black", "Other")),
      HISP2 = factor(HISP2, levels = c("Non-Hispanic", "Hispanic")),
      INS3 = factor(INS3, levels = c("Private", "Medicare", "Medicaid/OtherGov", "Uninsured")),
      MED_INC_QUAR = factor(MED_INC_QUAR, levels = c("1", "2", "3", "4"))
    )
}

make_complete_case_model_df <- function(df) {
  df %>%
    recode_model_variables() %>%
    drop_na(all_of(MODEL_VARS))
}

model_formula <- function() {
  Surv(time_months, event) ~
    rcs(AGE_NUM, 4) +
    rcs(tumor_size_mm, 4) +
    SEX3 + RACE3 + HISP2 + INS3 + MED_INC_QUAR + cdcc
}

set_datadist <- function(df) {
  dd <- rms::datadist(df[, MODEL_VARS, drop = FALSE])
  options(datadist = "dd")
  invisible(dd)
}

fmt_n_pct <- function(x) sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))

fmt_med_iqr <- function(x) {
  q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
  sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}

fmt_mean_sd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

make_reference_row <- function(df) {
  out <- data.frame(
    AGE_NUM = median(df$AGE_NUM, na.rm = TRUE),
    tumor_size_mm = median(df$tumor_size_mm, na.rm = TRUE),
    SEX3 = levels(df$SEX3)[1],
    RACE3 = levels(df$RACE3)[1],
    HISP2 = levels(df$HISP2)[1],
    INS3 = levels(df$INS3)[1],
    MED_INC_QUAR = levels(df$MED_INC_QUAR)[1],
    cdcc = median(df$cdcc, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  out$SEX3 <- factor(out$SEX3, levels = levels(df$SEX3))
  out$RACE3 <- factor(out$RACE3, levels = levels(df$RACE3))
  out$HISP2 <- factor(out$HISP2, levels = levels(df$HISP2))
  out$INS3 <- factor(out$INS3, levels = levels(df$INS3))
  out$MED_INC_QUAR <- factor(out$MED_INC_QUAR, levels = levels(df$MED_INC_QUAR))
  out
}

clean_term_label <- function(term, df) {
  if (grepl("^SEX3", term)) {
    lvl <- sub("^SEX3", "", term)
    return(paste0(lvl, " vs ", levels(df$SEX3)[1]))
  }
  if (grepl("^RACE3", term)) {
    lvl <- sub("^RACE3", "", term)
    return(paste0(lvl, " vs ", levels(df$RACE3)[1]))
  }
  if (grepl("^HISP2", term)) {
    lvl <- sub("^HISP2", "", term)
    return(paste0(lvl, " vs ", levels(df$HISP2)[1]))
  }
  if (grepl("^INS3", term)) {
    lvl <- sub("^INS3", "", term)
    return(paste0(lvl, " vs ", levels(df$INS3)[1]))
  }
  if (grepl("^MED_INC_QUAR", term)) {
    lvl <- sub("^MED_INC_QUAR", "", term)
    return(paste0("Income quartile ", lvl, " vs quartile ", levels(df$MED_INC_QUAR)[1]))
  }
  if (term == "cdcc") return("CDCC (per 1-point increase)")
  term
}

S_at <- function(sf, t) {
  ss <- summary(sf, times = t, extend = TRUE)
  if (length(ss$surv) == 0) return(NA_real_)
  as.numeric(ss$surv[length(ss$surv)])
}

calibration_bins <- function(time, event, pred_risk, horizon, groups = 5) {
  tibble(time = time, event = event, p = pred_risk) %>%
    filter(is.finite(p)) %>%
    mutate(bin = ntile(p, groups)) %>%
    group_by(bin) %>%
    summarise(
      n = n(),
      pred_mean = mean(p, na.rm = TRUE),
      obs_mean = {
        sf <- survfit(Surv(time, event) ~ 1, data = pick(everything()))
        1 - S_at(sf, horizon)
      },
      .groups = "drop"
    ) %>%
    mutate(horizon_m = horizon)
}

predict_risk_horizon_coxph <- function(fit, newdata, horizon) {
  sf <- survfit(fit, newdata = newdata)
  surv <- summary(sf, times = horizon, extend = TRUE)$surv
  as.numeric(1 - surv)
}

apply_factor_levels <- function(dat, factor_vars, factor_levels) {
  dat <- as.data.frame(dat)
  for (v in factor_vars) {
    dat[[v]] <- factor(as.character(dat[[v]]), levels = factor_levels[[v]])
  }
  dat
}