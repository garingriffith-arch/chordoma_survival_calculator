suppressPackageStartupMessages({
  library(shiny)
  library(survival)
  library(bslib)
  library(ggplot2)
  library(rms)
})

obj <- readRDS(file.path("..", "data", "processed", "chordoma_model_objects.rds"))
fit <- obj$cph_fit
df_ref <- obj$df2

model_n <- nrow(df_ref)

safe_factor <- function(val, levels) {
  if (is.null(levels) || length(levels) == 0) return(factor(val))
  if (!val %in% levels) val <- levels[1]
  factor(val, levels = levels)
}

surv_at <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)$surv
  as.numeric(s[1])
}

median_stats <- function(sf) {
  out <- list(med = NA_real_, lower = NA_real_, upper = NA_real_)
  
  q <- tryCatch(
    quantile(sf, probs = 0.5, conf.int = TRUE),
    error = function(e) NULL
  )
  
  if (!is.null(q)) {
    out$med <- suppressWarnings(as.numeric(q$quantile[1]))
    if (!is.null(q$lower)) out$lower <- suppressWarnings(as.numeric(q$lower[1]))
    if (!is.null(q$upper)) out$upper <- suppressWarnings(as.numeric(q$upper[1]))
    return(out)
  }
  
  if (!is.null(sf$time) && !is.null(sf$surv) && length(sf$time) > 0) {
    idx <- which(sf$surv <= 0.5)[1]
    if (!is.na(idx)) out$med <- as.numeric(sf$time[idx])
  }
  
  out
}

pretty_ins <- function(x) {
  x <- as.character(x)
  x <- gsub("OtherGov", "Other Government", x, fixed = TRUE)
  x <- gsub("/", " / ", x, fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

pretty_income <- function(x) {
  paste("Quartile", as.character(x))
}

ref_levels <- list(
  SEX3 = levels(df_ref$SEX3),
  RACE3 = levels(df_ref$RACE3),
  HISP2 = levels(df_ref$HISP2),
  INS3 = levels(df_ref$INS3),
  MED_INC_QUAR = levels(df_ref$MED_INC_QUAR)
)

ins_choices <- setNames(ref_levels$INS3, pretty_ins(ref_levels$INS3))
inc_choices <- setNames(ref_levels$MED_INC_QUAR, pretty_income(ref_levels$MED_INC_QUAR))

age_default <- round(median(df_ref$AGE_NUM, na.rm = TRUE))
tumor_default <- round(median(df_ref$tumor_size_mm, na.rm = TRUE))
cdcc_default <- round(median(df_ref$cdcc, na.rm = TRUE))

income_num <- suppressWarnings(as.numeric(as.character(df_ref$MED_INC_QUAR)))
income_default <- as.character(round(stats::median(income_num, na.rm = TRUE)))
if (!income_default %in% ref_levels$MED_INC_QUAR) income_default <- ref_levels$MED_INC_QUAR[1]

ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter"),
    primary = "#1f4e79",
    bg = "#f4f7fb",
    fg = "#243447"
  ),
  
  tags$head(
    tags$style(HTML("
      :root {
        --page-max: 1320px;
        --card-radius: 24px;
        --shadow-soft: 0 8px 28px rgba(31, 52, 73, 0.07);
        --border-soft: #e7edf5;
        --text-main: #243447;
        --text-muted: #5b6b7f;
        --bg-soft: #f4f7fb;
        --accent: #1f4e79;
      }

      body {
        background: var(--bg-soft);
      }

      .app-container {
        max-width: var(--page-max);
        margin: 0 auto;
        padding: 24px 22px 36px 22px;
      }

      .app-header {
        background: #ffffff;
        border-radius: 28px;
        padding: clamp(18px, 2.2vw, 30px);
        margin-bottom: 24px;
        box-shadow: var(--shadow-soft);
        border: 1px solid var(--border-soft);
      }

      .header-grid {
        display: grid;
        grid-template-columns: minmax(70px, 96px) 1fr;
        gap: 20px;
        align-items: center;
      }

      .logo-wrap {
        display: flex;
        align-items: center;
        justify-content: center;
      }

      .ohsu-logo {
        width: clamp(58px, 6vw, 92px);
        height: auto;
        display: block;
      }

      .header-title {
        margin: 0 0 8px 0;
        font-weight: 800;
        line-height: 1.04;
        font-size: clamp(2rem, 3.7vw, 3.4rem);
        color: var(--text-main);
        max-width: 900px;
      }

      .ohsu-subtitle {
        color: var(--text-muted);
        margin: 0 0 3px 0;
        font-size: 1.05rem;
      }

      .ohsu-dept {
        color: #738396;
        margin: 0;
        font-size: 0.98rem;
      }

      .input-card, .metric-card, .plot-card, .detail-card {
        background: #ffffff;
        border: 1px solid var(--border-soft) !important;
        border-radius: var(--card-radius) !important;
        box-shadow: var(--shadow-soft);
      }

      .metric-card .card-body,
      .plot-card .card-body,
      .detail-card .card-body {
        padding: 22px;
      }

      .input-card .card-body {
        padding: 18px 18px 16px 18px;
      }

      .sticky-panel {
        position: sticky;
        top: 24px;
      }

      .section-title {
        font-weight: 800;
        color: var(--text-main);
        margin-bottom: 14px;
        line-height: 1.06;
        font-size: clamp(1.55rem, 2vw, 2rem);
      }

      .plot-title {
        font-weight: 800;
        color: var(--text-main);
        margin-bottom: 10px;
        font-size: 1.15rem;
      }

      .form-label {
        font-weight: 650;
        color: #2f4257;
        margin-bottom: 5px;
        font-size: 0.97rem;
      }

      .shiny-input-container {
        margin-bottom: 10px;
      }

      .form-control, .form-select {
        border-radius: 14px !important;
        border: 1px solid #d4dde8 !important;
        min-height: 44px;
        box-shadow: none !important;
      }

      .form-control:focus, .form-select:focus {
        border-color: #8db1d5 !important;
        box-shadow: 0 0 0 0.16rem rgba(31, 78, 121, 0.10) !important;
      }

      .btn-primary {
        background-color: #245789 !important;
        border-color: #245789 !important;
        border-radius: 14px !important;
        font-weight: 750;
        min-height: 46px;
        margin-top: 6px;
      }

      .btn-primary:hover {
        background-color: #1d476f !important;
        border-color: #1d476f !important;
      }

      .metric-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 16px;
        margin-bottom: 18px;
      }

      .metric-card {
        min-height: 126px;
      }

      .metric-value {
        font-size: clamp(1.8rem, 2.4vw, 2.45rem);
        line-height: 1;
        font-weight: 800;
        color: var(--accent);
        margin-bottom: 10px;
      }

      .metric-label {
        font-size: 0.96rem;
        color: var(--text-muted);
        line-height: 1.35;
      }

      .detail-card h3 {
        font-size: 1.08rem;
        font-weight: 750;
        color: var(--text-main);
        margin-top: 0;
        margin-bottom: 0.8rem;
      }

      .detail-card ul {
        margin-bottom: 0;
        padding-left: 1.15rem;
      }

      .detail-card li {
        color: #425466;
        margin-bottom: 0.48rem;
        line-height: 1.5;
      }

      .block-gap {
        height: 18px;
      }

      .checkbox {
        margin-top: 6px;
        margin-bottom: 0;
      }

      .plot-card .shiny-plot-output {
        margin-top: 2px;
      }

      .detail-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 26px 38px;
      }

      .detail-section {
        min-width: 0;
      }

      @media (max-width: 1199px) {
        .metric-grid {
          grid-template-columns: repeat(2, minmax(0, 1fr));
        }
        .sticky-panel {
          position: static;
        }
      }

      @media (max-width: 767px) {
        .app-container {
          padding: 18px 14px 28px 14px;
        }
        .header-grid {
          grid-template-columns: 1fr;
          gap: 14px;
          text-align: center;
        }
        .logo-wrap {
          justify-content: center;
        }
        .metric-grid,
        .detail-grid {
          grid-template-columns: 1fr;
        }
        .input-card .card-body,
        .metric-card .card-body,
        .plot-card .card-body,
        .detail-card .card-body {
          padding: 18px;
        }
      }
    "))
  ),
  
  div(
    class = "app-container",
    
    div(
      class = "app-header",
      div(
        class = "header-grid",
        div(
          class = "logo-wrap",
          img(src = "ohsu_logo.png", class = "ohsu-logo")
        ),
        div(
          h1("Intracranial Chordoma Overall Survival Estimator", class = "header-title"),
          p("Oregon Health & Science University", class = "ohsu-subtitle"),
          p("Department of Neurological Surgery", class = "ohsu-dept")
        )
      )
    ),
    
    layout_columns(
      col_widths = c(4, 8),
      
      div(
        class = "sticky-panel",
        card(
          class = "input-card",
          card_body(
            h2("Patient characteristics", class = "section-title"),
            
            numericInput("age", "Age (years)", value = age_default, min = 18, max = 90, step = 1),
            
            numericInput(
              "tsize_mm",
              "Tumor size (mm)",
              value = tumor_default,
              min = 1,
              max = 70,
              step = 1
            ),
            
            selectInput(
              "sex3", "Sex",
              choices = ref_levels$SEX3,
              selected = ref_levels$SEX3[1],
              selectize = FALSE
            ),
            
            selectInput(
              "race3", "Race",
              choices = ref_levels$RACE3,
              selected = ref_levels$RACE3[1],
              selectize = FALSE
            ),
            
            selectInput(
              "hisp2", "Ethnicity",
              choices = ref_levels$HISP2,
              selected = ref_levels$HISP2[1],
              selectize = FALSE
            ),
            
            selectInput(
              "ins3", "Insurance",
              choices = ins_choices,
              selected = ref_levels$INS3[1],
              selectize = FALSE
            ),
            
            selectInput(
              "income", "Median income quartile",
              choices = inc_choices,
              selected = income_default,
              selectize = FALSE
            ),
            
            selectInput(
              "cdcc", "Charlson-Deyo score",
              choices = 0:3,
              selected = cdcc_default,
              selectize = FALSE
            ),
            
            actionButton("calc", "Estimate survival", class = "btn-primary w-100"),
            checkboxInput("show_ci", "Show confidence bands", value = FALSE)
          )
        )
      ),
      
      div(
        div(
          class = "metric-grid",
          
          card(
            class = "metric-card",
            card_body(
              div(textOutput("s5"), class = "metric-value"),
              div("5-year overall survival", class = "metric-label")
            )
          ),
          
          card(
            class = "metric-card",
            card_body(
              div(textOutput("s10"), class = "metric-value"),
              div("10-year overall survival", class = "metric-label")
            )
          ),
          
          card(
            class = "metric-card",
            card_body(
              div(textOutput("med"), class = "metric-value"),
              div(
                tags$span(
                  "Median predicted survival",
                  class = "metric-label",
                  title = "‘Not reached’ means the predicted survival curve does not fall below 50% within available follow-up, so the median time cannot be estimated."
                )
              )
            )
          )
        ),
        
        div(class = "block-gap"),
        
        card(
          class = "plot-card",
          card_body(
            h2("Estimated overall survival", class = "plot-title"),
            plotOutput("survplot", height = "560px")
          )
        )
      )
    ),
    
    div(class = "block-gap"),
    
    card(
      class = "detail-card",
      card_body(
        h2("Model details, analysis summary, and intended use", class = "section-title"),
        
        div(
          class = "detail-grid",
          
          div(
            class = "detail-section",
            h3("Model cohort and intended use"),
            tags$ul(
              tags$li(paste0("Model cohort: n = ", format(model_n, big.mark = ","), " adults with intracranial chordoma.")),
              tags$li("Intended use: this tool provides diagnosis-time, population-level overall survival estimates derived from the National Cancer Database."),
              tags$li("It is intended to support clinician-patient discussion and risk stratification and does not replace individualized clinical judgment."),
              tags$li("Predictions should be interpreted in the context of imaging, pathology, surgical planning, and multidisciplinary evaluation.")
            )
          ),
          
          div(
            class = "detail-section",
            h3("Cohort and variables"),
            tags$ul(
              tags$li("Data source: National Cancer Database (NCDB) bone tumor extract."),
              tags$li("Study population: adults with intracranial chordoma identified using ICD-O-3 histology/behavior codes 9370/3, 9371/3, and 9372/3 with primary site C41.0."),
              tags$li("Outcome: overall survival, measured in months from diagnosis."),
              tags$li("Predictors included in the model: age, tumor size, sex, race, Hispanic ethnicity, insurance category, median household income quartile, and Charlson-Deyo comorbidity score.")
            )
          ),
          
          div(
            class = "detail-section",
            h3("Statistical analysis"),
            tags$ul(
              tags$li("Model type: multivariable Cox proportional hazards regression."),
              tags$li("Age and tumor size were modeled flexibly using restricted cubic splines in the final model."),
              tags$li("Displayed outputs include predicted overall survival at 5 and 10 years, the estimated survival curve, and median predicted survival when estimable."),
              tags$li("Internal validation was performed with 300 bootstrap resamples, including global and horizon-specific performance assessment.")
            )
          ),
          
          div(
            class = "detail-section",
            h3("Performance and interpretation"),
            tags$ul(
              tags$li("Internal validation showed good discrimination, with an optimism-corrected Harrell C-index of 0.762 and optimism-corrected calibration slope of 0.912."),
              tags$li("Optimism-corrected AUCs were 0.782 at 5 years and 0.798 at 10 years."),
              tags$li("The calculator is based on registry data and does not incorporate extent of resection, recurrence, radiation details, or molecular markers."),
              tags$li("External validation is still needed before broad clinical application.")
            )
          )
        )
      )
    )
  )
)
server <- function(input, output, session) {
  
  observe({
    current_age <- suppressWarnings(as.numeric(input$age))
    if (!is.na(current_age) && current_age > 90) {
      updateNumericInput(session, "age", value = 90)
    }
    if (!is.na(current_age) && current_age < 18) {
      updateNumericInput(session, "age", value = 18)
    }
  })
  
  observe({
    current_val <- suppressWarnings(as.numeric(input$tsize_mm))
    if (!is.na(current_val) && current_val > 70) {
      updateNumericInput(session, "tsize_mm", value = 70)
    }
    if (!is.na(current_val) && current_val < 1) {
      updateNumericInput(session, "tsize_mm", value = 1)
    }
  })
  
  newdata <- eventReactive(input$calc, {
    age_val <- suppressWarnings(as.numeric(input$age))
    if (!is.finite(age_val)) age_val <- age_default
    age_val <- min(max(age_val, 18), 90)
    
    tum <- suppressWarnings(as.numeric(input$tsize_mm))
    if (!is.finite(tum)) tum <- tumor_default
    tum <- min(max(tum, 1), 70)
    
    data.frame(
      AGE_NUM = age_val,
      tumor_size_mm = tum,
      SEX3 = safe_factor(input$sex3, ref_levels$SEX3),
      RACE3 = safe_factor(input$race3, ref_levels$RACE3),
      HISP2 = safe_factor(input$hisp2, ref_levels$HISP2),
      INS3 = safe_factor(input$ins3, ref_levels$INS3),
      MED_INC_QUAR = safe_factor(input$income, ref_levels$MED_INC_QUAR),
      cdcc = as.numeric(input$cdcc),
      check.names = FALSE
    )
  })
  
  surv_obj <- eventReactive(input$calc, {
    req(newdata())
    survfit(fit, newdata = newdata())
  })
  
  output$s5 <- renderText({
    req(surv_obj())
    sprintf("%.1f%%", 100 * surv_at(surv_obj(), 60))
  })
  
  output$s10 <- renderText({
    req(surv_obj())
    sprintf("%.1f%%", 100 * surv_at(surv_obj(), 120))
  })
  
  output$med <- renderText({
    req(surv_obj())
    ms <- median_stats(surv_obj())
    
    if (!is.finite(ms$med)) return("Not reached")
    
    med_round <- as.integer(round(ms$med))
    
    if (is.finite(ms$lower) && is.finite(ms$upper)) {
      pm <- as.integer(round((ms$upper - ms$lower) / 2))
      if (is.finite(pm) && pm > 0) {
        return(sprintf("%d \u00B1 %d months", med_round, pm))
      }
    }
    
    sprintf("%d months", med_round)
  })
  
  output$survplot <- renderPlot({
    req(surv_obj())
    sf <- surv_obj()
    
    df_plot <- data.frame(
      time = sf$time,
      surv = sf$surv
    )
    
    if (!is.null(sf$lower) && !is.null(sf$upper)) {
      df_plot$lower <- sf$lower
      df_plot$upper <- sf$upper
    } else {
      df_plot$lower <- NA_real_
      df_plot$upper <- NA_real_
    }
    
    df_plot <- df_plot[df_plot$time <= 120, , drop = FALSE]
    
    if (nrow(df_plot) == 0) {
      df_plot <- data.frame(
        time = c(0, 120),
        surv = c(1, 1),
        lower = c(1, 1),
        upper = c(1, 1)
      )
    } else if (min(df_plot$time) > 0) {
      df_plot <- rbind(
        data.frame(time = 0, surv = 1, lower = 1, upper = 1),
        df_plot
      )
    }
    
    pts <- c(60, 120)
    s_main <- summary(sf, times = pts, extend = TRUE)
    
    pts_df <- data.frame(
      time = s_main$time,
      surv = s_main$surv
    )
    
    guide_df <- data.frame(
      time = c(60, 120)
    )
    
    ggplot(df_plot, aes(x = time, y = surv)) +
      {
        if (isTRUE(input$show_ci)) {
          geom_ribbon(
            aes(ymin = lower, ymax = upper),
            fill = "#8fb3d9",
            alpha = 0.22
          )
        }
      } +
      geom_vline(
        data = guide_df,
        aes(xintercept = time),
        linetype = "dashed",
        linewidth = 0.55,
        color = "#cbd6e2"
      ) +
      geom_step(
        color = "#1f6feb",
        linewidth = 1.5,
        direction = "hv"
      ) +
      geom_point(
        data = pts_df,
        aes(x = time, y = surv),
        inherit.aes = FALSE,
        color = "#1f6feb",
        size = 3.2
      ) +
      scale_x_continuous(
        limits = c(0, 121),
        breaks = seq(0, 120, by = 12),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      scale_y_continuous(
        limits = c(0, 1.04),
        breaks = seq(0, 1, by = 0.2),
        labels = function(x) sprintf("%.1f", x),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      labs(
        x = "Months",
        y = "Overall survival"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#e7edf5", linewidth = 0.7),
        axis.title = element_text(color = "#2f4257", face = "bold"),
        axis.text = element_text(color = "#425466"),
        plot.background = element_rect(fill = "#ffffff", color = NA),
        panel.background = element_rect(fill = "#ffffff", color = NA),
        plot.margin = margin(10, 10, 8, 8)
      )
  }, res = 120)
}

shinyApp(ui, server)
