suppressPackageStartupMessages({
  library(shiny)
  library(survival)
  library(bslib)
  library(ggplot2)
  library(rms)
  library(renv)
})

obj <- readRDS(file.path("..", "data", "processed", "chordoma_model_objects.rds"))
fit <- obj$cph_fit
df_ref <- obj$df2

safe_factor <- function(val, levels) {
  if (is.null(levels) || length(levels) == 0) {
    return(factor(val))
  }
  if (!val %in% levels) {
    val <- levels[1]
  }
  factor(val, levels = levels)
}

surv_at <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)$surv
  as.numeric(s[1])
}

# --- NEW: median + 95% CI extraction (months) ---
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
  
  # Fallback: first time S(t) <= 0.5
  if (!is.null(sf$time) && !is.null(sf$surv) && length(sf$time) > 0) {
    idx <- which(sf$surv <= 0.5)[1]
    if (!is.na(idx)) out$med <- as.numeric(sf$time[idx])
  }
  
  out
}

ref_levels <- list(
  SEX3 = levels(df_ref$SEX3),
  RACE3 = levels(df_ref$RACE3),
  HISP2 = levels(df_ref$HISP2),
  INS3 = levels(df_ref$INS3),
  MED_INC_QUAR = levels(df_ref$MED_INC_QUAR)
)

age_default <- round(median(df_ref$AGE_NUM, na.rm = TRUE))
tumor_default <- round(median(df_ref$tumor_size_mm, na.rm = TRUE))
cdcc_default <- round(median(df_ref$cdcc, na.rm = TRUE))

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
      .app-container {
        max-width: 1400px;
        margin: 0 auto;
        padding: 18px 18px 28px 18px;
      }

      .app-header {
        background: #ffffff;
        border-radius: 22px;
        padding: clamp(16px, 2vw, 26px);
        margin-bottom: 20px;
        box-shadow: 0 6px 24px rgba(31, 52, 73, 0.08);
      }

      .header-grid {
        display: grid;
        grid-template-columns: minmax(70px, 10vw) 1fr;
        gap: clamp(14px, 2vw, 26px);
        align-items: center;
      }

      .logo-wrap {
        display: flex;
        align-items: center;
        justify-content: center;
      }

      .ohsu-logo {
        width: clamp(58px, 7vw, 110px);
        height: auto;
        display: block;
      }

      .header-title {
        margin: 0 0 6px 0;
        font-weight: 800;
        line-height: 1.08;
        font-size: clamp(1.9rem, 3.8vw, 3.7rem);
        color: #243447;
      }

      .ohsu-subtitle {
        color: #5b6b7f;
        margin: 0;
        font-size: clamp(0.98rem, 1.4vw, 1.15rem);
      }

      .ohsu-dept {
        color: #738396;
        margin: 0;
        font-size: clamp(0.92rem, 1.2vw, 1.02rem);
      }

      .input-card, .metric-card, .plot-card, .info-card {
        background: #ffffff;
        border: 1px solid #e7edf5 !important;
        border-radius: 22px !important;
        box-shadow: 0 6px 24px rgba(31, 52, 73, 0.06);
      }

      .section-title {
        font-weight: 700;
        color: #243447;
        margin-bottom: 14px;
        line-height: 1.05;
      }

      .metric-card {
        min-height: 138px;
      }

      .metric-value {
        font-size: clamp(2rem, 2.8vw, 2.6rem);
        line-height: 1;
        font-weight: 800;
        color: #1f4e79;
        margin-bottom: 10px;
      }

      .metric-label {
        font-size: 0.98rem;
        color: #5b6b7f;
      }

      .plot-title {
        font-weight: 700;
        color: #243447;
        margin-bottom: 12px;
      }

      .form-label {
        font-weight: 600;
        color: #2f4257;
        margin-bottom: 6px;
      }

      .shiny-input-container {
        margin-bottom: 14px;
      }

      .form-control, .form-select {
        border-radius: 12px !important;
        border: 1px solid #d4dde8 !important;
        min-height: 46px;
      }

      .btn-primary {
        background-color: #245789 !important;
        border-color: #245789 !important;
        border-radius: 14px !important;
        font-weight: 700;
        min-height: 46px;
      }

      .btn-primary:hover {
        background-color: #1d476f !important;
        border-color: #1d476f !important;
      }

      .info-text {
        font-size: 0.98rem;
        color: #425466;
        margin-bottom: 0;
      }

      .info-label {
        font-weight: 700;
        color: #243447;
      }

      .checkbox {
        margin-top: 4px;
      }

      @media (max-width: 768px) {
        .section-title {
          font-size: 1.8rem;
        }

        .metric-card {
          min-height: 118px;
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
          img(
            src = "ohsu_logo.png",
            class = "ohsu-logo"
          )
        ),
        div(
          h1(
            "Intracranial Chordoma Overall Survival Estimator",
            class = "header-title"
          ),
          p("Oregon Health & Science University", class = "ohsu-subtitle"),
          p("Department of Neurological Surgery", class = "ohsu-dept")
        )
      )
    ),
    
    layout_columns(
      col_widths = c(4, 8),
      
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
            choices = ref_levels$INS3,
            selected = ref_levels$INS3[1],
            selectize = FALSE
          ),
          
          selectInput(
            "income", "Median income quartile",
            choices = ref_levels$MED_INC_QUAR,
            selected = ref_levels$MED_INC_QUAR[1],
            selectize = FALSE
          ),
          
          selectInput(
            "cdcc", "Charlson-Deyo score",
            choices = 0:3,
            selected = cdcc_default,
            selectize = FALSE
          ),
          
          div(style = "margin-top: 10px;"),
          
          actionButton(
            "calc",
            "Estimate survival",
            class = "btn-primary w-100"
          ),
          
          checkboxInput("show_ci", "Show confidence bands", value = FALSE)
        )
      ),
      
      div(
        # --- UPDATED: 3 metric cards (5y, 10y, median) ---
        layout_columns(
          col_widths = c(4, 4, 4),
          
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
        
        div(style = "height: 16px;"),
        
        card(
          class = "plot-card",
          card_body(
            h2("Estimated overall survival", class = "plot-title"),
            plotOutput("survplot", height = "450px")
          )
        ),
        
        div(style = "height: 16px;"),
        
        card(
          class = "info-card",
          card_body(
            p(
              class = "info-text",
              span("Intended use: ", class = "info-label"),
              "This tool provides population-level survival estimates derived from the National Cancer Database. It is intended to support clinician-patient discussion and does not replace individualized clinical judgment."
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
  
  # --- NEW: Big text "X ± Y months" (or "X months" / "Not reached") ---
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
        plot.margin = margin(12, 12, 10, 10)
      )
  }, res = 120)
}

shinyApp(ui, server)
