# Investigational 2024 NCDB C41.0 chordoma survival calculator.
# Deployment requires only this source and the two aggregate model CSV files.
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
})

find_root <- function() {
  choices <- c(".", "..")
  valid <- vapply(choices, function(d)
    file.exists(file.path(d, "model", "coefficients.csv")) &&
    file.exists(file.path(d, "model", "reference_survival.csv")),
    logical(1))
  if (!any(valid)) stop("Missing model/coefficients.csv or model/reference_survival.csv")
  normalizePath(choices[which(valid)[1]], winslash = "/", mustWork = TRUE)
}
root <- find_root()
coeff <- read.csv(file.path(root, "model", "coefficients.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
base <- read.csv(file.path(root, "model", "reference_survival.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)

expected <- c("AGE_NUM", "tumor_size_10mm", "SEX3Female", "CDCC21", "CDCC22+")
if (!identical(as.character(coeff$term), expected) ||
    !all(c("beta","reference_design_value") %in% names(coeff)) ||
    any(!is.finite(coeff$beta)) ||
    any(!is.finite(coeff$reference_design_value))) {
  stop("Model coefficients do not match the locked 2024 equation.")
}
if (!identical(as.integer(base$horizon_months), c(60L,120L)) ||
    any(!is.finite(base$reference_survival)) ||
    any(base$reference_survival <= 0 | base$reference_survival >= 1)) {
  stop("Reference survival file does not match 60/120 month model.")
}
if (!isTRUE(all.equal(coeff$reference_design_value,c(52,3.2,0,0,0),
                       tolerance=1e-8))) {
  stop("Unexpected locked equation reference profile.")
}

risk_for_profile <- function(age, size_mm, sex, cdcc) {
  stopifnot(length(age)==1, is.finite(age), age>=18, age<=110,
    length(size_mm)==1,is.finite(size_mm),size_mm>=1,size_mm<=100,
    length(sex)==1,sex %in% c("Male","Female"),
    length(cdcc)==1,cdcc %in% c("0","1","2+"))
  x <- c(age, size_mm/10, as.integer(sex=="Female"),
         as.integer(cdcc=="1"),as.integer(cdcc=="2+"))
  eta <- sum(coeff$beta*(x-coeff$reference_design_value))
  # Equivalent to 1 - S0^exp(eta), evaluated stably for small risks.
  mortality <- -expm1(log(base$reference_survival)*exp(eta))
  pmax(0,pmin(1,mortality))
}


# Presentation follows the oligodendroglioma survival estimator.
# Risk calculation above remains the published 2024 chordoma equation.
horizons <- c(60L, 120L)
repo_url <- "https://github.com/garingriffith-arch/chordoma_survival_calculator"
oligo_logo_url <- paste0(
  "https://raw.githubusercontent.com/garingriffith-arch/",
  "oligodendrioglioma_survival_calculator/main/shiny/www/ohsu_logo.png"
)

fmt_pct <- function(p) sprintf("%.1f%%", 100 * p)

theme <- bslib::bs_theme(
  version = 5, bootswatch = "flatly",
  base_font = bslib::font_google("Inter"),
  heading_font = bslib::font_google("Inter"),
  primary = "#1f4e79",
  bg = "#f4f7fb", fg = "#243447"
)

styles <- "
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
body { background: var(--bg-soft); }
.app-container { max-width: var(--page-max); margin: 0 auto; padding: 24px 22px 36px; }
.app-header {
  background: #fff; border-radius: 28px; padding: clamp(18px, 2.2vw, 30px);
  margin-bottom: 18px; box-shadow: var(--shadow-soft); border: 1px solid var(--border-soft);
}
.header-grid { display:grid; grid-template-columns:minmax(70px,96px) 1fr;gap:20px;align-items:center; }
.logo-wrap { display:flex;align-items:center;justify-content:center; }
.ohsu-logo { width:clamp(58px,6vw,92px);height:auto;display:block; }
.ohsu-logo-fallback {
  display:none; width:88px;height:88px;border-radius:22px;align-items:center;
  justify-content:center;background:#1f4e79;color:#fff;font-weight:900;letter-spacing:.06em;
}
.header-title {
  margin:0 0 8px;font-weight:800;line-height:1.04;
  font-size:clamp(1.9rem,3.4vw,3.2rem);color:var(--text-main);max-width:1000px;
}
.ohsu-subtitle { color:var(--text-muted); margin:0 0 3px; font-size:1.05rem; }
.ohsu-dept { color:#738396; margin:0; font-size:.98rem; }
.input-card,.metric-card,.plot-card,.detail-card {
  background:#fff; border:1px solid var(--border-soft) !important;
  border-radius:var(--card-radius) !important;box-shadow:var(--shadow-soft);
}
.metric-card .card-body,.plot-card .card-body,.detail-card .card-body {padding:22px;}
.input-card .card-body {padding:18px 18px 16px;}
.sticky-panel {
  position:sticky;top:24px;max-height:calc(100vh - 48px);overflow-y:auto;
  padding-right:4px;scrollbar-width:thin;
}
.section-title {
  font-weight:800;color:var(--text-main);margin-bottom:14px;line-height:1.06;
  font-size:clamp(1.55rem,2vw,2rem);
}
.plot-title {font-weight:800;color:var(--text-main);margin-bottom:10px;font-size:1.15rem;}
.form-label {font-weight:650;color:#2f4257;margin-bottom:5px;font-size:.97rem;}
.shiny-input-container {margin-bottom:12px;}
.form-control,.form-select {
  border-radius:14px !important;border:1px solid #d4dde8 !important;
  min-height:44px;box-shadow:none !important;
}
.btn-primary {
  background-color:#245789 !important;border-color:#245789 !important;
  border-radius:14px !important;font-weight:750;min-height:46px;margin-top:6px;
}
.metric-grid {
  display:grid;grid-template-columns:repeat(2,minmax(0,1fr));
  gap:14px;margin-bottom:18px;
}
.metric-card {min-height:118px;}
.metric-value {
  font-size:clamp(1.45rem,2vw,2.05rem);line-height:1;font-weight:800;
  color:var(--accent);margin-bottom:10px;
}
.metric-label {font-size:.92rem;color:var(--text-muted);line-height:1.35;}
.detail-card h3 {
  font-size:1.08rem;font-weight:750;color:var(--text-main);
  margin-top:0;margin-bottom:.8rem;
}
.detail-card ul {margin-bottom:0;padding-left:1.15rem;}
.detail-card li {color:#425466;margin-bottom:.48rem;line-height:1.5;}
.detail-grid {display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:26px 38px;}
.block-gap {height:18px;}
.foot-link {margin-top:18px;color:#607184;font-size:.88rem;}
.foot-link a {color:#245789;}
@media(max-width:1199px) {
  .sticky-panel {position:static;max-height:none;overflow-y:visible;padding-right:0;}
}
@media(max-width:767px) {
  .app-container {padding:18px 14px 28px;}
  .header-grid,.metric-grid,.detail-grid {grid-template-columns:1fr;}
  .header-grid {text-align:center;}
  .plot-card .shiny-plot-output {height:420px !important;}
}
"

logo_ui <- tags$div(
  tags$img(
    src = oligo_logo_url,
    class = "ohsu-logo",
    alt = "OHSU",
    onerror = "this.style.display='none';this.nextElementSibling.style.display='flex';"
  ),
  tags$div(class = "ohsu-logo-fallback", "OHSU")
)

metric_card <- function(h) {
  bslib::card(
    class = "metric-card",
    bslib::card_body(
      tags$div(textOutput(paste0("surv_", h)), class = "metric-value"),
      tags$div(paste0(h / 12, "-year overall survival"), class = "metric-label")
    )
  )
}

ui <- bslib::page_fluid(
  theme = theme,
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(shiny::HTML(styles))
  ),
  tags$main(
    class = "app-container",
    tags$div(
      class = "app-header",
      tags$div(
        class = "header-grid",
        tags$div(class = "logo-wrap", logo_ui),
        tags$div(
          tags$h1("Chordoma Overall Survival Estimator", class = "header-title"),
          tags$p("Oregon Health & Science University", class = "ohsu-subtitle"),
          tags$p("Department of Neurological Surgery", class = "ohsu-dept")
        )
      )
    ),

    bslib::layout_columns(
      col_widths = c(4, 8),

      tags$div(
        class = "sticky-panel",
        bslib::card(
          class = "input-card",
          bslib::card_body(
            tags$h2("Patient and tumor characteristics", class = "section-title"),
            numericInput("age", "Age at diagnosis (years)",
                         value = NA_real_, min = 18, max = 110, step = 1),
            selectInput("sex", "Sex",
                        choices = c("Select sex" = "", "Female" = "Female",
                                    "Male" = "Male"),
                        selected = "", selectize = FALSE),
            selectInput("cdcc", "Charlson-Deyo comorbidity score",
                        choices = c("Select score" = "", "0" = "0", "1" = "1",
                                    "2 or more" = "2+"),
                        selected = "", selectize = FALSE),
            numericInput("size", "Recorded tumor diameter (mm)",
                         value = NA_real_, min = 1, max = 100, step = 1),
            actionButton("calc", "Estimate overall survival",
                         class = "btn-primary w-100")
          )
        )
      ),

      tags$div(
        tags$div(class = "metric-grid", lapply(horizons, metric_card)),
        bslib::card(
          class = "plot-card",
          bslib::card_body(
            tags$h2("Estimated overall survival at 5 and 10 years",
                    class = "plot-title"),
            plotOutput("survplot", height = "560px")
          )
        )
      )
    ),

    tags$div(class = "block-gap"),

    bslib::card(
      class = "detail-card",
      bslib::card_body(
        tags$h2("Model details, validation, and intended use",
                class = "section-title"),
        tags$div(
          class = "detail-grid",
          tags$div(
            tags$h3("Cohort and intended use"),
            tags$ul(
              tags$li("The primary development cohort included 1,243 adults with 310 recorded deaths."),
              tags$li("Predictions are for 5- and 10-year overall survival measured from diagnosis. The endpoint is death from any cause, not necessarily chordoma-related mortality."),
              tags$li("This is a research-stage prognostic model. Clinical prognostic counseling or routine clinical use would require independent external validation, and clinical utility should be evaluated.")
            )
          ),
          tags$div(
            tags$h3("Data and predictors"),
            tags$ul(
              tags$li("Data source: 2024 National Cancer Database Bone/Joint Participant User File; diagnosis years 2004–2023."),
              tags$li("Study population: adults with chordoma recorded at primary site C41.0 (bones of the skull and face). This code does not establish a clival, intracranial, or other precise skull-base location."),
              tags$li("Predictors: age at diagnosis, sex, recorded tumor diameter, and Charlson-Deyo comorbidity category."),
              tags$li("Primary model eligibility required a recorded tumor diameter of 1–100 mm. All calculator inputs are mandatory; absent inputs are not automatically imputed.")
            )
          ),
          tags$div(
            tags$h3("Internal validation"),
            tags$ul(
              tags$li("The primary Cox model was developed using complete cases and internally validated by bootstrap resampling. The Shiny calculator evaluates the locked coefficients and reference survival without refitting."),
              tags$li("Optimism-corrected Harrell C-index: 0.743 (95% CI 0.717–0.775)."),
              tags$li("Optimism-corrected 5- and 10-year time-dependent AUCs: 0.760 and 0.760, respectively."),
              tags$li("Corresponding IPCW Brier scores: 0.128 and 0.177. Independent external validation has not been performed.")
            )
          ),
          tags$div(
            tags$h3("Interpretation and limitations"),
            tags$ul(
              tags$li("The model describes observational registry associations; it does not predict recurrence, treatment response, or the causal survival benefit of any intervention."),
              tags$li("Registry-recorded tumor diameter is not confirmed to have been measured before treatment or surgery. Results should not be extrapolated to recorded diameters outside 1–100 mm."),
              tags$li("The line connecting the displayed 0-, 5-, and 10-year survival points is a visual guide only; no intermediate-year survival probabilities are supplied by this calculator."),
              tags$li("Potential selection bias from complete-case analysis, differences between hospitals, and changes in treatment over time may affect transportability. The model has not undergone external validation or a prospective clinical-impact evaluation.")
            )
          )
        ),
        tags$p(
          class = "foot-link",
          tags$a("Model equation, analysis code, and data-access information",
                 href = repo_url, target = "_blank", rel = "noopener noreferrer"),
          " · Licensed NCDB participant-level data are not publicly distributed."
        )
      )
    )
  )
)

plot_survival <- function(s) {
  # Only 5- and 10-year model predictions are available.
  # The connecting segments are display guides, not an estimated continuous curve.
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mar = c(5, 6, 1, 2), bg = "white", family = "sans",
      cex.axis = 1.05, cex.lab = 1.1)
  plot(NA_real_, NA_real_, type = "n",
       xlim = c(0, 120), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  abline(h = seq(0, 1, 0.2), col = "#e7edf5", lwd = 1)
  abline(v = c(0, 60, 120), col = "#f0f4f9", lwd = 1)
  axis(1, at = c(0, 60, 120), labels = c("0", "60", "120"),
       tick = FALSE, col.axis = "#425466", line = 0.8)
  axis(2, at = seq(0, 1, 0.2),
       labels = paste0(seq(0, 100, 20), "%"),
       las = 1, tick = FALSE, col.axis = "#425466", line = 0.3)
  mtext("Months after diagnosis", side = 1, line = 3.4,
        font = 2, col = "#2f4257")
  mtext("Predicted overall survival", side = 2, line = 4.4,
        font = 2, col = "#2f4257")
  if (is.null(s)) {
    text(60, 0.55,
         "Enter patient characteristics and\nselect Estimate overall survival",
         col = "#718096", cex = 1.07)
  } else {
    points_x <- c(0, 60, 120)
    points_y <- c(1, unname(s))
    lines(points_x, points_y, lwd = 3.3, col = "#1f6feb")
    points(points_x[-1], points_y[-1], pch = 16, cex = 1.6, col = "#1f6feb")
  }
}

server <- function(input, output, session) {
  last_survival <- shiny::reactiveVal(NULL)

  # An estimate is not carried forward after any input is edited.
  observeEvent(list(input$age, input$sex, input$cdcc, input$size),
               { last_survival(NULL) }, ignoreInit = TRUE)

  observeEvent(input$calc, {
    age <- suppressWarnings(as.numeric(input$age))
    size <- suppressWarnings(as.numeric(input$size))
    sex <- as.character(input$sex)
    cdcc <- as.character(input$cdcc)

    if (length(age) != 1L || !is.finite(age) || age < 18 || age > 110) {
      showNotification("Enter a valid age (18–110 years).", type = "error")
      return(invisible(NULL))
    }
    if (length(size) != 1L || !is.finite(size) || size < 1 || size > 100) {
      showNotification("Enter a recorded tumor diameter between 1 and 100 mm.",
                       type = "error")
      return(invisible(NULL))
    }
    if (length(sex) != 1L || !sex %in% c("Male", "Female")) {
      showNotification("Select sex.", type = "error")
      return(invisible(NULL))
    }
    if (length(cdcc) != 1L || !cdcc %in% c("0", "1", "2+")) {
      showNotification("Select the Charlson-Deyo category.", type = "error")
      return(invisible(NULL))
    }
    predicted_mortality <- risk_for_profile(age, size, sex, cdcc)
    last_survival(1 - predicted_mortality)
  })

  for (h in horizons) {
    local({
      horizon <- h
      output[[paste0("surv_", horizon)]] <- renderText({
        s <- last_survival()
        if (is.null(s)) "—" else fmt_pct(s[which(horizons == horizon)])
      })
    })
  }

  output$survplot <- renderPlot({
    plot_survival(last_survival())
  }, res = 120)
}

shiny::shinyApp(ui, server)
