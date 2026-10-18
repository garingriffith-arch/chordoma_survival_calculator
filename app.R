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

theme <- bslib::bs_theme(
  version=5, bootswatch="flatly", primary="#245789",
  bg="#f4f7fb", fg="#253448",
  base_font="Arial, Helvetica, sans-serif",
  heading_font="Arial, Helvetica, sans-serif"
)
styles <- "
:root { --navy:#1f4e79; --ink:#26374b; --muted:#64748b;
        --line:#dfe6ee; --paper:#ffffff; --mist:#f4f7fb; }
body { background:var(--mist); color:var(--ink); }
.app-shell { max-width:1190px; margin:0 auto; padding:27px 14px 48px; }
.hero { background:#fff; border:1px solid var(--line);
  border-radius:22px; padding:26px 30px; margin-bottom:20px; }
.eyebrow { color:#527295; text-transform:uppercase; letter-spacing:.12em;
  font-size:.78rem; font-weight:750; margin-bottom:9px; }
.hero h1 { font-size:clamp(1.65rem,3vw,2.65rem); letter-spacing:-.035em;
  font-weight:750; line-height:1.15; margin:0 0 10px; }
.hero p { color:var(--muted); margin:0; line-height:1.55; max-width:900px; }
.panel { background:var(--paper); border:1px solid var(--line);
  border-radius:20px; padding:23px; margin-bottom:18px; }
.panel h2 { font-size:1.2rem; font-weight:750; margin:0 0 17px; }
.metric { background:#f8fafc; border:1px solid #e2e8f0;
  border-radius:16px; padding:19px; height:100%; }
.metric .period { color:#62768a; font-size:.9rem; margin-bottom:7px; }
.metric .value { color:var(--navy); font-size:2.5rem; line-height:1.15;
  font-weight:760; letter-spacing:-.04em; }
.metric .sub { color:var(--muted); margin-top:6px; font-size:.83rem; }
.track { height:9px; border-radius:99px; background:#dce5ef;
  margin-top:15px; overflow:hidden; }
.fill { height:9px; border-radius:99px; background:#376d9b; }
.explainer { border-top:1px solid var(--line); margin-top:23px;
  padding-top:15px; line-height:1.65; color:#475569; font-size:.92rem; }
.note { font-size:.86rem; color:#64748b; line-height:1.6; }
.form-control, .form-select { min-height:45px; border-radius:11px; }
.btn-primary { border-radius:12px; font-weight:700; min-height:46px; }
footer { color:#718096; font-size:.82rem; padding:8px 2px; line-height:1.6; }
.disclaimer { margin-top:12px; border-left:3px solid #7691ad;
  padding:9px 13px; background:#f7f9fc; font-size:.85rem;
  color:#536579; }
@media(max-width:767px){ .app-shell{padding:14px 8px 32px;}
  .hero,.panel{padding:19px;} }
"

metric_card <- function(label, p, sub) {
  tags$div(class="metric",
    tags$div(class="period",label),
    tags$div(class="value",sprintf("%.1f%%",100*p)),
    tags$div(class="sub",sub),
    tags$div(class="track",
      tags$div(class="fill", style=sprintf("width:%.3f%%;",100*p))))
}

ui <- bslib::page_fluid(
  theme=theme,
  tags$head(tags$style(HTML(styles))),
  tags$main(class="app-shell",
    tags$header(class="hero",
      tags$div(class="eyebrow","Research prediction model · 2024 NCDB"),
      tags$h1("Chordoma survival estimator"),
      tags$p("Exploratory five- and ten-year estimates of all-cause mortality for ",
             "adults with chordoma coded to bones of the skull and face (C41.0). ",
             "Internally validated; not intended for clinical decision-making.")),
    bslib::layout_columns(col_widths=c(4,8),
      tags$section(class="panel",
        tags$h2("Patient characteristics"),
        numericInput("age","Age at diagnosis (years)",value=NA,min=18,max=110,step=1),
        numericInput("size","Recorded tumor diameter (mm)",value=NA,min=1,max=100,step=1),
        selectInput("sex","Sex",choices=c("Select sex"="","Male"="Male","Female"="Female"),selected=""),
        selectInput("cdcc","Charlson–Deyo comorbidity category",
          choices=c("Select category"="","0 (none)"="0","1"="1","2 or more"="2+"),selected=""),
        actionButton("calculate","Calculate estimate",class="btn-primary",width="100%"),
        tags$div(class="disclaimer",
          "All fields are required. The recorded size must be 1–100 mm. ",
          "Registry size may reflect a surgical measurement; pretreatment timing is not verified.")),
      tags$section(class="panel",
        tags$h2("Estimated mortality"),
        uiOutput("risk_cards"),
        tags$div(class="explainer",
          tags$b("How to interpret these estimates"),
          tags$p("Mortality is death from any cause, not death necessarily ",
            "attributable to chordoma. Estimates reflect historical NCDB associations ",
            "and do not predict local recurrence, treatment benefit, or individual outcomes."),
          tags$p("The cohort's C41.0 code does not confirm clival or intracranial ",
            "disease. No independent external validation or prospective impact study ",
            "has been completed.")))),
    tags$footer(
      tags$p("Development cohort: 1,243 adults; 310 deaths. ",
        "Primary model: Cox proportional hazards, age, sex, recorded diameter, ",
        "and categorical Charlson–Deyo score. Estimates at 60 and 120 months."),
      tags$p(
        tags$a("Model equation and analysis code",
          href="https://github.com/garingriffith-arch/chordoma_survival_calculator",
          target="_blank",rel="noopener noreferrer"),
        " · National Cancer Database PUF data are not redistributed.")
    )
  )
)
server <- function(input, output, session) {
  prediction <- eventReactive(input$calculate, {
    validate(
      need(is.numeric(input$age) && length(input$age)==1L &&
             is.finite(input$age) && input$age>=18 && input$age<=110,
           "Enter a valid age between 18 and 110 years."),
      need(is.numeric(input$size) && length(input$size)==1L &&
             is.finite(input$size) && input$size>=1 && input$size<=100,
           "Enter a recorded tumor size between 1 and 100 mm."),
      need(input$sex %in% c("Male","Female"),"Select sex."),
      need(input$cdcc %in% c("0","1","2+"),"Select the comorbidity category.")
    )
    risk_for_profile(input$age,input$size,input$sex,input$cdcc)
  })
  output$risk_cards <- renderUI({
    p <- prediction()
    tagList(
      bslib::layout_columns(col_widths=c(6,6),
        metric_card("At 5 years",p[1],"Estimated all-cause mortality"),
        metric_card("At 10 years",p[2],"Estimated all-cause mortality")),
      tags$p(class="note",style="margin-top:18px",
        sprintf("Equivalent estimated overall survival: %.1f%% at 5 years and %.1f%% at 10 years.",
          100*(1-p[1]),100*(1-p[2])))
    )
  })
}
shiny::shinyApp(ui, server)
