# 07. Manuscript exhibits from the completed 2024 NCDB analysis.
# Run after stages 01-06. Outputs remain local and REQUIRE disclosure review.
source(file.path("R", "00_config.R"))

need <- function(path) {
  if (!file.exists(path)) stop("Missing upstream input: ", path, "\nRun its preceding stage first.")
  path
}
load_csv <- function(name) utils::read.csv(need(file.path(AGG_DIR, name)),
  check.names = FALSE, stringsAsFactors = FALSE)
require_fields <- function(x, fields, name) {
  absent <- setdiff(fields, names(x))
  if (length(absent)) stop(name, " missing columns: ", paste(absent, collapse = ", "))
}
models <- readRDS(need(file.path(PRIVATE_DIR, "models_INTERNAL_ONLY.rds")))
cohort <- readRDS(need(file.path(PRIVATE_DIR, "cohort_INTERNAL_ONLY.rds")))
A <- models$datasets$A_expanded
D <- models$datasets$D_common
eligible <- cohort$eligible
if (nrow(A) != 1243L || sum(A$event) != 310L ||
    nrow(D) != 946L || sum(D$event) != 242L ||
    nrow(eligible) != 1562L) stop("Locked cohort mismatch; abort manuscript exports.")
if (any(A$tumor_size_mm_raw < 1 | A$tumor_size_mm_raw > 100) ||
    any(D$tumor_size_mm_raw < 1 | D$tumor_size_mm_raw > 100))
  stop("Primary-model size eligibility was violated.")

coefficients <- load_csv("cox_coefficients.csv")
corrected <- load_csv("performance_optimism_corrected.csv")
intervals <- load_csv("nested_corrected_performance_CI.csv")
paired <- load_csv("nested_paired_D_minus_A_CI.csv")
mi <- load_csv("MI_pooled_coefficients.csv")
sensitivity <- load_csv("tumor_size_sensitivity.csv")
followup <- load_csv("Table_1_followup.csv")
calibration <- load_csv("Table_calibration_in_the_large.csv")
missing <- load_csv("predictor_missingness_LOCAL_DISCLOSURE_REVIEW.csv")
subgroup <- load_csv("subgroup_calibration_INTERNAL_DISCLOSURE_REVIEW.csv")
flow <- load_csv("cohort_flow.csv")
reference_coefficients <- utils::read.csv(need(file.path(ROOT, "model", "coefficients.csv")))
reference_survival <- utils::read.csv(need(file.path(ROOT, "model", "reference_survival.csv")))
need(file.path(AGG_DIR, "equation_check.csv"))
require_fields(coefficients, c("model","term","beta","SE","HR","low95","high95"), "Coefficients")
require_fields(corrected, c("model","metric","apparent","optimism_corrected"), "Validation")
require_fields(intervals, c("model","metric","lower95","upper95"), "Nested bootstrap")
require_fields(paired, c("metric","lower95","upper95"), "Paired intervals")
require_fields(sensitivity, c("model","max_mm","n","deaths","HR_per_10mm","lower95","upper95"), "Size sensitivity")
require_fields(subgroup, c("model","horizon_months","stratifier","group","n","at_risk",
                           "observed_mortality","mean_predicted_mortality"), "Subgroups")
if (!all(c(100,150,988) %in% sensitivity$max_mm)) stop("Incomplete size sensitivity.")
if (!all(c(60,120) %in% reference_survival$horizon_months)) stop("Baseline risk horizons missing.")
if (nrow(corrected) != 15L || nrow(intervals) != 15L) stop("Incomplete model metrics.")

# A single designated source: the successfully rerun public GitHub pipeline.
# Do not splice point estimates from an older bootstrap run into these CIs.
target <- file.path(PRIVATE_DIR, "manuscript_exhibits_DISCLOSURE_REVIEW")
tabdir <- file.path(target, "tables")
figdir <- file.path(target, "figures")
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
save_table <- function(data, name) {
  utils::write.csv(as.data.frame(data), file.path(tabdir, paste0(name, ".csv")),
                   row.names = FALSE, na = "")
}
fmt <- function(x, n = 3L) ifelse(is.finite(x), sprintf(paste0("%.", n, "f"), x), "NE")
fmt_ci <- function(lo, hi) paste0(fmt(lo), " to ", fmt(hi))
fmt_quartiles <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NE")
  z <- stats::quantile(x, c(.25,.5,.75), names=FALSE)
  paste0(fmt(z[2],1), " [",fmt(z[1],1),", ",fmt(z[3],1),"]")
}
read_med <- function(name) {
  d <- followup[followup$model == name, , drop=FALSE]
  if (nrow(d)!=1L) stop("Follow-up row missing: ",name)
  d
}
fA <- read_med("A_expanded"); fD <- read_med("D_common")
t1 <- data.frame(
  Characteristic=c("N","Deaths","Age at diagnosis, years: median [IQR]",
    "Recorded tumor size, mm: median [IQR]","Diagnosis year: median [IQR]",
    "Median reverse-KM follow-up, months","At risk at 60 months",
    "At risk at 120 months","Tumor-size eligibility"),
  Clinical_model=c(nrow(A),sum(A$event),fmt_quartiles(A$AGE_NUM),
    fmt_quartiles(A$tumor_size_mm_raw),fmt_quartiles(A$YEAR_NUM),
    fmt(fA$median_followup_reverse_KM_months,2),fA$at_risk_60,fA$at_risk_120,"1–100 mm"),
  Common_cohort=c(nrow(D),sum(D$event),fmt_quartiles(D$AGE_NUM),
    fmt_quartiles(D$tumor_size_mm_raw),fmt_quartiles(D$YEAR_NUM),
    fmt(fD$median_followup_reverse_KM_months,2),fD$at_risk_60,fD$at_risk_120,"1–100 mm"))
# Add core predictor distributions; withhold a whole block when a small
# category is present in either cohort to prevent reconstruction.
for (v in c("SEX3", "CDCC2", "RACE3", "HISP2", "INS3", "MED_INC_QUAR")) {
  lv <- levels(A[[v]])
  if (!length(lv)) stop("Expected factor variable missing: ", v)
  ct <- function(d) as.integer(table(factor(d[[v]], levels = lv)))
  nA <- ct(A); nD <- ct(D)
  misA <- sum(is.na(A[[v]])); misD <- sum(is.na(D[[v]]))
  # Both nested cohorts are displayed: small differences must also be suppressed.
  safe <- all(nA == 0L | nA >= 11L) &&
          all(nD == 0L | nD >= 11L) &&
          all((nA - nD) == 0L | (nA - nD) >= 11L) &&
          (misA == 0L || misA >= 11L) &&
          (misD == 0L || misD >= 11L) &&
          ((misA-misD) == 0L || (misA-misD) >= 11L)
  if (!safe) {
    t1 <- rbind(t1, data.frame(Characteristic=paste0(v, " (category distribution)"),
      Clinical_model="Withheld: disclosure review", Common_cohort="Withheld: disclosure review"))
  } else for (j in seq_along(lv)) {
    t1 <- rbind(t1, data.frame(Characteristic=paste0(v, ": ",lv[j]),
      Clinical_model=sprintf("%d (%.1f%%)",nA[j],100*nA[j]/nrow(A)),
      Common_cohort=sprintf("%d (%.1f%%)",nD[j],100*nD[j]/nrow(D))))
  }
}
save_table(t1,"Table_1_Cohort")

terms <- function(x) {
  labels <- c(AGE_NUM="Age, per year",tumor_size_10mm="Tumor size, per 10 mm",
     SEX3Female="Female vs male",CDCC21="Charlson-Deyo 1 vs 0",
     "CDCC22+"="Charlson-Deyo >=2 vs 0",
     RACE3Black="Black vs White",RACE3Other="Other race vs White",
     HISP2Hispanic="Hispanic vs non-Hispanic",
     INS3Medicare="Medicare vs private insurance",
     "INS3Medicaid/OtherGov"="Medicaid/other government vs private",
     INS3Uninsured="Uninsured vs private",
     MED_INC_QUAR2="Income quartile 2 vs 1",
     MED_INC_QUAR3="Income quartile 3 vs 1",
     MED_INC_QUAR4="Income quartile 4 vs 1")
  ifelse(x %in% names(labels),unname(labels[x]),x)
}
cox_tab <- function(model) {
  d <- coefficients[coefficients$model == model,,drop=FALSE]
  if (!nrow(d)) stop("Missing model coefficients: ",model)
  data.frame(Predictor=terms(d$term),Encoded_term=d$term,
    Beta=fmt(d$beta,6),SE=fmt(d$SE,6),HR=fmt(d$HR),
    Wald_95pct_CI=fmt_ci(d$low95,d$high95))
}
save_table(cox_tab("A_expanded"),"Table_2_Primary_Cox_Model")
save_table(cox_tab("D_common"),"Table_S2_Full_Comparator_Coefficients")

validated <- merge(corrected,intervals,by=c("model","metric"),all.x=TRUE,sort=FALSE)
if (nrow(validated)!=15L || anyNA(validated$lower95) ||
    anyNA(validated$upper95)) stop("Cannot match all corrected estimates to nested intervals.")
metric_labels <- c(Harrell_C="Harrell C-index",AUC_60="5-year AUC",
  AUC_120="10-year AUC",Brier_60="5-year IPCW Brier",
  Brier_120="10-year IPCW Brier")
validated$Metric <- unname(metric_labels[validated$metric])
save_table(data.frame(Model=validated$model,Metric=validated$Metric,
  Apparent=fmt(validated$apparent),
  Optimism_corrected=fmt(validated$optimism_corrected),
  Nested_bootstrap_95pct_CI=fmt_ci(validated$lower95,validated$upper95)),
  "Table_3_Internal_Validation")

pc <- corrected[corrected$model=="A_common",c("metric","optimism_corrected")]
pd <- corrected[corrected$model=="D_common",c("metric","optimism_corrected")]
names(pc)[2] <- "A_corrected";names(pd)[2] <- "D_corrected"
comparison <- merge(merge(pc,pd,by="metric"),paired,by="metric")
if(nrow(comparison)!=5L) stop("Expected five common-cohort performance contrasts.")
save_table(data.frame(Metric=unname(metric_labels[comparison$metric]),
  Clinical_corrected=fmt(comparison$A_corrected),
  Full_corrected=fmt(comparison$D_corrected),
  Full_minus_clinical=fmt(comparison$D_corrected-comparison$A_corrected),
  Paired_nested_95pct_CI=fmt_ci(comparison$lower95,comparison$upper95),
  Note="Same 946 patients; internally validated; percentile CI"),
  "Table_4_Common_Cohort_Comparison")

# Exclusion and missingness are descriptive, not mutually exclusive reasons.
require_fields(missing,c("variable","missing","denominator","missing_percent"),"Missingness")
save_table(data.frame(Predictor=missing$variable,
  Missing=ifelse(missing$missing > 0 & missing$missing < 11,
                 "Withheld: small cell",as.character(missing$missing)),
  Missing_percent=ifelse(missing$missing > 0 & missing$missing < 11,
                         "Withheld",fmt(missing$missing_percent,1)),
  Denominator=missing$denominator,
  Note="Predictor missingness overlaps; disclosure review required"),
  "Table_S1_Missingness")

require_fields(mi,c("model","term","estimate","conf.low","conf.high"),"MI")
miA <- mi[mi$model=="A_clinical",,drop=FALSE]
if(nrow(miA)!=5L) stop("Unexpected MI clinical-model terms.")
save_table(data.frame(Predictor=terms(miA$term),Encoded_term=miA$term,
  Pooled_HR=fmt(exp(miA$estimate)),
  Pooled_95pct_CI=fmt_ci(exp(miA$conf.low),exp(miA$conf.high)),
  Note="Outcome-informed MI: coefficient sensitivity only"),
  "Table_S3_Imputation_Sensitivity")

# Entire subgroup blocks are withheld unless EVERY category is eligible.
sublist <- list()
for (h in HORIZONS) for (v in c("HISP2","RACE3","INS3","SEX3")) {
  d <- subgroup[subgroup$model=="A_expanded" & subgroup$horizon_months==h &
                subgroup$stratifier==v,,drop=FALSE]
  levels_expected <- levels(A[[v]])
  if(length(levels_expected)<2L || nrow(d)!=length(levels_expected) ||
     !setequal(as.character(d$group),levels_expected) ||
     any(d$n<11L|d$at_risk<11L)) next
  # Unrecorded race/ethnicity/insurance does not make the model cohort smaller.
  # Permit a missing-category complement only if it is zero or at least 11.
  unclassified <- nrow(A)-sum(d$n)
  if(unclassified < 0L || (unclassified > 0L && unclassified < 11L)) next
  sublist[[length(sublist)+1L]] <- data.frame(
    Stratifier=v,Group=d$group,Years=h/12,N=d$n,At_risk=d$at_risk,
    Observed_mortality=fmt(100*d$observed_mortality,1),
    Mean_out_of_fold_prediction=fmt(100*d$mean_predicted_mortality,1),
    Difference_percentage_points=fmt(100*(d$observed_mortality-d$mean_predicted_mortality),1))
}
if (!length(sublist)) stop("No complete eligible subgroup blocks; review suppression rules.")
save_table(do.call(rbind,sublist),"Table_S4_Subgroup_Calibration")

# Threshold-specific denominators remain private: nested differences can disclose small cells.
save_table(data.frame(Model=sensitivity$model,Recorded_size_range_mm=
    paste0("1-",sensitivity$max_mm),N=sensitivity$n,Deaths=sensitivity$deaths,
    Adjusted_HR_per_10mm=fmt(sensitivity$HR_per_10mm),
    Wald_95pct_CI=fmt_ci(sensitivity$lower95,sensitivity$upper95),
    Nonlinearity_LRT_p=fmt(sensitivity$spline_vs_linear_p,4)),
  "Table_S5_Tumor_Size_Sensitivity_REVIEW_BEFORE_RELEASE")

save_table(data.frame(Predictor=terms(reference_coefficients$term),
  Encoded_term=reference_coefficients$term,
  Beta=fmt(reference_coefficients$beta,9),
  Reference_design_value=fmt(reference_coefficients$reference_design_value,5)),
  "Table_S6a_Prediction_Equation_Terms")
save_table(data.frame(Horizon_months=reference_survival$horizon_months,
  Reference_survival=fmt(reference_survival$reference_survival,9),
  Equation="P(death by t) = 1 - S_reference(t)^exp(sum(beta*(x - x_reference)))"),
  "Table_S6b_Prediction_Equation_Baseline")

# Retrospective sample-size scenarios; these are not a post-hoc proof of adequacy.
if(!requireNamespace("pmsampsize",quietly=TRUE))
  stop("Install required reporting package: install.packages('pmsampsize')")
sample_rows <- list()
for(nm in c("A_expanded","D_common")) {
  d <- models$datasets[[nm]]
  fit <- models$fits[[nm]]
  p <- length(stats::coef(fit))
  rate <- sum(d$event)/sum(d$time_months/12)
  mean_followup <- mean(d$time_months)/12
  for(h in c(5,10)) for(r2 in c(.03,.05,.10,.15,.20)) {
    z <- tryCatch(pmsampsize::pmsampsize(type="s",csrsquared=r2,
      parameters=p,rate=rate,timepoint=h,meanfup=mean_followup,
      shrinkage=.9),error=function(e)e)
    if(inherits(z,"error")) stop("pmsampsize failed for ",nm,": ",conditionMessage(z))
    sample_rows[[length(sample_rows)+1L]] <- data.frame(Model=nm,
      Horizon_years=h,Parameters=p,Assumed_Cox_Snell_R2=r2,
      Required_N=as.numeric(z$sample_size)[1],Observed_N=nrow(d),
      Status="Hypothetical anticipated R2; retrospective scenario only")
  }
}
sample_out <- do.call(rbind,sample_rows)
# pmsampsize's governing shrinkage criterion yields identical requirements at
# both horizons under these assumptions; avoid 10 redundant manuscript rows.
for (nm in unique(sample_out$Model))
  for (r2 in unique(sample_out$Assumed_Cox_Snell_R2)) {
    x <- sample_out[sample_out$Model==nm &
           sample_out$Assumed_Cox_Snell_R2==r2,,drop=FALSE]
    if (length(unique(x$Required_N)) != 1L)
      stop("Sample-size horizon dependence found: retain separate rows.")
  }
sample_out <- sample_out[sample_out$Horizon_years == 5,,drop=FALSE]
sample_out$Horizon_years <- NULL
sample_out$Status <- paste0(sample_out$Status,
    "; same estimate for 5- and 10-year calculations under these inputs")
save_table(sample_out,"Table_S7_Sample_Size_Scenarios")

included <- is.finite(eligible$tumor_size_mm_raw) &
  eligible$tumor_size_mm_raw>=1 & eligible$tumor_size_mm_raw<=100 &
  stats::complete.cases(eligible[,c("AGE_NUM","SEX3","CDCC2"),drop=FALSE])
if(sum(included)!=nrow(A) || sum(eligible$event[included])!=sum(A$event))
  stop("Included/excluded flags do not reproduce primary cohort.")
groups <- list(Included=eligible[included,,drop=FALSE],
               Excluded=eligible[!included,,drop=FALSE])
s8 <- list()
for(g in names(groups)) {
  d <- groups[[g]]
  s8[[length(s8)+1L]] <- data.frame(Group=g,
    Characteristic=c("N","Deaths","Age, median [IQR]","Diagnosis year, median [IQR]",
                     "Size observed, N","Size, median [IQR] among observed"),
    Value=c(nrow(d),sum(d$event),fmt_quartiles(d$AGE_NUM),
      fmt_quartiles(d$YEAR_NUM),sum(is.finite(d$tumor_size_mm_raw)),
      fmt_quartiles(d$tumor_size_mm_raw)))
}
for(v in c("SEX3","CDCC2","RACE3","HISP2","INS3")) {
  lev <- levels(eligible[[v]])
  for(g in names(groups)) for(z in lev) {
    n <- sum(as.character(groups[[g]][[v]])==z,na.rm=TRUE)
    # Suppress the entire variable block across BOTH groups when ANY level is small.
    counts <- sapply(groups,function(d) table(factor(d[[v]],levels=lev)))
    safe <- all(counts>=11) && all(vapply(groups,function(d) {
      m <- sum(is.na(d[[v]]));m==0L || m>=11L
    },logical(1)))
    s8[[length(s8)+1L]] <- data.frame(Group=g,
      Characteristic=paste0(v,": ",z),
      Value=if(safe) as.character(n) else "Withheld (whole variable block)")
  }
}
save_table(do.call(rbind,s8),"Table_S8_Included_vs_Excluded")
save_table(data.frame(Domain=c("Histology","Metastatic status"),
  Reporting_status=c(
    "Eligible histology codes fixed at 9370/3, 9371/3, 9372/3; small-cell review required before releasing subtype counts",
    "Not reconstructed in public cohort script; do not report subtype counts or adjusted effect until field-specific audit")),
  "QC_Histology_Metastasis_Awaiting_Stage08")

# Combine both 5- and 10-year smoothed calibration curves produced by stage 06.
curves <- do.call(rbind,lapply(HORIZONS,function(h) {
  d <- load_csv(paste0("calibration_curve_",h,"m.csv"))
  require_fields(d,c("predicted_risk","observed_smoothed"),"Calibration curve")
  data.frame(Years=h/12,predicted=d$predicted_risk,observed=d$observed_smoothed)
}))
p1 <- ggplot2::ggplot(curves,ggplot2::aes(x=predicted,y=observed))+
  ggplot2::geom_abline(slope=1,intercept=0,linetype="dashed",color="grey55")+
  ggplot2::geom_line(linewidth=.9,color="black")+
  ggplot2::facet_wrap(~Years,nrow=1,
      labeller=ggplot2::labeller(Years=function(x) paste0(x," years")))+
  ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  ggplot2::labs(x="Predicted all-cause mortality",y="IPCW-smoothed observed mortality",
    caption="Internal out-of-fold cross-fitting; smooth curves without pointwise confidence bands")+
  ggplot2::theme_bw(base_size=12)
sensitivity$Model_label <- factor(sensitivity$model,
    levels=c("A_clinical","D_full"),
    labels=c("Clinical model","Clinical + sociodemographic model"))
sensitivity$Size_limit <- factor(sensitivity$max_mm,
    levels=c(988,150,100), labels=c("988","150","100"))
p2 <- ggplot2::ggplot(sensitivity,
    ggplot2::aes(x=HR_per_10mm,y=Size_limit))+
  ggplot2::geom_vline(xintercept=1,linetype="dashed",color="grey55")+
  ggplot2::geom_errorbar(ggplot2::aes(xmin=lower95,xmax=upper95),
    orientation="y",width=.15)+
  ggplot2::geom_point(size=2)+
  ggplot2::facet_wrap(~Model_label,nrow=1)+
  ggplot2::labs(x="Adjusted tumor-size HR per 10 mm (95% Wald CI)",
    y="Maximum included size (mm)",
    caption="Different size thresholds change the eligible cohort; comparisons are descriptive")+
  ggplot2::theme_bw(base_size=12)
save_fig <- function(p,name,w=10,h=6) {
  ggplot2::ggsave(file.path(figdir,paste0(name,".png")),p,
                  width=w,height=h,dpi=300,bg="white")
  ggplot2::ggsave(file.path(figdir,paste0(name,".pdf")),p,
                  width=w,height=h,bg="white")
}
save_fig(p1,"Figure_1_Smooth_Calibration",10.5,5)
save_fig(p2,"Figure_2_Tumor_Size_Sensitivity",10,5)

# Risk-tertile Kaplan-Meier curves illustrate development-data separation only.
p120 <- as.numeric(riskRegression::predictRisk(models$fits$A_expanded,
                              newdata=as.data.frame(A),times=120))
if(length(p120)!=nrow(A) || any(!is.finite(p120))) stop("Invalid fitted risks.")
q <- as.numeric(stats::quantile(p120,c(1/3,2/3)))
if(q[1]>=q[2]) stop("Non-distinct risk-tertile cut points.")
A$risk_group <- factor(cut(p120,c(-Inf,q,Inf),
  labels=c("Low","Intermediate","High")),
  levels=c("Low","Intermediate","High"))
km <- survival::survfit(survival::Surv(time_months,event)~risk_group,data=A)
z <- summary(km)
kmd <- data.frame(month=z$time,overall_survival=z$surv,
  group=sub("^risk_group=","",as.character(z$strata)))
kmd <- rbind(data.frame(month=0,overall_survival=1,
  group=levels(A$risk_group)),kmd)
rt <- do.call(rbind,lapply(levels(A$risk_group),function(g) {
  sub <- A[A$risk_group==g,,drop=FALSE]
  data.frame(group=g,month=c(0,60,120),at_risk=vapply(c(0,60,120),
    function(t) sum(sub$time_months>=t),integer(1)))
}))
allowed <- c(0,60,120)[vapply(c(0,60,120),function(t)
  all(rt$at_risk[rt$month==t]>=11),logical(1))]
rt <- rt[rt$month %in% allowed,,drop=FALSE]
save_table(rt,"QC_Figure_S1_Numbers_at_Risk")
pS1 <- ggplot2::ggplot(kmd,ggplot2::aes(x=month,y=overall_survival,
            color=group,group=group))+
  ggplot2::geom_step(linewidth=.9)+
  ggplot2::coord_cartesian(xlim=c(0,120),ylim=c(0,1))+
  ggplot2::labs(x="Months from diagnosis",y="Overall survival",
    color="Fitted risk tertile",
    caption="Development-data KM curves; illustrative only, not external validation")+
  ggplot2::theme_bw(base_size=12)
if (!requireNamespace("patchwork",quietly=TRUE))
  stop("Install required figure package: install.packages('patchwork')")
risk_plot <- ggplot2::ggplot(rt,ggplot2::aes(x=month,y=group,label=at_risk))+
  ggplot2::geom_text(size=3.6)+
  ggplot2::scale_x_continuous(breaks=allowed,limits=c(0,120))+
  ggplot2::labs(x="Months from diagnosis",y="Number at risk")+
  ggplot2::theme_minimal(base_size=11)+
  ggplot2::theme(panel.grid=ggplot2::element_blank())
save_fig(patchwork::wrap_plots(pS1,risk_plot,ncol=1,heights=c(3.8,1.2)),
         "Figure_S1_Illustrative_Risk_Tertiles",10,7.5)

# Connected flow with non-overlapping exclusion totals (not separate cause counts).
if(nrow(flow)!=4L || any(diff(flow$n)>0))
  stop("Unexpected cohort-flow structure.")
flow$label <- paste0(flow$stage,"\nN = ",format(flow$n,big.mark=","))
flow$position <- -seq_len(nrow(flow))
attrition <- data.frame(position=c(-1.5,-2.5,-3.5),
  label=paste0("Not retained: n = ",
    format(-diff(flow$n),big.mark=",")))
pS2 <- ggplot2::ggplot()+
  ggplot2::geom_segment(data=attrition,
    ggplot2::aes(x=0,xend=0,y=position+.28,yend=position-.28),
    arrow=grid::arrow(length=grid::unit(.15,"inches")),linewidth=.6)+
  ggplot2::geom_label(data=flow,
    ggplot2::aes(x=0,y=position,label=label),size=3.15,
    fill="white",label.padding=grid::unit(.4,"lines"))+
  ggplot2::geom_text(data=attrition,
    ggplot2::aes(x=1.2,y=position,label=label),
    hjust=0,size=2.9)+
  ggplot2::coord_cartesian(xlim=c(-2.8,3.8),ylim=c(-4.4,-.6),clip="off")+
  ggplot2::theme_void(base_size=12)+
  ggplot2::labs(title="NCDB cohort selection",
    subtitle="Attrition reflects combined eligibility and completeness criteria")
save_fig(pS2,"Figure_S2_Cohort_Flow",11,6.5)

# Optional DOCX generation from CSVs; never hand-enter the manuscript estimates.
if (requireNamespace("officer",quietly=TRUE) &&
    requireNamespace("flextable",quietly=TRUE)) {
  for (f in list.files(tabdir,pattern="^Table_.*\\.csv$",full.names=TRUE)) {
    d <- utils::read.csv(f,check.names=FALSE,stringsAsFactors=FALSE)
    ft <- flextable::autofit(flextable::theme_booktabs(flextable::flextable(d)))
    flextable::save_as_docx(ft,path=file.path(tabdir,
      paste0(tools::file_path_sans_ext(basename(f)),".docx")))
  }
}
writeLines(c("Canonical performance source: this 2024 PUF GitHub pipeline rerun.",
 "Generated files remain PRIVATE pending NCDB disclosure review.",
 "Tables 1-4 and S1-S9; Figures 1-2 and S1-S2 generated.",
 "Figure S1 number-at-risk data are in QC_Figure_S1_Numbers_at_Risk.csv.",
 "Figure S1 KM is illustrative development-data stratification, not validation.",
 "Rare histology/metastasis results require additional source-field adjudication.",
 "Internal bootstrap/OOF evaluations are not external validation.",
 paste("R:",R.version.string),
 paste("Generated (UTC):",format(Sys.time(),tz="UTC",usetz=TRUE))),
 file.path(target,"REPORTING_MANIFEST.txt"))
cat("Reporting exhibits generated in ",target,"\n",
    "Disclosure review required before distributing output files.\n",sep="")
