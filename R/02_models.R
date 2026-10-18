# 02. Cox development, reference survival, and equation verification.
source(file.path("R","00_config.R"))
o<-readRDS(file.path(PRIVATE_DIR,"cohort_INTERNAL_ONLY.rds"))
A<-o$A;D<-o$D
fa<-survival::Surv(time_months,event)~AGE_NUM+tumor_size_10mm+SEX3+CDCC2
fd<-survival::Surv(time_months,event)~AGE_NUM+tumor_size_10mm+SEX3+CDCC2+
  RACE3+HISP2+INS3+MED_INC_QUAR
formulas<-list(A_expanded=fa,A_common=fa,D_common=fd)
datasets<-list(A_expanded=A,A_common=D,D_common=D)
fit_cox<-function(f,d) survival::coxph(f,data=d,ties="efron",x=TRUE,y=TRUE,
                                       model=TRUE,singular.ok=FALSE)
fits<-Map(fit_cox,formulas,datasets)
names(fits)<-names(formulas)
terms<-do.call(rbind,lapply(names(fits),function(nm){
  fit<-fits[[nm]]; beta<-coef(fit);se<-sqrt(diag(vcov(fit)))
  data.frame(model=nm,term=names(beta),beta=as.numeric(beta),SE=as.numeric(se),
             HR=exp(as.numeric(beta)),low95=exp(as.numeric(beta)-1.96*se),
             high95=exp(as.numeric(beta)+1.96*se))
}))
write_aggregate(terms,"cox_coefficients.csv")
reference<-data.frame(AGE_NUM=52,tumor_size_10mm=3.2,
                      SEX3=factor("Male",levels=levels(A$SEX3)),
                      CDCC2=factor("0",levels=levels(A$CDCC2)))
sf<-survival::survfit(fits$A_expanded,newdata=reference)
s<-summary(sf,times=HORIZONS,extend=TRUE)
baseref<-data.frame(horizon_months=HORIZONS,
                    reference_survival=as.numeric(s$surv))
write_aggregate(baseref,"baseline_survival_reference.csv")
fixed<-read.csv(file.path(ROOT,"model","coefficients.csv"),check.names=FALSE)
b<-coef(fits$A_expanded)
if(!identical(names(b),fixed$term)) stop("Primary Cox model design does not match released equation terms")
if(max(abs(unname(b)-fixed$beta))>1e-4) {
  stop("Locked coefficient parity FAILED; do not update calculator or manuscript")
}
fixed_s<-read.csv(file.path(ROOT,"model","reference_survival.csv"))
if(!identical(as.integer(baseref$horizon_months),as.integer(fixed_s$horizon_months))||
   max(abs(baseref$reference_survival-fixed_s$reference_survival))>1e-4) {
  stop("Reference survival parity FAILED; do not update calculator")
}
# The deployment equation is evaluated without a private patient-level RDS.
reference_lp<-as.numeric(stats::predict(fits$A_expanded,newdata=reference,type="lp"))
rows<-A[seq_len(min(50L,nrow(A))),,drop=FALSE]
cox_lp<-as.numeric(stats::predict(fits$A_expanded,newdata=rows,type="lp"))
X<-stats::model.matrix(stats::delete.response(stats::terms(fits$A_expanded)),rows)
X<-X[,names(b),drop=FALSE]
design<-c(52,3.2,0,0,0);names(design)<-names(b)
eq_lp<-as.numeric(sweep(X,2,design,"-")%*%b)
# Coxph may internally center the linear predictor; centered contrasts must agree.
delta<-max(abs((cox_lp-reference_lp)-eq_lp))
if(!is.finite(delta)||delta>1e-6) stop("Equation versus fitted-model linear predictor mismatch: ",delta)
write_aggregate(data.frame(test="reference-centered_linear_predictor",
                           maximum_absolute_error=delta),"equation_check.csv")
write_private(list(fits=fits,formulas=formulas,datasets=datasets,
                   reference=reference,baseline=baseref),"models_INTERNAL_ONLY.rds")
cat("Primary coefficients and baseline survival agree with locked public equation.\n")
