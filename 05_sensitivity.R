# 04. Reviewer-requested tumor-size and missing-data sensitivity.
# Neither analysis replaces the 100-mm complete-case primary prediction model.
source(file.path("R","00_config.R"))
o<-readRDS(file.path(PRIVATE_DIR,"cohort_INTERNAL_ONLY.rds"))
eligible<-o$eligible
clinical<-c("AGE_NUM","tumor_size_10mm","SEX3","CDCC2")
social<-c("RACE3","HISP2","INS3","MED_INC_QUAR")
fa<-survival::Surv(time_months,event)~AGE_NUM+tumor_size_10mm+SEX3+CDCC2
fd<-survival::Surv(time_months,event)~AGE_NUM+tumor_size_10mm+SEX3+CDCC2+
  RACE3+HISP2+INS3+MED_INC_QUAR
rows<-list()
for(maxmm in c(100L,150L,988L)){
  x<-eligible[is.finite(eligible$tumor_size_mm_raw)&
              eligible$tumor_size_mm_raw>=1&eligible$tumor_size_mm_raw<=maxmm,,drop=FALSE]
  for(nm in c("A_clinical","D_full")){
    need<-if(nm=="A_clinical")clinical else c(clinical,social)
    d<-x[complete.cases(x[,need,drop=FALSE]),,drop=FALSE]
    f<-survival::coxph(if(nm=="A_clinical")fa else fd,
                       data=d,ties="efron",model=TRUE,x=TRUE,y=TRUE)
    beta<-coef(f)["tumor_size_10mm"]
    se<-sqrt(vcov(f)["tumor_size_10mm","tumor_size_10mm"])
    ns<-survival::coxph(update(if(nm=="A_clinical")fa else fd,
        .~.-tumor_size_10mm+splines::ns(tumor_size_10mm,df=3)),
        data=d,ties="efron")
    linear_lr<-as.numeric(2*(ns$loglik[2]-f$loglik[2]))
    p_nonlinear<-stats::pchisq(max(0,linear_lr),df=2,lower.tail=FALSE)
    rows[[length(rows)+1L]]<-data.frame(model=nm,max_mm=maxmm,
      n=nrow(d),deaths=sum(d$event),HR_per_10mm=unname(exp(beta)),
      lower95=unname(exp(beta-1.96*se)),upper95=unname(exp(beta+1.96*se)),
      spline_vs_linear_p=p_nonlinear)
  }
}
write_aggregate(do.call(rbind,rows),"tumor_size_sensitivity.csv")

# MI is coefficient stability analysis only. It is not a validated deployable model.
# Exclude known size >100 mm, but retain missing-size records for imputation.
mi<-eligible[is.na(eligible$tumor_size_mm_raw)|
             (eligible$tumor_size_mm_raw>=1&eligible$tumor_size_mm_raw<=100),,drop=FALSE]
features<-c(clinical,social)
sf<-survival::survfit(survival::Surv(time_months,event)~1,data=mi)
haz<-cumsum(sf$n.event/pmax(sf$n.risk,1))
pos<-findInterval(mi$time_months,sf$time)
mi$na_hazard<-ifelse(pos==0,0,haz[pmax(pos,1)])
X<-as.data.frame(mi[,c(features,"YEAR_NUM","event","na_hazard","time_months"),drop=FALSE])
if(anyNA(X$event)|anyNA(X$time_months)) stop("MI outcomes cannot be missing")
methods<-mice::make.method(X);methods[]<-""
for(v in features){
  if(!anyNA(X[[v]])) next
  if(is.factor(X[[v]])){
    methods[v]<-if(nlevels(X[[v]])==2L)"logreg" else "polyreg"
  } else methods[v]<-"pmm"
}
pm<-mice::make.predictorMatrix(X);pm[,]<-0
pm[features,c(features,"YEAR_NUM","event","na_hazard")]<-1
diag(pm)<-0
pm[c("YEAR_NUM","event","na_hazard","time_months"),]<-0
if(anyNA(X$YEAR_NUM)) pm[,"YEAR_NUM"]<-0
checkpoint<-file.path(PRIVATE_DIR,paste0("imputation_m",MI_M,"_INTERNAL_ONLY.rds"))
mi_meta<-list(m=MI_M,n=nrow(X),seed=SEED+1193L,
              cohort_hash=unname(tools::md5sum(file.path(PRIVATE_DIR,"cohort_INTERNAL_ONLY.rds"))))
if(file.exists(checkpoint)){
  z<-readRDS(checkpoint)
  if(!identical(z$meta,mi_meta)) stop("Existing imputation checkpoint does not match")
  imp<-z$imp
} else {
  imp<-mice::mice(X,m=MI_M,maxit=if(QUICK)3L else 12L,
                  method=methods,predictorMatrix=pm,seed=mi_meta$seed,
                  printFlag=FALSE)
  saveRDS(list(meta=mi_meta,imp=imp),checkpoint)
}
pooled<-list()
for(nm in c("A_clinical","D_full")){
  f<-if(nm=="A_clinical")fa else fd
  ff<-lapply(seq_len(imp$m),function(i){
    d<-mice::complete(imp,i)
    if(anyNA(d[,features,drop=FALSE])) stop("Missing imputed predictors")
    if(any(d$tumor_size_10mm<0.1|d$tumor_size_10mm>10)) stop("MI size outside 1-100mm")
    survival::coxph(f,data=d,ties="efron")
  })
  summ<-as.data.frame(summary(mice::pool(mice::as.mira(ff)),conf.int=TRUE))
  summ$model<-nm
  pooled[[nm]]<-summ
}
write_aggregate(do.call(rbind,pooled),"MI_pooled_coefficients.csv")
write_aggregate(data.frame(M=MI_M,N=nrow(X),
  scope="Coefficient sensitivity only; imputations use survival outcomes, not a patient-level deployment method"),
  "MI_scope.csv")
cat("Tumor-size and multiple-imputation sensitivity analyses finished.\n")
