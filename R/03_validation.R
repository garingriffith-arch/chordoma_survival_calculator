# 03. Bootstrap optimism correction and internal held-out calibration.
# All resampling is internal to the NCDB development data.
source(file.path("R","00_config.R"))
oo<-readRDS(file.path(PRIVATE_DIR,"models_INTERNAL_ONLY.rds"))
fits<-oo$fits;ff<-oo$formulas;dd<-oo$datasets
fit_one<-function(f,d) survival::coxph(f,data=d,ties="efron",x=TRUE,y=TRUE,
                                       model=TRUE,singular.ok=FALSE)
cstat<-function(f,d) {
  lp<-as.numeric(stats::predict(f,newdata=d,type="lp"))
  as.numeric(survival::concordance(survival::Surv(d$time_months,d$event)~I(-lp))$concordance)
}
metrics<-function(f,d) {
  x<-riskRegression::Score(object=list(M=f),
    formula=Hist(time_months,event)~1,data=as.data.frame(d),times=HORIZONS,
    metrics=c("auc","brier"),null.model=FALSE,cens.model="km",
    conf.int=FALSE)
  score_values<-function(whichmetric,column){
    z<-as.data.frame(x[[whichmetric]]$score)
    z<-z[as.character(z$model)=="M",,drop=FALSE]
    tcol<-intersect(c("times","time"),names(z))[1]
    if(length(tcol)!=1L||is.na(tcol)) stop("Unexpected riskRegression::Score structure")
    out<-z[[column]][match(HORIZONS,as.numeric(z[[tcol]]))]
    if(length(out)!=2L||any(!is.finite(out))) stop("Incomplete Score metrics")
    as.numeric(out)
  }
  c(Harrell_C=cstat(f,d),
    AUC_60=score_values("AUC","AUC")[1],
    AUC_120=score_values("AUC","AUC")[2],
    Brier_60=score_values("Brier","Brier")[1],
    Brier_120=score_values("Brier","Brier")[2])
}
app<-lapply(names(fits),function(nm)metrics(fits[[nm]],dd[[nm]]))
names(app)<-names(fits)
write_aggregate(do.call(rbind,lapply(names(app),function(nm)
  data.frame(model=nm,metric=names(app[[nm]]),estimate=as.numeric(app[[nm]])))),
  "performance_apparent.csv")
checkpoint<-file.path(PRIVATE_DIR,paste0("bootstrap_B",B_BOOT,"_INTERNAL_ONLY.rds"))
cohort_hash<-unname(tools::md5sum(file.path(PRIVATE_DIR,"cohort_INTERNAL_ONLY.rds")))
meta<-list(B=B_BOOT,hash=cohort_hash,seed=SEED+310L)
if(file.exists(checkpoint)){
  ck<-readRDS(checkpoint)
  if(!identical(ck$meta,meta)) stop("Existing bootstrap checkpoint is for another data version")
  reps<-ck$reps
} else reps<-vector("list",B_BOOT)
set.seed(meta$seed)
seeds<-sample.int(.Machine$integer.max,B_BOOT)
for(b in seq_len(B_BOOT)){
  if(!is.null(reps[[b]])) next
  set.seed(seeds[b])
  ia<-sample.int(nrow(dd$A_expanded),nrow(dd$A_expanded),replace=TRUE)
  ic<-sample.int(nrow(dd$A_common),nrow(dd$A_common),replace=TRUE)
  idx<-list(A_expanded=ia,A_common=ic,D_common=ic)
  reps[[b]]<-tryCatch({
    out<-lapply(names(fits),function(nm){
      tr<-dd[[nm]][idx[[nm]],,drop=FALSE]
      m<-fit_one(ff[[nm]],tr)
      metrics(m,tr)-metrics(m,dd[[nm]])
    })
    names(out)<-names(fits)
    out
  },error=function(e)list(error=conditionMessage(e)))
  if(b%%10L==0L || b==B_BOOT){
    saveRDS(list(meta=meta,reps=reps),checkpoint)
    cat("Bootstrap",b,"/",B_BOOT,"\n")
  }
}
valid<-Filter(function(x)is.null(x$error),reps)
if(length(valid)<ceiling(0.80*B_BOOT)) stop("Too many bootstrap failures; inspect private checkpoint")
optimism<-lapply(names(fits),function(nm){
  z<-do.call(rbind,lapply(valid,function(x)x[[nm]]))
  colMeans(z)
})
names(optimism)<-names(fits)
corrected<-do.call(rbind,lapply(names(fits),function(nm)
  data.frame(model=nm,metric=names(app[[nm]]),
             apparent=as.numeric(app[[nm]]),
             mean_optimism=as.numeric(optimism[[nm]]),
             optimism_corrected=as.numeric(app[[nm]]-optimism[[nm]]))))
write_aggregate(corrected,"performance_optimism_corrected.csv")
write_aggregate(data.frame(requested=B_BOOT,successful=length(valid),
                           failed=B_BOOT-length(valid),
                           method="bootstrap apparent minus mean train-test optimism"),
                "bootstrap_qc.csv")

# Repeated K-fold cross-fitting yields predictions for patients excluded
# from each fold's model fitting; this is internal, not external validation.
oofcheckpoint<-file.path(PRIVATE_DIR,paste0("oof_R",CV_REPEATS,"_INTERNAL_ONLY.rds"))
cvmeta<-list(hash=cohort_hash,repeats=CV_REPEATS,folds=CV_FOLDS,
             horizon=HORIZONS,seed=SEED+490L)
if(file.exists(oofcheckpoint)){
  old<-readRDS(oofcheckpoint)
  if(!identical(old$meta,cvmeta)) stop("OOF checkpoint has different cohort/parameters")
  oof<-old$predictions
} else oof<-list()
set.seed(cvmeta$seed)
for(nm in c("A_expanded","A_common","D_common")){
  d<-dd[[nm]];n<-nrow(d)
  mat<-matrix(0,n,2)
  for(repno in seq_len(CV_REPEATS)){
    folds<-sample(rep(seq_len(CV_FOLDS),length.out=n))
    one<-matrix(NA_real_,n,2)
    for(k in seq_len(CV_FOLDS)){
      train<-d[folds!=k,,drop=FALSE]
      test<-d[folds==k,,drop=FALSE]
      f<-fit_one(ff[[nm]],train)
      p<-riskRegression::predictRisk(f,newdata=as.data.frame(test),times=HORIZONS)
      p<-as.matrix(p)
      if(!identical(dim(p),c(nrow(test),2L))) stop("Unexpected OOF risk matrix")
      one[folds==k,]<-p
    }
    if(any(!is.finite(one))) stop("Incomplete OOF predictions")
    mat<-mat+one/CV_REPEATS
  }
  oof[[nm]]<-mat
  saveRDS(list(meta=cvmeta,predictions=oof),oofcheckpoint)
}
# Grouped observed survival is KM, estimated only when follow-up reaches the horizon.
rows<-list()
for(nm in names(oof))for(j in seq_along(HORIZONS)){
  d<-dd[[nm]];h<-HORIZONS[j]
  for(groupvar in c("HISP2","INS3","RACE3","SEX3")){
    for(level in levels(d[[groupvar]])){
      ix<-which(!is.na(d[[groupvar]]) & as.character(d[[groupvar]])==level)
      if(length(ix)<11L) next
      sub<-d[ix,,drop=FALSE]
      if(sum(sub$event)<11L || sum(sub$time_months>=h)<11L) next
      km<-survival::survfit(survival::Surv(time_months,event)~1,data=sub)
      s<-summary(km,times=h,extend=FALSE)
      if(!length(s$surv)) next
      rows[[length(rows)+1L]]<-data.frame(model=nm,horizon_months=h,
        stratifier=groupvar,group=level,n=length(ix),
        at_risk=sum(sub$time_months>=h),
        observed_mortality=1-as.numeric(s$surv),
        mean_predicted_mortality=mean(oof[[nm]][ix,j]))
    }
  }
}
if(length(rows)) write_aggregate(do.call(rbind,rows),"subgroup_calibration_INTERNAL_DISCLOSURE_REVIEW.csv")
write_aggregate(data.frame(repeats=CV_REPEATS,folds=CV_FOLDS,
                           A_expanded_N=nrow(dd$A_expanded),
                           D_common_N=nrow(dd$D_common)),
                "crossfit_qc.csv")
cat("Bootstrap validation and cross-fitted predictions finished.\n")
