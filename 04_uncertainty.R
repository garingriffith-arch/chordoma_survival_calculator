# 04. Nested bootstrap uncertainty for optimism-corrected performance.
# Outer resamples estimate variability; inner resamples estimate optimism.
# The paired full-versus-clinical contrast uses the SAME common-cohort indices.
source(file.path("R","00_config.R"))
mod<-readRDS(file.path(PRIVATE_DIR,"models_INTERNAL_ONLY.rds"))
dlist<-mod$datasets;forms<-mod$formulas
outer_B<-if(QUICK)8L else 120L
inner_B<-if(QUICK)4L else 30L
cores<-parallel::detectCores(logical=TRUE)
target<-max(1L,floor(0.70*cores))
workers<-min(target,outer_B)
clim<-suppressWarnings(as.integer(Sys.getenv("CHORDOMA_MAX_WORKERS","0")))
if(is.finite(clim)&&clim>0) workers<-min(workers,clim)
meta<-list(B_outer=outer_B,B_inner=inner_B,
           cohort_hash=unname(tools::md5sum(file.path(PRIVATE_DIR,"cohort_INTERNAL_ONLY.rds"))),
           seed=SEED+1550L)
cpfile<-file.path(PRIVATE_DIR,paste0("nested_",outer_B,"x",inner_B,"_INTERNAL_ONLY.rds"))
if(file.exists(cpfile)) {
  old<-readRDS(cpfile)
  if(!identical(old$meta,meta)) stop("Nested checkpoint metadata mismatch")
  runs<-old$runs
} else runs<-vector("list",outer_B)
metric_one<-function(f,d){
  lp<-as.numeric(stats::predict(f,newdata=d,type="lp"))
  cc<-as.numeric(survival::concordance(survival::Surv(d$time_months,d$event)~I(-lp))$concordance)
  z<-riskRegression::Score(object=list(M=f),formula=Hist(time_months,event)~1,
       data=as.data.frame(d),times=c(60,120),
       metrics=c("auc","brier"),cens.model="km",null.model=FALSE,conf.int=FALSE)
  grab<-function(part,measure) {
    x<-as.data.frame(z[[part]]$score)
    x<-x[as.character(x$model)=="M",,drop=FALSE]
    tm<-intersect(c("time","times"),names(x))[1]
    as.numeric(x[[measure]][match(c(60,120),as.numeric(x[[tm]]))])
  }
  a<-grab("AUC","AUC"); b<-grab("Brier","Brier")
  result<-c(Harrell_C=cc,AUC_60=a[1],AUC_120=a[2],Brier_60=b[1],Brier_120=b[2])
  if(any(!is.finite(result))) stop("Non-finite performance metric")
  result
}
fit_one<-function(f,d) survival::coxph(f,data=d,ties="efron",
                                      model=TRUE,x=TRUE,y=TRUE,singular.ok=FALSE)
eval_outer<-function(b){
  set.seed(meta$seed+b)
  ia<-sample.int(nrow(dlist$A_expanded),nrow(dlist$A_expanded),replace=TRUE)
  ic<-sample.int(nrow(dlist$D_common),nrow(dlist$D_common),replace=TRUE)
  od<-list(A_expanded=dlist$A_expanded[ia,,drop=FALSE],
           A_common=dlist$D_common[ic,,drop=FALSE],
           D_common=dlist$D_common[ic,,drop=FALSE])
  nm<-names(od);app<-list()
  for(k in nm){
    f<-fit_one(forms[[k]],od[[k]])
    app[[k]]<-metric_one(f,od[[k]])
  }
  optimism<-lapply(nm,function(x)matrix(NA_real_,inner_B,5,
       dimnames=list(NULL,names(app[[x]]))))
  names(optimism)<-nm
  ok<-logical(inner_B)
  for(i in seq_len(inner_B)){
    ixA<-sample.int(nrow(od$A_expanded),nrow(od$A_expanded),replace=TRUE)
    ixC<-sample.int(nrow(od$D_common),nrow(od$D_common),replace=TRUE)
    inds<-list(A_expanded=ixA,A_common=ixC,D_common=ixC)
    z<-tryCatch(lapply(nm,function(k){
      train<-od[[k]][inds[[k]],,drop=FALSE]
      f<-fit_one(forms[[k]],train)
      metric_one(f,train)-metric_one(f,od[[k]])
    }),error=function(e)NULL)
    if(is.null(z))next
    names(z)<-nm
    for(k in nm)optimism[[k]][i,]<-z[[k]]
    ok[i]<-TRUE
  }
  if(sum(ok)<ceiling(.70*inner_B)) stop("Insufficient successful inner resamples")
  do.call(rbind,lapply(nm,function(k){
    cor<-app[[k]]-colMeans(optimism[[k]][ok,,drop=FALSE])
    data.frame(outer_iteration=b,model=k,metric=names(cor),
               internally_corrected=as.numeric(cor),
               inner_success=sum(ok))
  }))
}
set.seed(meta$seed)
pending<-which(vapply(runs,is.null,logical(1)))
if(length(pending)){
  if(workers==1L){
    for(b in pending){
      runs[[b]]<-tryCatch(list(ok=TRUE,data=eval_outer(b)),
             error=function(e)list(ok=FALSE,message=conditionMessage(e)))
      saveRDS(list(meta=meta,runs=runs),cpfile)
      cat("Nested outer",b,"/",outer_B,"\n")
    }
  } else {
    cl<-parallel::makeCluster(workers,type="PSOCK")
    parallel::clusterEvalQ(cl,{
      library(survival);library(prodlim);library(riskRegression);NULL
    })
    parallel::clusterExport(cl,c("dlist","forms","meta","outer_B","inner_B",
       "metric_one","fit_one","eval_outer"),envir=environment())
    for(k in split(pending,ceiling(seq_along(pending)/workers))){
      res<-parallel::parLapply(cl,k,function(b)
        tryCatch(list(ok=TRUE,data=eval_outer(b)),
          error=function(e)list(ok=FALSE,message=conditionMessage(e))))
      for(j in seq_along(k))runs[[k[j]]]<-res[[j]]
      saveRDS(list(meta=meta,runs=runs),cpfile)
      cat("Nested outer completed:",sum(!vapply(runs,is.null,logical(1))),"/",outer_B,"\n")
    }
    parallel::stopCluster(cl)
  }
}
ok<-which(vapply(runs,function(z)isTRUE(z$ok),logical(1)))
if(length(ok)<ceiling(.85*outer_B))stop("Too many failed outer resamples; no CI exported")
long<-do.call(rbind,lapply(ok,function(i)runs[[i]]$data))
ci<-do.call(rbind,lapply(split(long,list(long$model,long$metric),drop=TRUE),function(d){
  q<-stats::quantile(d$internally_corrected,c(.025,.975),na.rm=TRUE)
  data.frame(model=d$model[1],metric=d$metric[1],lower95=q[1],upper95=q[2])
}))
write_aggregate(ci,"nested_corrected_performance_CI.csv")
pair<-list()
for(metric in unique(long$metric)){
  da<-long[long$model=="A_common"&long$metric==metric,,drop=FALSE]
  db<-long[long$model=="D_common"&long$metric==metric,,drop=FALSE]
  x<-merge(da,db,by="outer_iteration",suffixes=c("_A","_D"))
  delta<-x$internally_corrected_D-x$internally_corrected_A
  q<-stats::quantile(delta,c(.025,.975),na.rm=TRUE)
  pair[[metric]]<-data.frame(metric=metric,
       point_estimate_outer_mean=mean(delta),
       lower95=q[1],upper95=q[2],successful_pairs=nrow(x))
}
write_aggregate(do.call(rbind,pair),"nested_paired_D_minus_A_CI.csv")
write_aggregate(data.frame(outer_requested=outer_B,outer_successful=length(ok),
                           inner_requested=inner_B,
                           interval_type="Approximate nested-bootstrap percentile interval"),
                "nested_bootstrap_QC.csv")
cat("Nested uncertainty completed; confidence intervals are internal only.\n")
