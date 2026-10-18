# 05. Figures and final aggregate reporting are generated from upstream results.
source(file.path("R","00_config.R"))
models<-readRDS(file.path(PRIVATE_DIR,"models_INTERNAL_ONLY.rds"))
o<-readRDS(file.path(PRIVATE_DIR,"cohort_INTERNAL_ONLY.rds"))
coefs<-read.csv(file.path(AGG_DIR,"cox_coefficients.csv"))
perf<-read.csv(file.path(AGG_DIR,"performance_optimism_corrected.csv"))
sens<-read.csv(file.path(AGG_DIR,"tumor_size_sensitivity.csv"))
# Model-design and survival reference tables are exported in 02_models.R.
primary<-coefs[coefs$model=="A_expanded",,drop=FALSE]
common<-coefs[coefs$model=="D_common",,drop=FALSE]
write_aggregate(primary,"Table_2_primary_model.csv")
write_aggregate(common,"Table_S2_full_comparator.csv")
write_aggregate(perf,"Table_3_validation_and_comparison.csv")
write_aggregate(sens,"Table_S5_size_sensitivity.csv")
# Reverse KM follow-up and horizon risk sets.
followup<-do.call(rbind,lapply(c("A_expanded","D_common"),function(nm){
  d<-models$datasets[[nm]]
  fk<-survival::survfit(survival::Surv(time_months,1L-event)~1,data=d)
  # survfit() does not supply a $table field. The median is in summary(fit)$table.
  tab <- summary(fk)$table
  if (is.matrix(tab)) {
    med <- if ("median" %in% colnames(tab)) as.numeric(tab[1L, "median"]) else NA_real_
  } else {
    med <- if ("median" %in% names(tab)) as.numeric(tab[["median"]]) else NA_real_
  }
  # NA is valid when reverse-KM follow-up never reaches 50%; numeric(0) is not.
  if (length(med) != 1L || !is.finite(med)) med <- NA_real_
  data.frame(model=nm,n=nrow(d),deaths=sum(d$event),
    median_followup_reverse_KM_months=med,
    at_risk_60=sum(d$time_months>=60),
    at_risk_120=sum(d$time_months>=120))
}))
write_aggregate(followup,"Table_1_followup.csv")
# Out-of-fold mean risk vs observed KM, with confidence bounds conditional on
# cross-fitted predictions rather than on the model-development process.
p<-file.path(PRIVATE_DIR,paste0("oof_R",CV_REPEATS,"_INTERNAL_ONLY.rds"))
ck<-readRDS(p)
d<-models$datasets$A_expanded;pred<-ck$predictions$A_expanded
cal<-list()
for(j in seq_along(HORIZONS)){
  h<-HORIZONS[j]
  km<-survival::survfit(survival::Surv(time_months,event)~1,data=d)
  observed<-1-as.numeric(summary(km,times=h,extend=FALSE)$surv)
  cal[[j]]<-data.frame(horizon_months=h,
      observed_KM_mortality=observed,
      mean_OOF_predicted_mortality=mean(pred[,j]),
      O_over_E=observed/mean(pred[,j]))
}
write_aggregate(do.call(rbind,cal),"Table_calibration_in_the_large.csv")

# Plot calibration using IPCW-weighted binary status at a fixed horizon.
# Confidence intervals here are not full model-refit bootstrap intervals.
for(j in seq_along(HORIZONS)){
  h<-HORIZONS[j];risk<-pred[,j]
  sf<-survival::survfit(survival::Surv(time_months,1L-event)~1,data=d)
  ctime<-sf$time
  survc<-sf$surv
  G<-function(x){
    i<-findInterval(x,ctime)
    ifelse(i==0,1,survc[pmax(i,1)])
  }
  early_event<-d$event==1 & d$time_months<=h
  observed_at_h<-d$time_months>=h
  usable<-early_event|observed_at_h
  weight<-ifelse(early_event,1/pmax(G(d$time_months-1e-8),.025),
           ifelse(observed_at_h,1/pmax(G(h),.025),0))
  target<-as.numeric(early_event)
  smooth_data_input<-data.frame(x=risk[usable],y=target[usable],
                                w=weight[usable])
  smooth<-stats::loess(y~x,data=smooth_data_input,weights=w,
                        span=.8,degree=1,control=stats::loess.control(surface="direct"))
  grid<-seq(max(0,quantile(risk,.01)),min(1,quantile(risk,.99)),length.out=125)
  yy<-as.numeric(stats::predict(smooth,newdata=data.frame(x=grid)))
  smooth_data<-data.frame(horizon_months=h,predicted_risk=grid,
                           observed_smoothed=pmin(1,pmax(0,yy)))
  write_aggregate(smooth_data,paste0("calibration_curve_",h,"m.csv"))
  png(file.path(FIG_DIR,paste0("Figure_calibration_",h,"m.png")),width=1600,height=1100,res=180)
  plot(grid,pmin(1,pmax(0,yy)),type="l",lwd=3,
       xlim=c(0,1),ylim=c(0,1),xlab="Predicted mortality",
       ylab="IPCW-smoothed observed mortality",
       main=paste0(h/12,"-year out-of-fold calibration"))
  abline(0,1,lty=2)
  mtext("Internal cross-fitting only; not external validation",side=1,line=4,cex=.8)
  dev.off()
}
# Effect-size sensitivity plot without manually transcribed manuscript numbers.
png(file.path(FIG_DIR,"Figure_tumor_size_sensitivity.png"),width=1700,height=1050,res=180)
r<-sens[sens$model=="A_clinical",]
plot(r$HR_per_10mm,seq_len(nrow(r)),xlim=range(c(r$lower95,r$upper95)),
     yaxt="n",ylab="",xlab="Adjusted hazard ratio per 10 mm",pch=19,
     main="Tumor-size exclusion sensitivity")
axis(2,at=seq_len(nrow(r)),labels=paste0("Up to ",r$max_mm," mm"),las=1)
segments(r$lower95,seq_len(nrow(r)),r$upper95,seq_len(nrow(r)),lwd=2)
abline(v=1,lty=2)
dev.off()
writeLines(c(paste0("R version: ",R.version.string),
  "Data: 2024 NCDB Bone/Joint PUF, licensed locally.",
  "All tables and figures generated from local fitted model objects.",
  "Cohort, imputation, and fitted-object RDS remain private.",
  "Internal cross-fitting and bootstraps are not external validation.",
  paste0("Cohort RDS md5: ",unname(tools::md5sum(file.path(PRIVATE_DIR,"cohort_INTERNAL_ONLY.rds"))))),
  file.path(AGG_DIR,"run_manifest.txt"))
cat("Aggregate reporting complete: ",AGG_DIR,"\n")
