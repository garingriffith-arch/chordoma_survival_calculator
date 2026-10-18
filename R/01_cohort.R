# 01. Import, cohort definition, baseline coding, and audit.
# The NCDB file is read only. Individual records are never written to aggregate/.
source(file.path("R","00_config.R"))
set.seed(SEED)
raw <- readr::read_csv(PUF_FILE, col_types=readr::cols(.default=readr::col_character()),
                       na=c("","NA","N/A"),show_col_types=FALSE)
names(raw) <- toupper(trimws(sub("^\ufeff","",names(raw))))
if(anyDuplicated(names(raw))) stop("Duplicated PUF header fields")
required <- c("PRIMARY_SITE","HISTOLOGY","BEHAVIOR","YEAR_OF_DIAGNOSIS",
              "AGE","SEX","RACE","SPANISH_HISPANIC_ORIGIN","INSURANCE_STATUS",
              "CDCC_TOTAL_BEST","DX_LASTCONTACT_DEATH_MONTHS","PUF_VITAL_STATUS")
missing_fields<-setdiff(required,names(raw))
if(length(missing_fields)) stop("Required PUF variables missing: ",paste(missing_fields,collapse=", "))
n_all<-nrow(raw)
field<-function(x) if(x %in% names(raw)) trimws(as.character(raw[[x]])) else rep(NA_character_,n_all)
number<-function(x) suppressWarnings(as.numeric(x))
coalesce_era<-function(choices,accepted=NULL) {
  result<-rep(NA_character_,n_all)
  for(nm in choices) {
    v<-field(nm); ok<-is.na(result)&!is.na(v)&nzchar(v)
    if(!is.null(accepted)) ok<-ok & v%in%accepted
    result[ok]<-v[ok]
  }
  result
}
site<-toupper(gsub("[^A-Z0-9]","",field("PRIMARY_SITE")))
site[!is.na(site)&grepl("^[0-9]{3}$",site)] <- paste0("C",site[!is.na(site)&grepl("^[0-9]{3}$",site)])
site[!is.na(site)&grepl("^C[0-9]{2}$",site)] <- paste0(site[!is.na(site)&grepl("^C[0-9]{2}$",site)],"0")
hist0<-field("HISTOLOGY"); beh0<-field("BEHAVIOR")
hist<-ifelse(grepl("^[0-9]{4}/[0-9]$",hist0),substr(hist0,1,4),gsub("[^0-9]","",hist0))
beh<-ifelse(is.na(beh0)&grepl("^[0-9]{4}/[0-9]$",hist0),
            substr(hist0,6,6),gsub("[^0-9]","",beh0))
morph<-paste0(hist,"/",beh)
age<-number(field("AGE")); year<-number(field("YEAR_OF_DIAGNOSIS"))
time<-number(field("DX_LASTCONTACT_DEATH_MONTHS"))
event<-ifelse(field("PUF_VITAL_STATUS")=="0",1L,
              ifelse(field("PUF_VITAL_STATUS")=="1",0L,NA_integer_))
parse_size<-function(x) {
  z<-number(x); ifelse(is.finite(z)&z==floor(z)&z>=1&z<=988,z,NA_real_)
}
old_size<-coalesce_era(c("CS_TUMOR_SIZE","TUMOR_SIZE"))
modern<-field("TUMOR_SIZE_SUMMARY_2016")
is_modern<-!is.na(year)&year>=2016
original_size<-ifelse(is_modern,modern,old_size)
source_size<-ifelse(is_modern,"TUMOR_SIZE_SUMMARY_2016","CS_TUMOR_SIZE_or_TUMOR_SIZE")
# An out-of-era field is used only when that era's source column is absent entirely.
if(!("TUMOR_SIZE_SUMMARY_2016"%in%names(raw))) original_size[is_modern]<-old_size[is_modern]
if(!any(c("CS_TUMOR_SIZE","TUMOR_SIZE")%in%names(raw))) original_size[!is_modern]<-modern[!is_modern]
size<-parse_size(original_size)
sex<-ifelse(field("SEX")=="1","Male",ifelse(field("SEX")=="2","Female",NA_character_))
race<-ifelse(field("RACE")%in%c("01","1"),"White",
      ifelse(field("RACE")%in%c("02","2"),"Black",
      ifelse(field("RACE")%in%c("99","",NA),NA_character_,"Other")))
hisp<-ifelse(field("SPANISH_HISPANIC_ORIGIN")=="0","Non-Hispanic",
       ifelse(field("SPANISH_HISPANIC_ORIGIN")%in%as.character(1:8),"Hispanic",NA_character_))
ins<-ifelse(field("INSURANCE_STATUS")=="1","Private",
     ifelse(field("INSURANCE_STATUS")=="3","Medicare",
     ifelse(field("INSURANCE_STATUS")%in%c("2","4"),"Medicaid/OtherGov",
     ifelse(field("INSURANCE_STATUS")=="0","Uninsured",NA_character_))))
income<-coalesce_era(c("MED_INC_QUAR_2020","MED_INC_QUAR_2016",
                       "MED_INC_QUAR_12","MED_INC_QUAR_00"),accepted=as.character(1:4))
cdcc<-number(field("CDCC_TOTAL_BEST"))
cc<-ifelse(is.na(cdcc)|!(cdcc%in%0:9),NA_character_,
           ifelse(cdcc>=2,"2+",as.character(cdcc)))
dat<-data.frame(
  YEAR_NUM=year,AGE_NUM=age,PRIMARY_SITE_STD=site,HIST_BEH=morph,
  time_months=time,event=event,
  tumor_size_mm_raw=size,tumor_size_code=original_size,
  tumor_size_source=source_size,
  SEX3=factor(sex,levels=c("Male","Female")),
  CDCC2=factor(cc,levels=c("0","1","2+")),
  RACE3=factor(race,levels=c("White","Black","Other")),
  HISP2=factor(hisp,levels=c("Non-Hispanic","Hispanic")),
  INS3=factor(ins,levels=c("Private","Medicare","Medicaid/OtherGov","Uninsured")),
  MED_INC_QUAR=factor(income,levels=as.character(1:4)),
  stringsAsFactors=FALSE)
eligible<-dat[!is.na(site)&site=="C410"&morph%in%c("9370/3","9371/3","9372/3")&
             is.finite(age)&age>=18&is.finite(time)&time>=0&!is.na(event),,drop=FALSE]
eligible$tumor_size_10mm<-eligible$tumor_size_mm_raw/10
clinical<-c("AGE_NUM","tumor_size_10mm","SEX3","CDCC2")
social<-c("RACE3","HISP2","INS3","MED_INC_QUAR")
keep_size<-is.finite(eligible$tumor_size_mm_raw)&eligible$tumor_size_mm_raw>=1&
           eligible$tumor_size_mm_raw<=MAX_MM
A<-eligible[keep_size&complete.cases(eligible[,clinical,drop=FALSE]),,drop=FALSE]
D<-eligible[keep_size&complete.cases(eligible[,c(clinical,social),drop=FALSE]),,drop=FALSE]
A<-droplevels(A);D<-droplevels(D)
flow<-data.frame(
  stage=c("All bone/joint PUF rows","Survival-eligible adult C41.0 chordoma",
          "Complete clinical predictors, size 1–100 mm",
          "Complete clinical and social predictors, size 1–100 mm"),
  n=c(n_all,nrow(eligible),nrow(A),nrow(D)))
write_aggregate(flow,"cohort_flow.csv")
miss<-vapply(c(clinical,social),function(v)sum(is.na(eligible[[v]])),numeric(1))
write_aggregate(data.frame(variable=names(miss),missing=unname(miss),
                           denominator=nrow(eligible),
                           missing_percent=100*unname(miss)/nrow(eligible)),
                "predictor_missingness_LOCAL_DISCLOSURE_REVIEW.csv")
write_aggregate(data.frame(model=c("A_expanded","A_common","D_common"),
  n=c(nrow(A),nrow(D),nrow(D)),
  deaths=c(sum(A$event),sum(D$event),sum(D$event))),
  "model_cohorts.csv")
if(!QUICK && (nrow(A)!=1243L || sum(A$event)!=310L ||
               nrow(D)!=946L || sum(D$event)!=242L)) {
  stop("2024 cohort parity check FAILED. Verify field coding, PUF version, and exclusions.")
}
# Fitted objects embedding records and patient-level data must never be committed.
write_private(list(eligible=eligible,A=A,D=D,flow=flow),"cohort_INTERNAL_ONLY.rds")
cat("Cohort:",nrow(A),"primary,",nrow(D),"common-cohort patients.\n")
