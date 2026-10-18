# Stage 08: histology and clinical M status; leaves locked models unchanged.
source(file.path("R","00_config.R"))
co <- readRDS(file.path(PRIVATE_DIR,"cohort_INTERNAL_ONLY.rds"))
mod <- readRDS(file.path(PRIVATE_DIR,"models_INTERNAL_ONLY.rds"))
el <- co$eligible; a <- mod$datasets$A_expanded
if(nrow(el)!=1562L || nrow(a)!=1243L || sum(a$event)!=310L)
  stop("Locked population mismatch")
out <- file.path(PRIVATE_DIR,"hist_mets_RESTRICTED_REVIEW")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
w <- function(x,p) utils::write.csv(x,file.path(out,p),row.names=FALSE)
r <- readr::read_csv(PUF_FILE,col_types=readr::cols(.default=readr::col_character()),
                    na=c("","NA","N/A"),show_col_types=FALSE)
names(r) <- toupper(trimws(sub("^\ufeff","",names(r))))
f <- function(n) if(n %in% names(r)) as.character(r[[n]]) else rep(NA_character_,nrow(r))
num <- function(x) suppressWarnings(as.numeric(x))
site <- toupper(gsub("[^A-Z0-9]","",f("PRIMARY_SITE")))
q <- !is.na(site) & grepl("^[0-9]{3}$",site);site[q] <- paste0("C",site[q])
q <- !is.na(site) & grepl("^C[0-9]{2}$",site);site[q] <- paste0(site[q],"0")
h <- f("HISTOLOGY");b <- f("BEHAVIOR")
m <- paste0(gsub("[^0-9]","",h),"/",gsub("[^0-9]","",b))
ix <- which(site=="C410" & m %in% c("9370/3","9371/3","9372/3") &
  is.finite(num(f("AGE"))) & num(f("AGE"))>=18 &
  is.finite(num(f("DX_LASTCONTACT_DEATH_MONTHS"))) &
  num(f("DX_LASTCONTACT_DEATH_MONTHS"))>=0 &
  f("PUF_VITAL_STATUS") %in% c("0","1"))
if(length(ix)!=nrow(el))stop("PUF eligibility mismatch")
if(!isTRUE(all.equal(num(f("AGE"))[ix],el$AGE_NUM)) ||
   !isTRUE(all.equal(num(f("DX_LASTCONTACT_DEATH_MONTHS"))[ix],el$time_months)) ||
   !identical(m[ix],as.character(el$HIST_BEH)))
  stop("PUF row-order mismatch")
r <- r[ix,,drop=FALSE]; f <- function(n) if(n %in% names(r))
  trimws(as.character(r[[n]])) else rep(NA_character_,nrow(r))
clinical <- c("AGE_NUM","tumor_size_10mm","SEX3","CDCC2")
ids <- which(is.finite(el$tumor_size_mm_raw) &
 el$tumor_size_mm_raw>=1 & el$tumor_size_mm_raw<=100 &
 stats::complete.cases(el[,clinical,drop=FALSE]))
if(length(ids)!=nrow(a))stop("Primary-cohort index mismatch")
for(v in c("AGE_NUM","YEAR_NUM","time_months","event","tumor_size_mm_raw"))
 if(!isTRUE(all.equal(el[[v]][ids],a[[v]],check.attributes=FALSE)))
   stop("Primary row order mismatch: ",v)
if(!all(c("TNM_CLIN_M","AJCC_TNM_CLIN_M") %in% names(r)))
  stop("Missing clinical M era variables")
year <- el$YEAR_NUM
clin <- ifelse(year>=2018,f("AJCC_TNM_CLIN_M"),f("TNM_CLIN_M"))
recode <- function(x) {
 z <- toupper(gsub("[[:space:]]+","",trimws(as.character(x))))
 z[is.na(z)] <- "";z <- sub("^(C|P)?M","",z)
 out <- rep("Unknown",length(z))
 out[z=="0"] <- "M0"
 out[grepl("^1[A-D]?$",z)] <- "M1"
 out
}
stopifnot(identical(recode(c(NA,"","MX","M0","cM1","M1a")),
  c("Unknown","Unknown","Unknown","M0","M1","M1")))
cm <- recode(clin)
fields <- c("METS_AT_DX_BONE","METS_AT_DX_BRAIN","METS_AT_DX_LIVER",
 "METS_AT_DX_LUNG","METS_AT_DX_DISTANT_LN","METS_AT_DX_OTHER")
available <- intersect(fields,names(r))
sm <- rep("Unknown",nrow(el))
if(length(available)){
 z <- sapply(available,f);if(is.null(dim(z)))z <- matrix(z,ncol=1)
 pos <- z=="1";pos[is.na(pos)] <- FALSE
 if("METS_AT_DX_OTHER" %in% colnames(z))
   pos[,"METS_AT_DX_OTHER"] <- z[,"METS_AT_DX_OTHER"] %in% c("1","2")
 yes <- rowSums(pos)>0
 zero <- if(length(available)==length(fields))
   rowSums(is.na(z)|z!="0")==0 else rep(FALSE,nrow(el))
 sm[year>=2016 & yes] <- "M1"
 sm[year>=2016 & zero] <- "M0"
}
conflict <- (cm=="M0" & sm=="M1")|(cm=="M1" & sm=="M0")
met <- rep("Unknown",length(cm))
met[cm=="M0"|sm=="M0"] <- "M0"
met[cm=="M1"|sm=="M1"] <- "M1"
met[conflict] <- "Conflicting"
stopifnot(!anyNA(met))
a$histology_code <- factor(el$HIST_BEH[ids],levels=c("9370/3","9371/3","9372/3"))
a$mets_status <- factor(met[ids],levels=c("M0","M1","Unknown","Conflicting"))
base <- stats::coef(mod$fits$A_expanded)
fit <- function(d,f) survival::coxph(f,data=d,ties="efron",singular.ok=FALSE)
refit <- function(select,name) {
 d <- droplevels(a[select,,drop=FALSE])
 if(nrow(d)<100L || sum(d$event)<50L)stop("Too sparse: ",name)
 b <- stats::coef(fit(d,mod$formulas$A_expanded))
 if(!identical(names(b),names(base)))stop("Parameterization mismatch")
 data.frame(sensitivity=name,term=names(base),HR=exp(unname(b)),
  primary_HR=exp(unname(base)),ratio_HR_to_primary=exp(unname(b-base)))
}
changes <- rbind(
 refit(as.character(a$histology_code)!="9372/3","Exclude 9372/3"),
 refit(!as.character(a$mets_status)%in%c("M1","Conflicting"),
       "Exclude documented M1 and conflicting cases"),
 refit(as.character(a$mets_status)=="M0","Restrict to documented M0"))
w(changes,"08_PRIVATE_sensitivity_coefficient_comparison.csv")
d <- droplevels(a[as.character(a$histology_code)%in%c("9370/3","9371/3"),])
ct <- table(d$histology_code)
ev <- tapply(d$event,d$histology_code,sum)
if(length(ct)!=2L||any(ct<20L)||any(ev<10L))
 stop("Two-code histology comparison not estimable")
d$histology_code <- factor(as.character(d$histology_code),
   levels=c("9370/3","9371/3"))
f0 <- fit(d,mod$formulas$A_expanded)
f1 <- fit(d,update(mod$formulas$A_expanded,.~.+histology_code))
k <- grep("^histology_code",names(stats::coef(f1)))
if(length(k)!=1L)stop("Unexpected histology term")
b <- unname(stats::coef(f1)[k]);s <- sqrt(diag(stats::vcov(f1)))[k]
hist <- data.frame(comparison="9371/3 versus 9370/3",
 adjusted_HR=exp(b),low95=exp(b-1.96*s),high95=exp(b+1.96*s),
 LRT_p=stats::pchisq(max(0,2*(f1$loglik[2]-f0$loglik[2])),1,lower.tail=FALSE))
w(hist,"08_PRIVATE_two_code_adjusted_histology.csv")
known <- droplevels(a[as.character(a$mets_status)%in%c("M0","M1"),])
ct <- table(known$mets_status)
ev <- tapply(known$event,known$mets_status,sum)
if(length(ct)!=2L||any(ct<20L)||any(ev<10L)) {
  status <- "Not estimable: sparse metastatic events"
  result <- data.frame(status=status,
    comparison="M1 versus M0 in documented-status cohort")
} else {
  status <- "Exploratory adjusted Cox sensitivity"
  known$mets_status <- factor(as.character(known$mets_status),levels=c("M0","M1"))
  fm <- fit(known,update(mod$formulas$A_expanded,.~.+mets_status))
  result <- data.frame(status=status,
    comparison="M1 versus M0 in documented-status cohort")
}
w(result,"08_PRIVATE_adjusted_metastatic_status.csv")
writeLines(c("08 Histology / metastatic status: stage completed.",
 "Source records matched exactly to locked eligible and primary model cohorts.",
 "Histology adjusted two-code model: FIT",
 "Exclusion-of-9372/3 coefficient sensitivity: FIT.",
 "Exclusion-of-documented-M1/conflicts coefficient sensitivity: FIT.",
 "M0-only coefficient sensitivity: FIT.",
 paste("M1 versus M0 adjusted model:",status),
 "Private outputs require NCDB disclosure review.",
 "Locked Cox equation and bootstraps unchanged."),
 file.path(out,"08_STATUS_FOR_REVIEW_NO_SMALL_CELLS.txt"))
cat("Stage 08 completed. Private results at ",out,"\n")
