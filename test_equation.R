# No patient-level data or Shiny packages required for this equation check.
stopifnot(file.exists("model/coefficients.csv"),
          file.exists("model/reference_survival.csv"))
b <- read.csv("model/coefficients.csv",check.names=FALSE)
s <- read.csv("model/reference_survival.csv",check.names=FALSE)
required_terms <- c("AGE_NUM","tumor_size_10mm","SEX3Female","CDCC21","CDCC22+")
stopifnot(identical(b$term,required_terms),
          identical(as.integer(s$horizon_months),c(60L,120L)),
          all(is.finite(b$beta)),all(is.finite(b$reference_design_value)),
          all(s$reference_survival>0 & s$reference_survival<1))
risk <- function(age,mm,sex,cdcc) {
  stopifnot(is.finite(age),age>=18,age<=110,
            is.finite(mm),mm>=1,mm<=100,
            sex%in%c("Male","Female"),cdcc%in%c("0","1","2+"))
  x <- c(age,mm/10,as.integer(sex=="Female"),
         as.integer(cdcc=="1"),as.integer(cdcc=="2+"))
  eta <- sum(b$beta*(x-b$reference_design_value))
  -expm1(log(s$reference_survival)*exp(eta))
}
ref <- risk(52,32,"Male","0")
stopifnot(max(abs(ref-(1-s$reference_survival)))<1e-12,
          all(diff(ref)>0),
          all(risk(62,32,"Male","0")>ref),
          all(risk(52,50,"Male","0")>ref),
          all(risk(52,32,"Male","1")>ref),
          all(risk(52,32,"Male","2+")>ref))
cat("Locked 2024 equation tests passed. Reference 5-/10-year mortality:",
    paste(sprintf("%.5f%%",100*ref),collapse="; "), "\n")
