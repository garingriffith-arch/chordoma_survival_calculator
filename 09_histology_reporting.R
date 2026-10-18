# Stage 09. Build S9 from the protected Stage 08 estimates.
source(file.path("R","00_config.R"))
dir08 <- file.path(PRIVATE_DIR,"hist_mets_RESTRICTED_REVIEW")
getcsv <- function(name) {
 path <- file.path(dir08,name)
 if(!file.exists(path))stop("Run Stage 08 first: ",name)
 utils::read.csv(path,stringsAsFactors=FALSE)
}
h <- getcsv("08_PRIVATE_two_code_adjusted_histology.csv")
m <- getcsv("08_PRIVATE_adjusted_metastatic_status.csv")
b <- getcsv("08_PRIVATE_sensitivity_coefficient_comparison.csv")
terms <- c("AGE_NUM","tumor_size_10mm","SEX3Female","CDCC21","CDCC22+")
scenarios <- c("Exclude 9372/3","Exclude documented M1 and conflicting cases",
               "Restrict to documented M0")
if(nrow(h)!=1L||nrow(b)!=15L||!setequal(unique(b$term),terms))
 stop("Unexpected Stage 08 table contents")
if(!grepl("Not estimable",m$status[1],fixed=TRUE))
 stop("Review M1 estimate before producing manuscript table")
fmt <- function(z)sprintf("%.3f",as.numeric(z))
A <- data.frame(Analysis=c("9371/3 vs 9370/3","9372/3 separate model",
 "M1 vs M0 adjusted comparison"),
 HR=c(fmt(h$adjusted_HR),"Not estimated","Not estimated"),
 CI=c(paste0(fmt(h$low95)," to ",fmt(h$high95)),"—","—"),
 Note=c(paste0("Likelihood-ratio p = ",fmt(h$LRT_p)),
  "Very sparse category","Sparse metastatic events"))
labels <- c("Age, per year","Size, per 10 mm","Female vs male",
 "Charlson-Deyo 1 vs 0","Charlson-Deyo >=2 vs 0")
B <- do.call(rbind,lapply(seq_along(terms),function(i){
 x<-b[b$term==terms[i],,drop=FALSE]
 if(nrow(x)!=3L)stop("Missing coefficient sensitivity: ",terms[i])
 v<-vapply(scenarios,function(s){
  a<-x$HR[x$sensitivity==s]
  if(length(a)!=1L)stop("Missing scenario ",s)
  a
 },numeric(1))
 data.frame(Predictor=labels[i],Primary=fmt(x$primary_HR[1]),
  Exclude_9372=fmt(v[1]),Exclude_M1=fmt(v[2]),Known_M0=fmt(v[3]))
}))
dest<-file.path(PRIVATE_DIR,"manuscript_exhibits_DISCLOSURE_REVIEW","tables")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(A,file.path(dest,"Table_S9a_Histology_Mets.csv"),row.names=FALSE)
utils::write.csv(B,file.path(dest,"Table_S9b_Coefficient_Stability.csv"),row.names=FALSE)

if(requireNamespace("officer",quietly=TRUE) && requireNamespace("flextable",quietly=TRUE)) {
  doc<-officer::read_docx()
  doc<-officer::body_add_par(doc,"Supplementary Table S9. Histology and baseline metastatic-status sensitivity",style="heading 1")
  doc<-officer::body_add_par(doc,"Panel A. Exploratory adjusted comparisons",style="heading 2")
  doc<-flextable::body_add_flextable(doc,flextable::autofit(flextable::flextable(A)))
  doc<-officer::body_add_par(doc,"Panel B. Stability of primary-model coefficients",style="heading 2")
  doc<-flextable::body_add_flextable(doc,flextable::autofit(flextable::flextable(B)))
  doc<-officer::body_add_par(doc,paste(
    "Sensitivity models retain the same five clinical predictors.",
    "These are associations, not additional validated calculators.",
    "Rare-category counts remain withheld; M0 restriction changes the population."))
  print(doc,target=file.path(dest,"Table_S9_Histology_Mets_Sensitivity.docx"))
}

cat("Stage 09 complete; private S9 CSVs in ",dest,"\n")
