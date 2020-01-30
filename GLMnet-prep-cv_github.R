# R script to prepare covariate data for training glmnet model (cleaning up variables)

library(tidyr)
library(tidyverse)
library(glmnet)
library(data.table)
library(methods)


timestamp()

#read in the covariates:
load("/nfs/rprdata/ALOFT/AL/counts/GC/gene_counts/AL-covariates_waveA_119.Rd")
##### 2/8/2019 convert pincme to numeric:
cv$pincme <- as.factor(cv$pincme)
levels(cv$pincme)[levels(cv$pincme)=="$0-$7,825 per year"] <- 0
levels(cv$pincme)[levels(cv$pincme)=="$7,826- $31,850 per year"] <- 1
levels(cv$pincme)[levels(cv$pincme)=="$31,851- $64,250 per year"] <- 2
levels(cv$pincme)[levels(cv$pincme)=="$64,251- $97,925 per year"] <- 3
levels(cv$pincme)[levels(cv$pincme)=="$97,926- $174,850 per year"] <- 4
cv$pincme <- as.numeric(cv$pincme)
# drop cssdh value > 16 (outliers):
cv$cssdh <- as.numeric(cv$cssdh)
cv[!is.na(cv$cssdh) & cv$cssdh>16, "cssdh"] <- NA
##### 2/10/2019 convert pedu to numeric:
cv$pedu <- as.factor(cv$pedu)
levels(cv$pedu)[levels(cv$pedu)=="No schooling completed"] <- 0
levels(cv$pedu)[levels(cv$pedu)=="9th grade"] <- 1
levels(cv$pedu)[levels(cv$pedu)=="10th grade"] <- 2
levels(cv$pedu)[levels(cv$pedu)=="11th grade"] <- 3
levels(cv$pedu)[levels(cv$pedu)=="12th grade, No diploma"] <- 4
levels(cv$pedu)[levels(cv$pedu)=="High school graduate- high school diploma or equivalent (for"] <- 5
levels(cv$pedu)[levels(cv$pedu)=="Some college credit, but less than 1 year"] <- 6
levels(cv$pedu)[levels(cv$pedu)=="1 or more years of college, no degree"] <- 7
levels(cv$pedu)[levels(cv$pedu)=="Associates degree (AA, AS)"] <- 8
levels(cv$pedu)[levels(cv$pedu)=="Bachelor's degree (BA, AB, BS)"] <- 9
levels(cv$pedu)[levels(cv$pedu)=="Master's degree (MA, MS, MEng, MEd, MSW, MBA)"] <- 10
levels(cv$pedu)[levels(cv$pedu)=="Doctorate degree (PhD, EdD)"] <- 11
cv$pedu <- as.numeric(cv$pedu)

# save the cleaned covariate:
save(cv,file="/nfs/rprdata/ALOFT/AL/counts/GC/gene_counts/AL-covariates_waveA_119_clean.Rd")
