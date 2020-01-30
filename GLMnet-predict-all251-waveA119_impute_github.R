# R script to take gene expression residuals for the 251 participants and impute their phenotypic data based on glmnet models built on the 119 Wave A participants

library(tidyr)
library(tidyverse)
library(glmnet)
library(data.table)
library(methods)

myvar <- "FEVPP"

args = commandArgs(trailingOnly=TRUE)

if (length(args)>0){
    myvar <- args[1]
  }

cores            <- as.integer(Sys.getenv("NCPUS"))
ParallelSapply <- function(...,mc.cores=cores){
    simplify2array(mclapply(...,mc.cores=mc.cores))
  }

timestamp()

myvar
cat("##Processing ",myvar,"\n")

# read in the gene expression residuals:
Res <- read.table("../counts/GC/gene_counts/residuals/residuals_251.txt", sep="\t", header=TRUE)
colnames(Res) <- gsub("[.]","-",colnames(Res))

#read in the corresponding covariates:
cv1 <- read.table("../covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19.txt", sep="\t", header=T, stringsAsFactors=F, quote='"')
# keep only samples selected for eQTL analysis and order by dbgap:
cv2 <- cv1[cv1$expanded_eQTL=="TRUE",]
cv <- cv2[order(cv2$dbgap.ID),]
rownames(cv) <- cv$barcode

# check that samples in GE data and covariate are in the same order:
identical(colnames(Res),cv$barcode)

# load the appropriate model:
load(paste0("/nfs/rprdata/ALOFT/AL/GLMnet/alpha0.1-LOO_119_rmXY/models/model_", myvar, "_.Rd"))

#Get the model parameters:
ind = which(mymodel$lambda==mymodel$lambda.min)
coefs <- mymodel$glmnet.fit$beta[,ind]
nzero <- coefs[which(coefs!=0)]
lambda <- mymodel$lambda[ind]

genes <- as.character(names(nzero))

intercept <- coef(mymodel)[1]
enstsweights <- coefs[genes]

# check if all genes in the model are found in residuals:
stopifnot(length(genes)==sum(genes %in% rownames(Res)))

## # make predictions on all 251:
vardata <- Res[genes,]
varpredictions <- enstsweights %*% as.matrix(vardata) + intercept
varpredictions <- t(data.frame(varpredictions))
rownames(varpredictions) <- gsub("[.]", "-", rownames(varpredictions))
identical(rownames(varpredictions), cv$barcode)

# check correlation of prediction:
myvar <- gsub("_1","_0",myvar)
cor.test(varpredictions[,1], as.numeric(cv[,myvar]))
corr <- cor.test(varpredictions[,1], as.numeric(cv[,myvar]))
# and plot:
pdf(paste0("./plots/model_fit/", myvar,".pdf"))
plot(x=varpredictions[,1],y=as.numeric(cv[,myvar]), main=myvar)
dev.off()
# and save corr
tab <- t(c(myvar, corr$estimate,corr$p.value))
write.table(tab, file="./GLMnet-correlations_all.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

# save the imputed values:
colnames(varpredictions) <- paste0(myvar, "pred")
write.table(varpredictions, file=paste0("./predictions/all/", myvar, "pred.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

sessionInfo()

### END
