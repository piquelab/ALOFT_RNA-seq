# R script to take gene expression residuals and train elastic net models for phenotype data using glmnet

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

cores <- as.integer(Sys.getenv("NCPUS"))
ParallelSapply <- function(...,mc.cores=cores){
    simplify2array(mclapply(...,mc.cores=mc.cores))
  }

timestamp()

myvar
cat("##Processing ",myvar,"\n")

# read in the residuals:
load("/nfs/rprdata/ALOFT/AL/counts/GC/gene_counts/AL-residuals_waveA_119_expanded_baseline.Rd")

#read in the corresponding covariates:
load("/nfs/rprdata/ALOFT/AL/counts/GC/gene_counts/AL-covariates_waveA_119.Rd")

# make sure the cv and GE data are in the same order:
stopifnot(identical(colnames(data), cv$barcode))

# keep the full covariate and GE data before subsetting:
cvfull <- cv
datafull <- data

# if Ear variable, subset to :
EAR <- c("EYell_1","EPMres_1","ECPA_1","ECNA_1")
if(myvar %in% EAR){
  cv <-cv[-which(as.numeric(cv$Enumtf_1)<30),]
}

# subset data to only samples that have y:
cv <- cv[!is.na(cv[,myvar]),]
cv[,myvar] <- as.numeric(cv[,myvar])
# subset data make sure they're in the same order:
data <- data[,colnames(data) %in% cv$barcode]
identical(colnames(data),cv$barcode)

# train glmnet model:
y <- cv[,myvar]
N <- length(y)
x <- t(data)
colnames(x) <- rownames(data)
set.seed(12)
mymodel <- cv.glmnet(x,y,family="gaussian", alpha=0.1, type.measure="mse", nfold=length(y))

#Get lambda from model
ind = which(mymodel$lambda==mymodel$lambda.min)
coefs <- mymodel$glmnet.fit$beta[,ind]
nzero <- coefs[which(coefs!=0)]
lambda <- mymodel$lambda[ind]

ensts <- names(nzero)
genes <- as.character(names(nzero))

intercept <- coef(mymodel)[1]
enstsweights <- coefs[ensts]

cvm <- mymodel$cvm[ind]
cvm_null <- mymodel$cvm[1]
improve <- 1-cvm/cvm_null

# save the glmnet model:
save(mymodel,file=paste0("./models/model_", myvar, "_.Rd"))

# make predictions on all, including dropouts and check correlation:
vardata <- datafull[ensts,]
varpredictions <- enstsweights %*% as.matrix(vardata) + intercept
varpredictions <- t(data.frame(varpredictions))
rownames(varpredictions) <- gsub("[.]", "-", rownames(varpredictions))
identical(rownames(varpredictions), cvfull$barcode)
cor.test(varpredictions[,1], as.numeric(cvfull[,myvar]))
corr <- cor.test(varpredictions[,1], as.numeric(cvfull[,myvar]))

# add official variable name and save correlation:
labels <- read.table("./ALOFT_variables_for_GE_analysis_labels.txt", sep="\t", header=TRUE)
rownames(labels) <- as.character(labels[,1])
labels[,2] <- as.character(labels[,2])
tab <- t(c(labels[myvar,2], myvar, corr$estimate,corr$p.value, improve, N))
write.table(tab, file="./GLMnet-correlations.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

cat("Predictive:\t", myvar, "\t", length(genes), "\t")

# save the names of genes that went into the model:
if (length(genes)>0){
write.table(genes, file=paste0("./predictive_genes/waveA/", myvar, "-pred-genes.txt"), sep="/n", quote=FALSE, row.names=FALSE, col.names=FALSE)
# and the predictions:
write.table(varpredictions, file=paste0("./predictions/waveA/", myvar, "pred.txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")

# plot the correlation between actual and predicted:
pdf(paste0("./plots/waveA/", myvar,"_expanded_baseline_res_corr.pdf"))
plot(x=cvfull[,myvar], y=varpredictions[,1], xlab=paste0(myvar), ylab=paste0("GLMnet-predicted ", myvar), main=paste0("Correlation: ", myvar, " vs. GLMnet-predicted ", myvar, " on ", length(enstsweights), " genes"))
lmfit <- lm(varpredictions[,1]~as.numeric(cvfull[,myvar]))
abline(a=lmfit$coefficients[1], b=lmfit$coefficients[2], col="red")
dev.off()
# simple plot:
pdf(paste0("./plots/waveA/GLMnet-model-", myvar, ".pdf"))
plot(mymodel, main=paste0(length(genes), " genes best predict ", myvar))
dev.off()

sessionInfo()

### END

