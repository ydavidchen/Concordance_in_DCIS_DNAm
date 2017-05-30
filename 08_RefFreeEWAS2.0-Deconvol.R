######################################################################################################
# Cell-Type Deconvolution via Cluster Computing
# Script author: David (Youdinghuan) Chen <youinghuan.chen.gr@dartmouth.edu>
# Start date: 12/17/2016
# Notes:
# 1. Load the object for RPMM, which is the single-letter-coded (completely de-identified) beta matrix
# 2. Run this script on Discovery cluster
# 3. Upload this R script, companion PBS script, RData into globa/scratch/USER/WORKING_DIR
# 4. Reference: Dr. K Johnson's RefFreeEWAS deconvolution scripts for Komen normal breast methylation
######################################################################################################

## Load packages and data:
if("RefFreeEWAS" %in% installed.packages()) {library(RefFreeEWAS)} else {
    install.packages("RefFreeEWAS", repos='http://cran.us.r-project.org')}

library(RefFreeEWAS)
library(matrixStats)

setwd("/global/scratch/ydchen/DCIS_RefFreeEWAS")
load("121316_DCIS_betas2_for_RPMM.RData")

## Select the most variable 15K CpGs for deconvolution:
n <- 1.5 * 10^4
Y_shortened <- betas2[order(rowVars(betas2), decreasing=TRUE)[1:n], ]

## Iteratively compute cell type:
myCellMixArray <- RefFreeCellMixArray(Y_shortened, Klist=2:10, iters=25)

## Determine optimal number of cell types K by computing a bootstrapped vector
## accounting for cell mixture effect via Mu %*% t(Omega)
myDevianceBoots <- RefFreeCellMixArrayDevianceBoots(
  rfArray = myCellMixArray,
  Y = Y_shortened,
  R = 1000, #default: 5
  bootstrapIterations = 10 #default: 5
)

## Export relevant object(s) as RData into working directory:
save(list=c("myCellMixArray", "myDevianceBoots"), file="122516_DCIS_RefFreeEWAS.RData", compress=T)

sessionInfo()
