######################################################################################################
# RPMM Classification Part I: Execute RPMM
######################################################################################################

library(RPMM)
library(matrixStats)
library(pheatmap)
library(doParallel); registerDoParallel(detectCores() - 1)

## Set directory based on host name:
setwd("Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/")

## Load beta value matrix:
load("121316_DCIS_betas2_for_RPMM.RData") #test node only

## Input matrix Y needs to be n specimens * J CpGs
## Here, BLC tree does not split unless J is approximately 13000. 
for(n in c(1.5*10^4)) {
  ## Matirx with J most variable CpGs:
  Y <- t(betas2[order(rowVars(betas2), decreasing=TRUE)[1:n], ])
  
  ## Specimen:
  assign(x=paste0("specimens_", n, "_top_CpGs"), value=rownames(Y), envir=.GlobalEnv)
  
  ## RPMM object:
  blc_rpmm <- blcTree(Y, verbose=1, maxlevel=1) #maxlevel 
  assign(x=paste0("rpmm_", n, "_top_CpGs"), value=blc_rpmm, envir=.GlobalEnv)
}

## Save RPMM objects to working directory
save(list=c("rpmm_15000_top_CpGs", "specimens_15000_top_CpGs"), file="122316_RPMM_objects_with_maxLev.RData")

#-----------------------------Images: Heatmaps & Dendrograms-----------------------------
{
  pdf(file=paste0("123116_RPMM_with_maxlevels.pdf"), paper="a4r")
  
  ## Plot heat maps:
  plot(rpmm_15000_top_CpGs); title("RPMM: top 15k most variable CpGs with 1 split")
  
  ## Plot tree dendrograms:
  labelFunction <- function(u, digits) table(as.character(specimens_15000_top_CpGs[u$index]))
  plotTree.blcTree(rpmm_15000_top_CpGs, labelFunction=labelFunction, cex=0.5)
  mtext(paste("DCIS: Class Membership for 1.5k most variable CpGs"), outer=TRUE, side=3)
  
  dev.off()
}

#-----------------------------Class Membership-----------------------------
write.csv(table(blcTreeLeafClasses(rpmm_15000_top_CpGs), specimens_15000_top_CpGs),
          file="122316_RPMM_with_1.5k_CpGs_preset_2_classes.csv", row.names=T, quote=F)

## Identify samples with same class membership:
odd <- function(x) x %% 2 != 0 #for indexing of paired specimens
summ15k <- as.data.frame(table(blcTreeLeafClasses(rpmm_15000_top_CpGs), specimens_15000_top_CpGs))
colnames(summ15k) <- c("Class", "Specimen", "Membership")
summ15k <- summ15k[summ15k$Membership == 1, ]
summ15k <- summ15k[order(summ15k$Specimen), ]

for(i in seq(nrow(summ15k))[odd(1:nrow(summ15k))]) {
  if(summ15k$Class[i] == summ15k$Class[i+1]) print(paste(summ15k$Specimen[i],summ15k$Specimen[i+1],sep=", "))
}

## Identify samples with different class membership:
for(i in seq(nrow(summ15k))[odd(1:nrow(summ15k))]) {
  if(summ15k$Class[i] != summ15k$Class[i+1]) print(paste(summ15k$Specimen[i],summ15k$Specimen[i+1],sep=", "))
}

######################################################################################################
# RPMM Classification Part II: Query RPMM objects
######################################################################################################

#-----------------------------Images: Heatmaps & Dendrograms-----------------------------
pdf(file=paste0("122116_RPMM_with_4_iter.pdf"), paper="a4r")

## Plot heat maps: 
plot(rpmm_15000_top_CpGs); title("RPMM: top 15k most variable CpGs")
plot(rpmm_17000_top_CpGs); title("RPMM: top 17k most variable CpGs")
plot(rpmm_19000_top_CpGs); title("RPMM: top 19k most variable CpGs")
plot(rpmm_20000_top_CpGs); title("RPMM: top 20k most variable CpGs")

## Plot tree dendrograms:
labelFunction <- function(u, digits) table(as.character(specimens_15000_top_CpGs[u$index]))
plotTree.blcTree(rpmm_15000_top_CpGs, labelFunction=labelFunction, cex=0.5)
mtext(paste("DCIS: Class Membership for 1.5k most variable CpGs"), outer=TRUE, side=3)

labelFunction <- function(u, digits) table(as.character(specimens_17000_top_CpGs[u$index]))
plotTree.blcTree(rpmm_17000_top_CpGs, labelFunction=labelFunction, cex=0.5)
mtext(paste("DCIS: Class Membership for 1.7k most variable CpGs"), outer=TRUE, side=3)

labelFunction <- function(u, digits) table(as.character(specimens_19000_top_CpGs[u$index]))
plotTree.blcTree(rpmm_19000_top_CpGs, labelFunction=labelFunction, cex=0.5)
mtext(paste("DCIS: Class Membership for 1.9k most variable CpGs"), outer=TRUE, side=3)

labelFunction <- function(u, digits) table(as.character(specimens_20000_top_CpGs[u$index]))
plotTree.blcTree(rpmm_20000_top_CpGs, labelFunction=labelFunction, cex=0.5)
mtext(paste("DCIS: Class Membership for 2.0k most variable CpGs"), outer=TRUE, side=3)

dev.off()

#-----------------------------CSV Export of Class Membership-----------------------------
write.csv(table(blcTreeLeafClasses(rpmm_15000_top_CpGs), specimens_15000_top_CpGs),
          file="122116_RPMM_with_1.5k_CpGs.csv", row.names=T, quote=F)
write.csv(table(blcTreeLeafClasses(rpmm_17000_top_CpGs), specimens_17000_top_CpGs),
          file="122116_RPMM_with_1.7k_CpGs.csv", row.names=T, quote=F)
write.csv(table(blcTreeLeafClasses(rpmm_19000_top_CpGs), specimens_19000_top_CpGs),
          file="122116_RPMM_with_1.9k_CpGs.csv", row.names=T, quote=F)
write.csv(table(blcTreeLeafClasses(rpmm_20000_top_CpGs), specimens_20000_top_CpGs),
          file="122116_RPMM_with_2k_CpGs.csv", row.names=T, quote=F)

#-----------------------------Summary Table for Same/Nearby/Distant Class Membership-----------------------------
odd <- function(x) x %% 2 != 0 #for indexing of paired specimens

## 4-class structure: 15k:
summ15k <- as.data.frame(table(blcTreeLeafClasses(rpmm_15000_top_CpGs), specimens_15000_top_CpGs))
colnames(summ15k) <- c("Class", "Specimen", "Membership")
summ15k <- summ15k[summ15k$Membership == 1, ]
summ15k <- summ15k[order(summ15k$Specimen), ]
## Samples with same class membership:
for(i in seq(nrow(summ15k))[odd(1:nrow(summ15k))]) {
  if(summ15k$Class[i] == summ15k$Class[i+1]) print(paste(summ15k$Specimen[i],summ15k$Specimen[i+1],sep=", "))
}

## 4-class structure: 17k:
summ17k <- as.data.frame(table(blcTreeLeafClasses(rpmm_15000_top_CpGs), specimens_15000_top_CpGs))
colnames(summ17k) <- c("Class", "Specimen", "Membership")
summ17k <- summ17k[summ17k$Membership == 1, ]
summ17k <- summ17k[order(summ17k$Specimen), ]
## Samples with same class membership:
for(i in seq(nrow(summ17k))[odd(1:nrow(summ17k))]) {
  if(summ17k$Class[i] == summ17k$Class[i+1]) print(paste(summ17k$Specimen[i],summ17k$Specimen[i+1],sep=", "))
}

## 2-class structure: 19k
summ19k <- as.data.frame(table(blcTreeLeafClasses(rpmm_19000_top_CpGs), specimens_19000_top_CpGs))
colnames(summ19k) <- c("Class", "Specimen", "Membership")
summ19k <- summ19k[summ19k$Membership == 1, ]
summ19k <- summ19k[order(summ19k$Specimen), ]
## Samples with same class membership:
for(i in seq(nrow(summ19k))[odd(1:nrow(summ19k))]) {
  if(summ19k$Class[i] == summ19k$Class[i+1]) print(paste(summ19k$Specimen[i],summ19k$Specimen[i+1],sep=", "))
}

## 2-class structure: 20k
summ20k <- as.data.frame(table(blcTreeLeafClasses(rpmm_19000_top_CpGs), specimens_19000_top_CpGs))
colnames(summ20k) <- c("Class", "Specimen", "Membership")
summ20k <- summ20k[summ20k$Membership == 1, ]
summ20k <- summ20k[order(summ20k$Specimen), ]
## Samples with same class membership:
for(i in seq(nrow(summ20k))[odd(1:nrow(summ20k))]) {
  if(summ20k$Class[i] == summ20k$Class[i+1]) print(paste(summ20k$Specimen[i],summ20k$Specimen[i+1],sep=", "))
}
