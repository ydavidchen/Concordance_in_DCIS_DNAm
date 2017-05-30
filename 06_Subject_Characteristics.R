######################################################################################################
# Study Population Characteristics (Table 1)
# Script author: David (Youdinghuan) Chen <youinghuan.chen.gr@dartmouth.edu>
# Date: 12/08/2016
######################################################################################################

library(matrixStats)
library(tableone)
library(readxl)
library(doParallel); registerDoParallel(detectCores() - 1)

load("/Users/DavidKevinChen/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/113016-stratified-beta-matrices.RData")

## Check dimensions & matched samples:
dim(master)
dim(betas)
colnames(master)
table(master$bcid)
table(master$casepr)
all(master$key %in% colnames(betas))

## Remove levels that are not needed:
master <- droplevels(master)

## Update the values of certain varaibles:
master$tissue_descrip <- gsub("core_dcis", "Biopsy", master$tissue_descrip)
master$tissue_descrip <- gsub("exc_dcis", "Surgical", master$tissue_descrip)

master$Necrosis <- gsub(0, "Absent", master$Necrosis)
master$Necrosis <- gsub(1, "Present", master$Necrosis)

master$Sclerosis <- gsub(0, "Absent", master$Sclerosis)
master$Sclerosis <- gsub(1, "Present", master$Sclerosis)

master$Inflammation <- gsub(0, "Absent", master$Inflammation)
master$Inflammation <- gsub(1, "Present", master$Inflammation)

master$Calcs <- gsub(0, "Absent", master$Calcs)
master$Calcs <- gsub(1, "Present", master$Calcs)

master$mfamhx <- gsub(0, "Absent", master$mfamhx)
master$mfamhx <- gsub(1, "Present", master$mfamhx)

master$Grade <- gsub(1, "Low", master$Grade)
master$Grade <- gsub(2, "Intermediate", master$Grade)
master$Grade <- gsub(3, "High", master$Grade)

## Convert tumor size to numerical:
master$Tumor.size <- gsub(">", "", master$Tumor.size)
master$Tumor.size <- gsub("cm", "", master$Tumor.size)
master$Tumor.size <- gsub(" x.*", "", master$Tumor.size)
master$Tumor.size <- gsub(" largest.*", "", master$Tumor.size)
master$Tumor.size[master$Tumor.size == "\"Small focus\""] <- NA
master$Tumor.size <- as.numeric(master$Tumor.size)

## Make col.name abbreviations understandable:
colnames(master)[colnames(master) == "Calcs"] <- "Calcifications"
colnames(master)[colnames(master) == "mfamhx"] <- "Family history"
colnames(master)[colnames(master) == "agepr"] <- "Age"
colnames(master)[colnames(master) == "Tumor.size"] <- "Tumor size"
colnames(master)[colnames(master) == "tissue_descrip"] <- "Specimen"

## Assign categorical variables:
factorVars <- c("Grade", "Necrosis", "Sclerosis", "Inflammation", "Calcifications", "Family history", "Pattern")
all(factorVars %in% colnames(master)) #check

## Assign continuous random variables:
vars <- c(
  "Age", "Tumor size", #continuous
  "Grade", "Necrosis", "Sclerosis", "Inflammation", "Calcifications", "Family history", "Pattern" #discrete
  )
all(vars %in% colnames(master)) #check

TableOneForAll <- CreateTableOne(
  vars = vars, 
  strata = "Specimen", 
  data = master, 
  factorVars = factorVars,
  test = FALSE,
  smd = FALSE
)

x <- print(TableOneForAll, showAllLevels = TRUE)
write.csv(x, "~/repos/DCIS-Pathology-Paper-2016/Analysis-Results-Figures/120816D-Study-Population.csv", quote=F)
