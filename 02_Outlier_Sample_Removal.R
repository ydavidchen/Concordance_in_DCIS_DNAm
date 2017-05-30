####################################################################################################
# Checking & removing 2 pairs of samples classified as poor-performing outliers
# Script author: David Chen <youinghuan.chen.gr@dartmouth.edu>
# Date: 12/07/2016
####################################################################################################

# Load the complete minfi workspace (2 .RData objects, including the renamed):
load("/Users/DavidKevinChen/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/DCIS-minfi-workspace.RData")
load("/Users/DavidKevinChen/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/113016-stratified-beta-matrices-OLD.RData")

# Check extreme outliers:
plotQC(Mset2$qc)

# Query the names of the samples:
failed_samples <- c()
for(i in c(22,25)) failed_samples <- c(failed_samples, paste(sheet$Slide[i], sheet$Array[i], sep="_"))
failed_samples

# Remove failed samples and their companion pairs:
failed_bcid <- master$bcid[master$key %in% failed_samples]
failed_bcid

master <- master[-which(master$bcid %in% failed_bcid), ]
table(as.character(master$bcid))
length(unique(master$bcid))

# Select beta values:
betas <- getBeta(GrsetFUNnormfiltered)
betas <- betas[ , -which(!colnames(betas) %in% master$key)]

# Keys for subsetting:
key.corePun <- master[master$tissue_descrip == "core_dcis", ]
key.surgExc <- master[master$tissue_descrip == "exc_dcis", ]

# IMPORTANT: Make sure unique (identifying) columns are pairwise identical:
identical(key.corePun$bcid, key.surgExc$bcid) ##TRUE expected
identical(key.corePun$casepr, key.surgExc$casepr) ##TRUE expected

# Subsetting:
betas.corePun <- betas[ , key.corePun$key]
betas.surgExc <- betas[ , key.surgExc$key]

# IMPORTANT: Check dimensions of all objects to export
for(x in list(betas, betas.corePun, betas.surgExc, covars, manifest, master)) print(dim(x))

# Save as before: 
save(list=c("betas", "betas.corePun", "betas.surgExc", "covars", "manifest", "master"), 
     file = "~/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/113016-stratified-beta-matrices.RData")
