################################################################################################
# Loading and Pre-Processing of DCIS .idat files by minfi
# Script author: David Chen <youdinghuan.chen.gr@dartmouth.edu>
################################################################################################

library(minfi)
library(doParallel); registerDoParallel(detectCores() - 1)
rm(list=ls()) #ensure workspace cleared

parent.dir <- "/Users/DavidKevinChen/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper" #for ChAMP only
idat.dir <- "/Users/DavidKevinChen/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/DCIS-idats-only"
setwd(idat.dir)

# Make sure idat.dir has only .idat files and one .csv sample sheet.
# Specify samples of interest in the sample sheet.
# Here, the sample sheet was exported from R after sample selection.
sheet <- read.metharray.sheet(getwd()) #this line is optional for ChAMP but required for minfi
RGset <- read.metharray.exp(targets=sheet) # S4 object 27K, 450K or EPIC arrays are imported automatically 
names <- pData(RGset)$Sample_Name # or equivalent
groups <- pData(RGset)$Sample_Group # or equivalent

# Density plots for 1st QC:
par(mar=c(5,8,4,2)) #c(bottom, left, top, right)
densityBeanPlot(RGset, sampNames=names, sampGroups=groups, main="Raw intensities")
densityPlot(RGset, sampGroups=groups, main="Raw intensities")

# Generate the simple QC report and export as .PDF into a specified directory:
# Note: Sometimes the density or bean plots are not properly generated
qcReport(rgSet=RGset, pdf=paste0(parent.dir, "/qcReport.pdf"))

# Sample mix-ups by heat map:
library(pheatmap)
pheatmap(getSnpBeta(RGset))

# Convert to a MethylSet: 
Mset <- preprocessRaw(RGset)

# Fix outliers: 
Mset@annotation #proper annotation required
Mset2 <- minfiQC(Mset, fixOutliers=TRUE, verbose=TRUE)

# Check extreme outliers:
plotQC(Mset2$qc)

# Repeat the density plots:
# Use getBeta, since the object is Mset now:
densityPlot(getBeta(Mset), sampGroups=groups, main="Preprocessed raw data (Mset)")
densityBeanPlot(getBeta(Mset), sampNames = names, sampGroups=groups, main="Preprocessed raw data (Mset)")

# Convert into a mapToGenome object:
GRset <- mapToGenome(Mset)

# Normalization (Funnorm), background (methylumi.noob) & dye correction
# Update accordingly
GrsetFUNnorm <- preprocessFunnorm(RGset, nPCs=2, sex = NULL, bgCorr=TRUE,
                                  dyeCorr=TRUE, verbose=TRUE)

# More QC:
densityBeanPlot(getBeta(GrsetFUNnorm), sampNames=names, sampGroups=groups, main="Funnorm and background+dye correction")
densityPlot(getBeta(GrsetFUNnorm), sampGroups=groups, main="Funnorm and background+dye correction")

# Check phenotype:
# pData(GrsetFUNnorm)

# Quality filter:
detP <- detectionP(RGset) #get p-value
failed.05 <- detP > 0.05  #threshold
summary(failed.05)
dim(failed.05)

# List failed probes:
failedProbes <- rownames(failed.05)[rowMeans(failed.05) > 0.2] # list of probes that failed in more than 20% of the sample
sum(rowMeans(failed.05) > 0.2) # how many probes failed in more than 20% of samples 

# Remove low call-rate probes:
GrsetFUNnormfiltered <- GrsetFUNnorm[!rownames(GrsetFUNnorm) %in% failedProbes] 

# Add the QC in the phenotypes:
addQC(GrsetFUNnormfiltered, qc=Mset2$qc)

# Check dimensions of objects:
dim(GrsetFUNnorm)
dim(GrsetFUNnormfiltered) 

# Check filtered Funnorm quality:
densityBeanPlot(getBeta(GrsetFUNnormfiltered), sampNames=names, sampGroups=groups, main="Post-filtering Funnorm (exported)")
densityPlot(getBeta(GrsetFUNnormfiltered), sampGroups=groups, main="Post-filtering Funnorm (exported)") #OK to set `legend=F`

# Save the object:
save(list=ls(), file=paste(parent.dir, "DCIS-minfi-workspace.RData", sep="/"))
