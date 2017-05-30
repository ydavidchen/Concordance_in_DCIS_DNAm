################################################################################################
# Loading and pre-processing of TCGA .idat files by minfi
# Script author: David Chen <youdinghuan.chen.gr@dartmouth.edu>
# Initial processing script date: 11/28/2016
# Date adpated to TCGA data set: 03/16/2017
################################################################################################

library(minfi)
library(pheatmap)
library(splitstackshape)
library(doParallel); registerDoParallel(detectCores() - 1)

#------------------------------------------Assemble a CSV sample sheet------------------------------------------
load("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/RNAseq&450kSampleAnnot.Rdata")

parent.dir <- "~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/" 
idat.dir <- "~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/Unpacked_450k_IDATs/"
setwd(idat.dir)

all(HM450_meta_filtered$file_name %in% list.files()) #checkpoint: TRUE required
all(list.files() %in% HM450_meta_filtered$file_name) 
## [1] FALSE
## Note: more sample idats were downloaded
## so we need to use the sample sheet to select only the relevant ones (51 pairs of idats)
nrow(HM450_meta_filtered) / 2 #number of TCGA subjects
sum(duplicated(substr(HM450_meta_filtered$cases, 1, 16))) #checkpoint: 51 expected
## [1] 51

sample_sheet <- data.frame(
  Sample_Name  = unique(substr(HM450_meta_filtered$cases, 1, 16)),
  Sample_Well  = NA,
  Sample_Plate = NA,
  Sample_Group = NA,
  Pool_ID      = NA,
  Sentrix_ID   = unique(substr(HM450_meta_filtered$file_name, 1, 17))
)
sample_sheet <- cSplit(sample_sheet, splitCols=ncol(sample_sheet), sep="_", drop=TRUE)
class(sample_sheet) <- "data.frame" #required
colnames(sample_sheet)[(ncol(sample_sheet)-1):ncol(sample_sheet)] <- c("Sentrix_ID", "Sentrix_Position")
dim(sample_sheet)
## [1] 51  7
View(sample_sheet)
# write.csv(sample_sheet, file="TCGA_idats_minfi_samp_sheets.csv", row.names=F, quote=F)

#------------------------------------------Pre-processing & QC------------------------------------------
# Make sure idat.dir has only .idat files and one .csv sample sheet.
# Specify samples of interest in the sample sheet.
# Here, the sample sheet was exported from R after sample selection.
getwd()
sheet <- read.metharray.sheet(getwd()) #this line is optional for ChAMP but required for minfi
RGset <- read.metharray.exp(targets=sheet) # S4 object 27K, 450K or EPIC arrays are imported automatically 
names <- pData(RGset)$Sample_Name # or equivalent
# groups <- pData(RGset)$Sample_Group # or equivalent

# Density plots for 1st QC:

densityBeanPlot(RGset, sampNames=names, main="Raw intensities")
par(mar=c(5,4,4,2)) #default
densityPlot(RGset, main="Raw intensities")

# Generate the simple QC report and export as .PDF into a specified directory:
# Note: Sometimes the density or bean plots are not properly generated
qcReport(rgSet=RGset, pdf=paste0(parent.dir, "/20170316_TCGA_stage1ERpos_qcReport.pdf"))

# Sample mix-ups by heat map:
pheatmap(getSnpBeta(RGset), border_color=NA)

# Convert to a MethylSet: 
Mset <- preprocessRaw(RGset)

# Fix outliers: 
Mset@annotation #proper annotation required
## array                     annotation 
## "IlluminaHumanMethylation450k"                  "ilmn12.hg19" 
Mset2 <- minfiQC(Mset, fixOutliers=TRUE, verbose=TRUE)

# Check extreme outliers:
plotQC(Mset2$qc)

# Repeat the density plots:
# Use getBeta, since the object is Mset now:

par(mar=c(5,15,4,2))
densityBeanPlot(getBeta(Mset), sampNames = names, main="Preprocessed raw data (Mset)")
par(mar=c(5,4,4,2))
densityPlot(getBeta(Mset), main="Preprocessed raw data (Mset)")


# Convert into a mapToGenome object:
GRset <- mapToGenome(Mset)

# Normalization (Funnorm), background (methylumi.noob) & dye correction
# Update accordingly:
GrsetFUNnorm <- preprocessFunnorm(RGset, nPCs=2, sex=NULL, bgCorr=TRUE, dyeCorr=TRUE, verbose=TRUE)

# More QC:
par(mar=c(5,4,4,2))
densityPlot(getBeta(GrsetFUNnorm), main="Funnorm and background+dye correction")
par(mar=c(5,15,4,2))
densityBeanPlot(getBeta(GrsetFUNnorm), sampNames=names, main="Funnorm and background+dye correction")

## Temporarily extract beta values for confirming samples to drop:
tempBeta <- getBeta(GrsetFUNnorm)
densityBeanPlot(getBeta(GrsetFUNnorm), main="Funnorm and background+dye correction") #see probe name
bad_samples <- c("6042316157_R02C02", "6057833142_R01C02") #found by examining densityBeanPlot
which(colnames(tempBeta) %in% bad_samples)
densityPlot(tempBeta[ ,which(colnames(tempBeta) %in% bad_samples)], main="Samples considered for removal")
legend("topleft", legend=bad_samples)
densityPlot(tempBeta[ ,-which(colnames(tempBeta) %in% bad_samples)], main="Funnorm and background+dye correction, after bad samples removed")
rm(list=c("tempBeta", "bad_samples"))

## Check phenotype:
# pData(GrsetFUNnorm)

## Quality filter:
detP <- detectionP(RGset) #get p-value
failed.05 <- detP > 0.05  #threshold
summary(failed.05)
dim(failed.05)

## List failed probes:
failedProbes <- rownames(failed.05)[rowMeans(failed.05) > 0.20] # list of probes that failed in more than 20% of the sample
sum(rowMeans(failed.05) > 0.20) # num. probes failed in more than 20% of samples 
## [1] 1176

## Remove low call-rate probes:
GrsetFUNnormfiltered <- GrsetFUNnorm[!rownames(GrsetFUNnorm) %in% failedProbes] 

## Add the QC in the phenotypes:
addQC(GrsetFUNnormfiltered, qc=Mset2$qc)

## Check dimensions of objects:
dim(GrsetFUNnorm)
dim(GrsetFUNnormfiltered) 

## Check filtered Funnorm quality:
densityPlot(getBeta(GrsetFUNnormfiltered), main="Post-filtering Funnorm (to be exported for analysis)") #OK to set `legend=F`
densityBeanPlot(getBeta(GrsetFUNnormfiltered), sampNames=names, main="Post-filtering Funnorm (to be exported for analysis)")

#------------------------------------------Export------------------------------------------
## Save workspace:
rm(list=c("HM450_meta", "RNAseq_meta"))
save(list=ls(), file=paste(parent.dir, "20170316_GDC_TCGA_stage1ERpos_completeWorkspace.RData", sep="/"), compress=TRUE)

## Retrieve beta values
betas_stage1ERpos <- getBeta(GrsetFUNnormfiltered)
save(list=c("betas_stage1ERpos", "sample_sheet", "HM450_meta_filtered", "RNAseq_meta_filtered"), 
     file=paste(parent.dir, "20170316_TCGA_stage1ERpos_betas.RData", sep="/"))
