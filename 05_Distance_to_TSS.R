######################################################################################################
# For Manuscript: Curves and Heat Matrices for Distance to TSS
# Script author: David (Youdinghuan) Chen <youinghuan.chen.gr@dartmouth.edu>
# Date: 12/08/2016
######################################################################################################

library(ggplot2)
library(reshape)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(matrixStats)
library(genomation)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(doParallel); registerDoParallel(detectCores() - 1)

load("/Users/DavidKevinChen/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/113016-stratified-beta-matrices.RData")

#--------------------------------------Heat Matrices for Summary Measure--------------------------------------
# Collapse the stratified beta values by median:
## IMPORTANT: the data column, `beta`, needs to be named identically across the two data matrices
betas.corePun <- data.frame(row.names=rownames(betas.corePun), beta=rowMedians(betas.corePun))
betas.surgExc <- data.frame(row.names=rownames(betas.surgExc), beta=rowMedians(betas.surgExc))

# Extract a +/- 2kb window around annotated TSS:
library(TxDb.Hsapiens.UCSC.hg19.knownGene) #annotation data package: txdb object
promoter_regions <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream=2000, downstream=2000)
head(promoter_regions)

# Generate score matrix & plot
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) #annotation data package
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot.450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot.450k <- as.data.frame(annot.450k@listData)
annot.450k <- data.frame(
  chr = annot.450k$chr,
  start = as.numeric(annot.450k$pos),
  end = as.numeric(annot.450k$pos),
  cgID = annot.450k$Name
)
rownames(annot.450k) <- annot.450k$cgID

# Make GRanges objects:
GR_corePun <- merge(annot.450k, betas.corePun, by="row.names")
GR_corePun <- makeGRangesFromDataFrame(GR_corePun, keep.extra.columns = TRUE)
GR_corePun@elementMetadata@listData$Row.names <- NULL
GR_corePun@elementMetadata@listData$cgID <- NULL
GR_corePun@elementMetadata@listData$UCSCgene <- NULL
head(GR_corePun)

GR_surgExc <- merge(annot.450k, betas.surgExc, by="row.names")
GR_surgExc <- makeGRangesFromDataFrame(GR_surgExc, keep.extra.columns = TRUE)
GR_surgExc@elementMetadata@listData$Row.names <- NULL
GR_surgExc@elementMetadata@listData$cgID <- NULL
GR_surgExc@elementMetadata@listData$UCSCgene <- NULL
head(GR_surgExc)

# Curves:
## `window` argument: pre-defined region around TSS or others of interest w/ EQUAL length
## `is.noCovNA` argument: required for CpG array data where coverage of bases is incomplete
sml <- ScoreMatrixList(
  targets = list(Biopsy = GR_corePun, Surgical = GR_surgExc), 
  windows = promoter_regions, 
  bin.num = 75,
  strand.aware = FALSE,
  weight.col = "beta",
  is.noCovNA = TRUE, 
  cores = detectCores()-1
)

plotMeta(
  sml, 
  line.col = brewer.pal(n = 9, name = "Set1"), 
  xcoords = c(-2000, 2000),
  main = "Matched DCIS specimens (n=13)",
  xlab = "Distance to TSS (bp)",
  ylab = "Median beta (fraction of methylated alleles)"
)
legend("bottomright", c("Biopsy", "Surgical"), lty=c(1, 1), lwd=c(2.5, 2.5), col=brewer.pal(n = 9, name = "Set1"))

# Heat Matrices
pdf("~/repos/DCIS-Pathology-Paper-2016/Analysis-Results-Figures/120916C-Distance-to-TSS-Heat-Matrices.PDF")
multiHeatMatrix(
  sml,
  col = colorRampPalette(colors = c("yellow", "blue"))(1024),
  xcoords = c(-2000, 2000), 
  matrix.main = c("Biopsy", "Surgical"), 
  legend.name = c("Median beta", "Median beta"),
  xlab = "Distance to TSS (bp)"
)
dev.off()

#--------------------------------------Individual Plots (For Supplemental)--------------------------------------
rm(list=ls())
load("/Users/DavidKevinChen/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/113016-stratified-beta-matrices.RData")

# Create a specimen key:
master$specimenKey <- paste(master$bcid, master$tissue_descrip, sep="-")

# Iteratively change the col.names of betas (`key`) to respective `specimenKey`:
# Match names:
betas <- betas[ , match(master$key, colnames(betas)) ] #reorder the columns
identical(colnames(betas), master$key) #TRUE expected
colnames(betas) <- master$specimenKey

# Letter code the paired samples iteratively
# Each pair of the sample will get the same name
master$hmID <- NA
for(i in 1:(nrow(master)/2) ) {
  x <- unique(as.character(master$bcid))[i]
  master$hmID[grep(x, master$bcid)] <- c(LETTERS[i], LETTERS[i])
}

# Reassign column names iteratively using the master$hmID column:
for(i in 1:(nrow(master)/2) ) {
  index <- grep(colnames(betas.corePun)[i], master$key) #search
  colnames(betas.corePun)[i] <- paste(master$hmID[index], "biopsy", sep=":") #reassign
}

# Reassign column names iteratively using the master$hmID column:
for(i in 1:(nrow(master)/2) ) {
  index <- grep(colnames(betas.surgExc)[i], master$key) #search
  colnames(betas.surgExc)[i] <- paste(master$hmID[index], "surgical", sep=":") #reassign
}

# Extract genome-wide TSS:
promoter_regions <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream=2000, downstream=2000)
head(promoter_regions)

# Generate score matrix & plot
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot.450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot.450k <- as.data.frame(annot.450k@listData)
annot.450k <- data.frame(
  chr = annot.450k$chr,
  start = as.numeric(annot.450k$pos),
  end = as.numeric(annot.450k$pos),
  cgID = annot.450k$Name
)
rownames(annot.450k) <- annot.450k$cgID

# Iteratively generate GRanges objects and make composite plots:
identical(gsub("\\:.*", "", colnames(betas.corePun)), gsub("\\:.*", "", colnames(betas.surgExc))) #check
identical(gsub("\\:.*", "", colnames(betas.corePun)), LETTERS[1:13])

par(mfrow=c(4, 4))
for(i in 1:13) {
  GR_corePun <- merge(annot.450k, betas.corePun[ , i, drop=F], by="row.names")
  colnames(GR_corePun)[6] <- "beta"
  GR_corePun <- makeGRangesFromDataFrame(GR_corePun, keep.extra.columns = TRUE)
  GR_corePun@elementMetadata@listData$Row.names <- NULL
  GR_corePun@elementMetadata@listData$cgID <- NULL
  GR_corePun@elementMetadata@listData$UCSCgene <- NULL
  
  GR_surgExc <- merge(annot.450k, betas.surgExc[ , i, drop=F], by="row.names")
  colnames(GR_surgExc)[6] <- "beta"
  GR_surgExc <- makeGRangesFromDataFrame(GR_surgExc, keep.extra.columns = TRUE)
  GR_surgExc@elementMetadata@listData$Row.names <- NULL
  GR_surgExc@elementMetadata@listData$cgID <- NULL
  GR_surgExc@elementMetadata@listData$UCSCgene <- NULL
  
  sml <- ScoreMatrixList(
    targets = list(Biopsy = GR_corePun, Surgical = GR_surgExc), 
    windows = promoter_regions, 
    bin.num = 75,
    strand.aware = FALSE,
    weight.col = "beta",
    is.noCovNA = TRUE, 
    cores = detectCores()-1
  )
  
  plotMeta(
    sml, 
    line.col = brewer.pal(n = 9, name = "Set1"), 
    xcoords = c(-2000, 2000),
    main = paste("Subject", LETTERS[i], sep=" "),
    ylim = c(0.15, 0.55),
    xlab = "Distance to TSS (bp)",
    ylab = "Beta"
  )
  #legend("bottomright", c("Biopsy", "Surgical"), lty=c(1, 1), lwd=c(2.5, 2.5), col=brewer.pal(n = 9, name = "Set1"))
}
