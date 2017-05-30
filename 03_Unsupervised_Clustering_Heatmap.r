######################################################################################################
# For Manuscript: Heat Maps with Clustering by Manhattan Distance & Average Linkage Metric
# Script author: David (Youdinghuan) Chen <youinghuan.chen.gr@dartmouth.edu>
# Start Date: 12/08/2016
# Notes:
######################################################################################################

library(matrixStats)
library(pheatmap); library(RColorBrewer)
library(doParallel);registerDoParallel(detectCores() - 1)

load("/Users/DavidKevinChen/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/113016-stratified-beta-matrices.RData")

## Check dimensions & matched samples:
dim(master)
dim(betas)
table(master$bcid)
table(master$casepr)
all(master$key %in% colnames(betas)) #checkpoint: TRUE required

## Create a specimen key:
master$specimenKey <- paste(master$bcid, master$tissue_descrip, sep="-")

## Iteratively change the col.names of betas (`key`) to respective `specimenKey`:
## Match names:
betas <- betas[ , match(master$key, colnames(betas)) ] #reorder the columns
identical(colnames(betas), master$key) #TRUE expected
colnames(betas) <- master$specimenKey

## Letter code the paired samples iteratively. ## NOTE: This Operation will differ from preliminary results ##
## Each pair of the sample will get the same name
master$hmID <- NA
for(i in 1:(nrow(master)/2)) {
  x <- unique(as.character(master$bcid))[i]
  master$hmID[grep(x, master$bcid)] <- c(LETTERS[i], LETTERS[i])
}

# Update the values of certain varaibles:
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

## Make a copy of beta data.frame for plotting:
betas2 <- betas #initialize
identical(colnames(betas2), master$specimenKey) #TRUE required
for(i in 1:nrow(master)){
  colnames(betas2)[i] <- paste(master$hmID[i], master$tissue_descrip[i], sep=":") #complete de-identify
  rownames(master)[i] <- paste(master$hmID[i], master$tissue_descrip[i], sep=":") #update master accordingly
}
identical(rownames(master), colnames(betas2)) #check before plotting

## Annotation colors:
ann_cols <- list(
  names = rownames(master),
  Pattern = c(Cribriform="orange", Papillary="purple", Solid="cyan", Comedo="#1B9E77"),
  Grade = c(Intermediate="violet", High="royalblue3"), ##Low grade is not present after dropping samples
  Necrosis = c(Absent="blue", Present="red"),
  Sclerosis = c(Absent="blue", Present="red"),
  Inflammation = c(Absent="blue", Present="red"),
  Calcifications = c(Absent="blue", Present="red"),
  Specimen = c(Biopsy="black", Surgical="gray")
)

## Annotation matrix:
heat_annot <- data.frame(
  row.names = rownames(master),
  Pattern = as.factor(master$Pattern),
  Grade = as.factor(master$Grade),
  Necrosis = as.factor(master$Necrosis),
  Sclerosis = as.factor(master$Sclerosis),
  Inflammation = as.factor(master$Inflammation),
  Calcifications = as.factor(master$Calcs),
  Specimen = as.factor(master$tissue_descrip)
)

## Avoid rotating labels of heatmap by 270 degrees by running the line below:
## Avoid rotating labels of heatmap by 270 degrees by running the line below:
## Update the following arguments: `vjust = 1`, `rot = 0`
fixInNamespace("draw_colnames","pheatmap")
## When exiting and restarting R, the Namespace will reset.

## Hierarchical clustering with multiple CpGs:
for(n in c(1000, 2500, 5000, 1.5*10^4, 2*10^4)) {
  pheatmap(betas2[order(rowVars(betas2), decreasing=TRUE)[1:n], ], #most variable n*CpGs
           labels_row = rep("", nrow(betas2)),
           color = colorRampPalette(c("yellow", "blue"))(1024),
           labels_col = gsub("\\:.*", "",colnames(betas2)),
           clustering_distance_rows = "manhattan",
           clustering_distance_cols = "manhattan",
           clustering_method = "average",
           treeheight_row = 0, #NOTE: Do NOT set cluster_rows = FALSE, 
           annotation_col = heat_annot, 
           annotation_colors = ann_cols,
           border_color = NA
           #main=paste("Top", n, "most variable CpGs in matched DCIS specimens") #Omit when generating for paper
  )
}

save("betas2", file="~/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/121316_DCIS_betas2_for_RPMM.RData")

## Variance distribution as justification of using 1,000 CpGs
load("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/121316_DCIS_betas2_for_RPMM.RData")
Index <- 1:nrow(betas2)
plot(
  Index,
  sort(rowVars(betas2), decreasing=TRUE),
  col = ifelse(Index <= 1000, "red", "black"),
  bty="l", pch=16, cex=0.5, ylab="Variance in methylation beta",
  main="Distribution of sample variance (n=13 pairs)"
)

for(n in c(1000, 2500, 5000)){
  var_beta <- sort(rowVars(betas2), decreasing=TRUE)[n]
  print(paste("The sample variance is", var_beta, "when at the number of CpG of", n))
  abline(h=var_beta, col="darkgray", lty=2)
  text(2e05, var_beta+0.002, labels=paste(n, "most variable CpGs"))
}
