######################################################################################################
# Summary of Delta betas by specimen types
# Script author: David (Youdinghuan) Chen <youinghuan.chen.gr@dartmouth.edu>
# Start Date: 01/19/2017
# Notes:
######################################################################################################

library(biomaRt)
library(matrixStats)
library(splitstackshape)
library(KEGGprofile)
library(doParallel); registerDoParallel(detectCores() - 1)

setwd("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/")

## Load data:
Delta_betas_summary <- read.csv("20170119_Delta_betas_summary_stats.csv", header=TRUE, row.names=1)
sum(Delta_betas_summary$Direction == "Up")
sum(Delta_betas_summary$Direction == "Down")

## Filter by Direction & genomic context:
# bool_mask <- ( (Delta_betas_summary$Relation_to_Island == "Island" & Delta_betas_summary$Direction == "Up") |
#   (grepl("Body", Delta_betas_summary$UCSC_RefGene_Group) & Delta_betas_summary$Direction == "Down") )
# Del_b_sele <- Delta_betas_summary[bool_mask, ]

Del_b_Up <- Delta_betas_summary[Delta_betas_summary$Direction == "Up", ]
UniqueGeneMatrix_Up <- Del_b_Up[ , "UCSC_RefGene_Name", drop=FALSE]
UniqueGeneMatrix_Up <- cSplit(UniqueGeneMatrix_Up, splitCols ="UCSC_RefGene_Name", drop=TRUE, sep=";") #Text-to-Column equivalent
rownames(UniqueGeneMatrix_Up) <- rownames(Del_b_Up)
## For each row, set duplicated elements (genes) to NA:
for(i in 1:nrow(UniqueGeneMatrix_Up)) {
  ## Find column indices for duplicated genes by row:
  col_indices <- which(duplicated(as.character(unlist(UniqueGeneMatrix_Up[i, ]))))
  ## Set any duplicated genes to NA:
  UniqueGeneMatrix_Up[i, col_indices] <- NA
}
## IMPORTANT: Iterative checkpoint:
iter_checkpt <- c()
for(row in 1:nrow(UniqueGeneMatrix_Up)) iter_checkpt <- c(iter_checkpt, anyDuplicated(na.omit(unlist(UniqueGeneMatrix_Up[row, ]))))
all(iter_checkpt == 0) # TRUE required
upGenes <- unique(unlist(UniqueGeneMatrix_Up))
upGenes <- upGenes[!is.na(upGenes)]

Del_b_Down <- Delta_betas_summary[Delta_betas_summary$Direction == "Down", ]
UniqueGeneMatrix_Down <- Del_b_Down[ , "UCSC_RefGene_Name", drop=FALSE]
UniqueGeneMatrix_Down <- cSplit(UniqueGeneMatrix_Down, splitCols ="UCSC_RefGene_Name", drop=TRUE, sep=";") #Text-to-Column equivalent
rownames(UniqueGeneMatrix_Down) <- rownames(Del_b_Down)
## For each row, set duplicated elements (genes) to NA:
for(i in 1:nrow(UniqueGeneMatrix_Down)) {
  ## Find column indices for duplicated genes by row:
  col_indices <- which(duplicated(as.character(unlist(UniqueGeneMatrix_Down[i, ]))))
  ## Set any duplicated genes to NA:
  UniqueGeneMatrix_Down[i, col_indices] <- NA
}
## IMPORTANT: Iterative checkpoint:
iter_checkpt <- c()
for(row in 1:nrow(UniqueGeneMatrix_Down)) iter_checkpt <- c(iter_checkpt, anyDuplicated(na.omit(unlist(UniqueGeneMatrix_Down[row, ]))))
all(iter_checkpt == 0) # TRUE required
downGenes <- unique(unlist(UniqueGeneMatrix_Down))
downGenes <- downGenes[!is.na(downGenes)]

#--------------------------------------Pathway Analysis by KEGG--------------------------------------
## Replace `method` argument in pvalueAdj <- p.adjust(pvalue, method = "BH")` with `fdr` on line 33:
fixInNamespace("find_enriched_pathway", "KEGGprofile")

## Convert hgnc.symbol to entrez gene IDs (Reference: `KEGGprofile::convertId`): 
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
newIdTable <- getBM(values=downGenes, attributes=c("hgnc_symbol", "entrezgene"), filters="hgnc_symbol", mart=ensembl)
newIdTable <- newIdTable[which(newIdTable[ , 1] != "" & newIdTable[ , 2] != ""), ]
myEntrezIDs <- as.character(newIdTable$entrezgene)

## Perform KEGG enrichment analysis:
myKEGG <- find_enriched_pathway(myEntrezIDs, species="hsa", returned_adjpvalue=0.05, download_latest=TRUE)
View(myKEGG[[1]])
# download_KEGGfile(pathway_id=c("04010", "04020"), species='hsa', target_dir="~/Downloads/") #optional download to view pathway

## Visualization of significantly enriched KEGG pathways:
myKEGG <- myKEGG[[1]]
myKEGG <- myKEGG[order(myKEGG$pvalueAdj, decreasing=FALSE), ]
myKEGG$`-log10(FDR)` <- -log10(myKEGG$pvalueAdj)

## Barplot:
barplot(
  myKEGG$`-log10(FDR)`, 
  names.arg = rownames(myKEGG), 
  #xaxt = "n",
  las=2, ylim=c(0,5), border=NA, col="darkblue", 
  ylab = expression("-log"[10]~"(FDR)"), 
  main = expression("Enriched KEGG pathways with" ~ Delta ~ "beta" <= -0.15 ~ "and FDR" <= 0.05)
)
abline(h=-log10(0.05), lty=4, lwd=2, col="red")
text(x=34, y=-log10(0.05)+0.17, labels="FDR = 0.05", col="red")
abline(h=-log10(0.01), lty=4, lwd=2, col="red")
text(x=34, y=-log10(0.01)+0.1, labels="FDR = 0.01", col="red")
legend("topright", bty="n", legend=paste(rownames(myKEGG), myKEGG$Pathway_Name, sep=": "), cex=0.7)

## Export data:
# write.csv(myKEGG, "~/repos/DCIS-Pathology-Paper-2016/Analysis-Results-Figures/Down_group_KEGG_analysis.csv", row.names=TRUE, quote=FALSE)
