######################################################################################################
# For Manuscript: Beta-Value Distributions
# Script author: David (Youdinghuan) Chen <youinghuan.chen.gr@dartmouth.edu>
# Date: 12/08/2016
######################################################################################################

library(ggplot2); library(reshape); library(grid); library(gridExtra);  library(RColorBrewer)
library(matrixStats)
library(doParallel);registerDoParallel(detectCores() - 1)

load("/Users/DavidKevinChen/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/113016-stratified-beta-matrices.RData")

## Create a specimen key:
master$specimenKey <- paste(master$bcid, master$tissue_descrip, sep="-")

## Iteratively change the col.names of betas (`key`) to respective `specimenKey`:
## Match names:
betas <- betas[ , match(master$key, colnames(betas)) ] #reorder the columns
identical(colnames(betas), master$key) #checkpoint: TRUE required
colnames(betas) <- master$specimenKey

## Letter code the paired samples iteratively
## Each pair of the sample will get the same name
master$hmID <- NA
for(i in 1:(nrow(master)/2) ) {
  x <- unique(as.character(master$bcid))[i]
  master$hmID[grep(x, master$bcid)] <- c(LETTERS[i], LETTERS[i])
}

## Reassign column names iteratively using the master$hmID column:
for(i in 1:(nrow(master)/2) ) {
 index <- grep(colnames(betas.corePun)[i], master$key) #search
 colnames(betas.corePun)[i] <- paste(master$hmID[index], "biopsy", sep=":") #reassign
}

## Reassign column names iteratively using the master$hmID column:
for(i in 1:(nrow(master)/2) ) {
  index <- grep(colnames(betas.surgExc)[i], master$key) #search
  colnames(betas.surgExc)[i] <- paste(master$hmID[index], "surgical", sep=":") #reassign
}

identical(rownames(betas.corePun), rownames(betas.surgExc)) #check
betas_relabeled <- cbind(betas.corePun, betas.surgExc)
betas_relabeled <- betas_relabeled[ , order(colnames(betas_relabeled))]
colnames(betas_relabeled) #check

## Break columns into pairs of 2's and store each pair in a list object:
paired_list <- list()
odd <- function(x) x %% 2 != 0
for(i in seq(nrow(master))[odd(1:nrow(master))]) {
  name <- gsub("\\:.*", "", colnames(betas_relabeled)[i])
  paired_list[[name]] <- betas_relabeled[ , i:(i+1)]
}

## Store each "multiplot object" in a separate R list object:
plots <- list()  # new empty list
for(i in 1:13) {
  pi <- ggplot(
    melt(paired_list[[i]], varnames = c("CpG","Sample")), aes(value, fill=Sample, colour=Sample), expand=c(0,0)) +
    geom_density(alpha=0.2) + 
    scale_x_continuous(expand = c(0, 0)) +
    labs(title=paste("Subject", names(paired_list)[i], sep=" ")) + 
    guides(fill=guide_legend(keywidth = 1, keyheight = 1)) + 
    theme(legend.title=element_text(size=0.5), panel.background = element_blank(), axis.text.x=element_blank()
)
  plots[[i]] <- pi  # add each plot into plot list
}
source("~/repos/DCIS-Pathology-Paper-2016/multi_ggplot2.R") ## R code from R Cookbook online
multiplot(plotlist=plots, layout=matrix(c(1:13,NA,NA,NA), nrow=4, byrow=TRUE))

#-----------------------------(Optional) All on one page-----------------------------
ggplot(melt(betas_relabeled, varnames = c("CpG","Sample")), aes(value, fill=Sample, colour=Sample), expand=c(0,0)) +
  geom_density(alpha=0.2) + 
  scale_x_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(keywidth = 1, keyheight = 1)) + 
  labs(title="DCIS") +
  theme(legend.title=element_text(size=0.5))

corePun_resh <- melt(betas.corePun, varnames = c("CpG","Sample")) ##row, column
ggplot(corePun_resh, aes(value, fill=Sample, colour=Sample), expand=c(0,0)) +
  geom_density(alpha=0.2) + 
  scale_x_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(keywidth = 1, keyheight = 1)) + 
  labs(title="DCIS, biopsies") +
  theme(legend.title=element_text(size=0.5)) + 
  facet_wrap(~Sample, ncol=3) + theme(legend.position="none")

surgExc_resh <- melt(betas.surgExc, varnames = c("CpG","Sample")) ##row, column
ggplot(surgExc_resh, aes(value, fill=Sample, colour=Sample), expand=c(0,0)) +
  geom_density(alpha=0.2) + 
  scale_x_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(keywidth = 1, keyheight = 1)) + 
  labs(title="DCIS, sugrical excision") +
  theme(legend.title=element_text(size=0.5)) + 
  facet_wrap(~Sample, ncol=3) + theme(legend.position="none")
