####################################################################################################################
# Permutation-based comparison of within- vs. between-subject variation with iterations
# Script author: David (Youdinghuan) Chen <youinghuan.chen.gr@dartmouth.edu>
# Date: February 14, 2017
####################################################################################################################

library(doParallel); registerDoParallel(detectCores() - 1)
library(matrixStats)
library(scales)
library(ggplot2)
library(grid)
library(gridBase)
library(gridExtra)

load("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/121316_DCIS_betas2_for_RPMM.RData")

#--------------------------------------Distribution 1: Delta betas for matched specimens--------------------------------------
## Subset by specimen type:
biopsies <- betas2[ , grep("Biopsy", colnames(betas2))]
colnames(biopsies) <- gsub("\\:.*", "", colnames(biopsies))

surgicals <- betas2[ , grep("Surgical", colnames(betas2))]
colnames(surgicals) <- gsub("\\:.*", "", colnames(surgicals))

## Check labels:
identical(colnames(biopsies), colnames(surgicals)) #TRUE required
identical(rownames(biopsies), rownames(surgicals)) #TRUE required

## Compute Delta_beta matrix:
Delta_betas <- biopsies - surgicals

#--------------------------------------Distribution 2: Permuted distribution with n iterations--------------------------------------
identical(dimnames(biopsies), dimnames(surgicals)) #checkpoint: TRUE required
identical(colnames(biopsies), LETTERS[1:13]) #checkpoint: TRUE required

## Iterate through all:
## Do this n times by a giant for-loop
n <- 4

## Iteratively initialize n empty matrices in a single list
perm_distribution_list <- list() #empty list

## The "giant" for-loop:
for(iter in 1:n){
  ## Initialize empty matrices within the list:
  perm_distribution_list[[iter]] <- matrix(NA, nrow=nrow(betas2), ncol=13, dimnames=list(rownames(betas2), LETTERS[1:13])) #initialize
  
  ## Nested-for loop for each matrix (distribution) in the list object:
  for(j in 1:13){
    surg_samp <- surgicals[ , -j] #matrix without self to sample from
    for(i in 1:nrow(betas2)){
      perm_distribution_list[[iter]][i, j] <-  sample(surg_samp[i, ], 1) ##update
    }
  }

}

## Random checkpoint:
perm_surgical_i <- perm_distribution_list[[sample(1:n, 1)]]
row <- sample(1:nrow(perm_surgical_i), 1, replace=T)
col <- sample(1:13, 1, replace=T)
c(row, col)
perm_surgical_i[row, col] %in% surgicals[row, c(1:13)[-col]] #TRUE required

## Complete calculating the second distribution:
for(iter in 1:n) perm_distribution_list[[iter]] <- biopsies - perm_distribution_list[[iter]]

# save(list=c("Delta_betas", "perm_distribution_list"), file="20170214_Permuted_Distributions_with_iterations.RData", compress=TRUE)

#--------------------------------------Compute Summary Measures--------------------------------------
load("/Users/DavidKevinChen/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/20170214_Permutation with Iterations/20170214_Permuted_Distributions_with_iterations.RData")
## Load Dr. Johnson's progression-related CpGs (already exported as CSV from XLS)
progressTab <- read.csv("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS Pathology Paper/13148_2015_94_MOESM4_ESM.csv", skip = 1, row.names = 1)
progressCpGs <- rownames(progressTab)

## Compute summary measure:
median_perm_distr <- apply(simplify2array(perm_distribution_list), 1:2, median)
mean_perm_distr <- apply(simplify2array(perm_distribution_list), 1:2, mean)

#--------------------------------------Data Visualization (exploratory)--------------------------------------
## Initialize summary table:
pval_table <- NULL

## Median in all CpGs
test_res_by_subj <- setNames(rep(NA, 13), LETTERS[1:13]) #empty P-value vector
op <- par(oma=c(2,0,3,0), mfrow=c(3,5)) #graphic param
for(j in 1:13) {
  x <- Delta_betas[ , j]
  y <- median_perm_distr[ , j]
  
  plot(density(x), main=paste("Subject", LETTERS[j]), col=alpha('red',0.5), xlab=expression("Median"~Delta~"beta"), xlim=c(-1,1), ylim=c(0,15), lwd=2, bty="n")
  lines(density(y), col=alpha("skyblue",0.5), lwd=2)
  
  ## Perform test and store P-values:
  ks <- ks.test(x, y, alternative = "two.sided")
  test_res_by_subj[j] <- ks$p
}
par(op)
op <- par(usr=c(3,0,3,0), xpd=NA) #reset coordinate
legend(0.75, 3, fill=c(alpha('red',0.5), alpha('skyblue', 0.5)), cex=0.75,
       legend=c("biopsy - matched surgical", "biopsy - permuted surgical"), bty = "n")
test_res_by_subj
pval_table <- cbind(pval_table, median_all=test_res_by_subj) #first update


## Median across 641 CpGs
test_res_by_subj <- setNames(rep(NA, 13), LETTERS[1:13]) #empty P-value vector
op <- par(oma=c(0,0,3,0), mfrow=c(3,5)) #graphic param
for(j in 1:13) {
  x <- Delta_betas[rownames(Delta_betas) %in% progressCpGs, j]
  y <- median_perm_distr[rownames(median_perm_distr) %in% progressCpGs, j]

  plot(density(x), main=paste("Subject", LETTERS[j]), col=alpha('red',0.5), xlab=expression("Median"~Delta~"beta"), xlim=c(-1,1), ylim=c(0,4), lwd=2)
  lines(density(y), col=alpha("skyblue",0.5), lwd=2)
  
  ## Perform test and store P-values:
  ks <- ks.test(x, y, alternative = "two.sided")
  test_res_by_subj[j] <- ks$p
  
  ## Add ks P-value to plot
  if(ks$p.value == 0){
    text(0.5, 0.5, "P < 2.20e-16")} else{
      text(0.5, 0.5, paste("P =", signif(ks$p.value, 3)))
    }
}
par(op)
op <- par(usr=c(3,0,3,0), xpd=NA) #reset coordinate
legend(0.75, 3, fill=c(alpha('red',0.5), alpha('skyblue', 0.5)), cex=0.75,
       legend=c("biopsy - matched surgical", "biopsy - permuted surgical"), bty = "n")
test_res_by_subj
pval_table <- cbind(pval_table, median_641=test_res_by_subj) #second update


## Mean across all CpGs
test_res_by_subj <- setNames(rep(NA, 13), LETTERS[1:13]) #empty P-value vector
op <- par(oma=c(0,0,3,0), mfrow=c(3,5)) #graphic param
for(j in 1:13) {
  x <- Delta_betas[ , j]
  y <- mean_perm_distr[ , j]
  
  plot(density(x), main=paste("Subject", LETTERS[j]), col=alpha('red',0.5), xlab=expression("Mean"~Delta~"beta"), xlim=c(-1,1), ylim=c(0,15), lwd=2)
  lines(density(y), col=alpha("skyblue",0.5), lwd=2)
  
  ## Perform test and store results
  ks <- ks.test(x, y, alternative = "two.sided")
  test_res_by_subj[j] <- ks$p
  
  ## Add ks P-value to plot
  if(ks$p.value == 0){
    text(0.5, 0.5, "P < 2.20e-16")} else{
      text(0.5, 0.5, paste("P =", signif(ks$p.value, 3)))
    }
}
par(op)
op <- par(usr=c(3,0,3,0), xpd=NA) #reset coordinate
legend(0.75, 3, fill=c(alpha('red',0.5), alpha('skyblue', 0.5)), cex=0.75,
       legend=c("biopsy - matched surgical", "biopsy - permuted surgical"), bty = "n")
test_res_by_subj
pval_table <- cbind(pval_table, mean_all=test_res_by_subj) #third update


## Mean across 641 CpGs
test_res_by_subj <- setNames(rep(NA, 13), LETTERS[1:13]) #empty P-value vector
op <- par(oma=c(0,0,3,0), mfrow=c(3,5)) #graphic param
for(j in 1:13) {
  x <- Delta_betas[rownames(Delta_betas) %in% progressCpGs, j]
  y <- mean_perm_distr[rownames(mean_perm_distr) %in% progressCpGs, j]
  
  plot(density(x), main=paste("Subject", LETTERS[j]), col=alpha('red',0.5), xlab=expression("Mean"~Delta~"beta"), xlim=c(-1,1), ylim=c(0,4), lwd=2)
  lines(density(y), col=alpha("skyblue",0.5), lwd=2)
  
  ## Perform test and store results
  ks <- ks.test(x, y, alternative = "two.sided")
  test_res_by_subj[j] <- ks$p
  
  ## Add ks P-value to plot
  if(ks$p.value == 0){
    text(0.5, 0.5, "P < 2.20e-16")} else{
      text(0.5, 0.5, paste("P =", signif(ks$p.value, 3)))
    }
}
par(op)
op <- par(usr=c(3,0,3,0), xpd=NA) #reset coordinate
legend(0.75, 3, fill=c(alpha('red',0.5), alpha('skyblue', 0.5)), cex=0.75,
       legend=c("biopsy - matched surgical", "biopsy - permuted surgical"), bty = "n")
test_res_by_subj

#--------------------------------------Data Visualization (final): 641 CpGs--------------------------------------
## Clean up the P-value table:
pval_table <- cbind(pval_table, mean641=test_res_by_subj) #fourth update
rownames(pval_table) <- paste("Subject", rownames(pval_table))
pval_table[pval_table == 0] <- 2.20e-16
pval_table
colnames(pval_table) <- c("Median P-value in all CpGs", "Median P-value in 641 CpGs", "Mean P-value in all CpGs", "Mean P-value in 641 CpGs")
pval_table <- signif(pval_table, digits = 3)
pval_table <- pval_table[ ,-c(1,3)]
pval_table

## Keep the following code! **DO NOT DELETE**
# my_ggplots <- list()
# for(j in 1:13) {
#   temp_df <- data.frame(
#     db = Delta_betas[rownames(Delta_betas) %in% progressCpGs, j],
#     pdb = median_perm_distr[rownames(median_perm_distr) %in% progressCpGs, j]
#   )
#   
#   my_ggplots[[j]] <- ggplot(data = temp_df) + 
#     geom_density(aes(x=db), colour="red", fill="red", alpha=0.25) + 
#     geom_density(aes(x=pdb), colour="skyblue", fill="skyblue", alpha=0.25) +
#     xlim(c(-1,1)) +
#     ylim(c(0,4)) +
#     xlab(label = expression("Median"~Delta~"beta")) + 
#     ggtitle(paste("Subject", LETTERS[j])) +
#     theme_classic()
# }
# #my_ggplots[[length(my_ggplots)+1]] <- tableGrob(pval_table)
# 
# source("~/repos/DCIS-Pathology-Paper-2016/multi_ggplot2.R")
# multiplot(plotlist=my_ggplots, layout=matrix(c(1:13,"tab", NA, NA), ncol=5, byrow=T))

grobTheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex=0.3), bg_params = list(fill=NA, col=1)),
  colhead = list(fg_params=list(cex=0.3), bg_params = list(fill=NA, col=1)),
  rowhead = list(fg_params=list(cex=0.3), bg_params = list(fill=NA, col=1))
)

## First add the median distributions for 13 subjects:
op <- par(oma=c(0,0,2.3,0), mfrow=c(3,5))
for(j in 1:13) {
  x <- Delta_betas[rownames(Delta_betas) %in% progressCpGs, j]
  y <- median_perm_distr[rownames(median_perm_distr) %in% progressCpGs, j]
  
  plot(density(x), main=paste("Subject", LETTERS[j]), col=alpha('red',0.5), xlab=expression("Median"~Delta~"beta"), xlim=c(-1,1), ylim=c(0,4), lwd=2, bty="l")
  lines(density(y), col=alpha("skyblue",0.5), lwd=2)
  
  # ks <- ks.test(x, y, alternative = "two.sided")
  # if(ks$p.value == 0){
  #   text(0.5, 0.5, "P < 2.20E-16")} else{
  #     text(0.5, 0.5, paste("P =", signif(ks$p.value, 3)))
  #   }
}
# write.csv(pval_table, file="~/Downloads/20170305_pval_table_641CpGs.csv", quote=F)
par(op)
op <- par(usr=c(2,0,3,0), xpd=NA) #reset coordinate
legend(0.75, 3, fill=c(alpha('red',0.5), alpha('skyblue', 0.5)), cex=0.7,
       legend=c("biopsy - matched surgical", "biopsy - permuted surgical"), bty = "n")

## Add in the tableGrob:
op <- par(oma=c(0,0,2.3,0), mfrow=c(3,5))
#frame()
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
grid.draw(tableGrob(pval_table, theme=grobTheme))

# ## Add Common legend at top **DO NOT DELETE**
# reset <- function() {
#   par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
#   plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
# }
# reset()
# legend("top", fill=c(alpha('red',0.5), alpha('skyblue', 0.5)), cex=0.7, horiz=TRUE,
#        legend=c(expression(Delta~"beta"), expression("Permuted"~Delta~"beta")), bty = "n")

