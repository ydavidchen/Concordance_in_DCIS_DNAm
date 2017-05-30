################################################################################################
# Visualizing TCGA 450k data & merge clinical covariates
# Script author: David Chen <youdinghuan.chen.gr@dartmouth.edu>
# Start date: 03/17/2017
################################################################################################

library(matrixStats)
library(pheatmap)

load("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/20170316_TCGA_stage1ERpos_betas.RData")
load("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/RNAseq&450kSampleAnnot.Rdata")

#------------------------------------------Prepare complete sample sheets for methylation data------------------------------------------
sample_sheet$Sample_Name <- gsub("-01[A-Z]", "", sample_sheet$Sample_Name)
all(sample_sheet$Sample_Name %in%  clin_stageI_ERpos$bcr_patient_barcode) #checkpoint: TRUE required

sample_sheet$Sentrix <- paste(sample_sheet$Sentrix_ID, sample_sheet$Sentrix_Position, sep="_")
identical(sample_sheet$Sentrix, colnames(betas_stage1ERpos)) # TRUE required to proceed; if false, run commented line below:
# betas_stage1ERpos <- betas_stage1ERpos[ , match(sample_sheet$Sentrix, colnames(betas_stage1ERpos))]

colnames(betas_stage1ERpos) <- sample_sheet$Sample_Name #update sample names from chip barcode to TCGA barcode

master_clincov <- clin_stageI_ERpos #initialize
master_clincov$Sample_Name <- master_clincov$bcr_patient_barcode
master_clincov <- merge(master_clincov, sample_sheet, by="Sample_Name")
nrow(master_clincov) #checkpoint: 51 expected

## Check & drop the duplicated columns in methylation annotation:
col_indices <- which(colnames(master_clincov) %in% colnames(RNAseq_meta_filtered))
col_indices
master_clincov <- master_clincov[ , -col_indices]

RNAseq_meta_filtered$Sample_Name <- substring(RNAseq_meta_filtered$cases, 1, 12)
master_clincov <- merge(master_clincov, RNAseq_meta_filtered, by="Sample_Name")
nrow(master_clincov) #checkpoint: 51 expected also

# write.csv(master_clincov, "~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/20170317_clin_covar_with_Sentrix.csv", row.names=F, quote=F)

#------------------------------------------Subset protocadherins genes and visualize methylation------------------------------------------
# Protocad CpGs of interest:
protocad_CpGs <- strsplit("cg00978427,cg27553119,cg23445461,cg26890354,cg15949044,cg02119152,cg13972793,cg02022808", split=",", fixed=T)[[1]]
protocad_CpGs

protocad_methyl_TCGA <- betas_stage1ERpos[rownames(betas_stage1ERpos) %in% protocad_CpGs, ]
pheatmap(
  protocad_methyl_TCGA, 
  border_color=NA, 
  color = colorRampPalette(c("blue", "yellow"))(1024),
  main="TCGA minfi-processed 450k methylation beta (n=51 ER+ stage I/Ia/Ib)"
)

## Variance barplot in order of genomic position:
barplot(rowVars(protocad_methyl_TCGA), names.arg=rownames(protocad_methyl_TCGA), las=2, border=F, ylim=c(0, 0.08), ylab="Sample variance in beta", main="TCGA ER+ stage I/Ia/Ib")


#------------------------------------------Scatter plots (11 mRNA * 8 CpGs = 96 plots total)------------------------------------------
## Start with the 3 sites covered by all transcripts
load("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/20170315_TCGA_lev3_RNAseq.RData")
all(colnames(RNAseqV2_pcdh) %in% colnames(betas_stage1ERpos)) #checkpoint: TRUE required
## [1] TRUE
identical(colnames(RNAseqV2_pcdh), colnames(betas_stage1ERpos)) #if FALSE, use `match` and re-run the line
## [1] TRUE

## Consider saving the two objects for plotting as a separate RData for sharing/easy loading
# save(list=c("master_clincov", "protocad_methyl_TCGA", "RNAseqV2_pcdh"), file="~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/20170320_meth_mRNA_dfs.RData")

## Check data structure:
View(RNAseqV2_pcdh)
View(protocad_methyl_TCGA)

dim(RNAseqV2_pcdh)
dim(protocad_methyl_TCGA)

## Visualize CpG::expression one-by-one:
# pdf("~/Downloads/20170320_one_by_one_comparisons.pdf", paper="a4r")
for(cg in 1:nrow(protocad_methyl_TCGA)){
  par(mfrow=c(3, 4)) #reset graphics per iteration of CpG
  for(t in 1:nrow(RNAseqV2_pcdh)){
    plot(
      RNAseqV2_pcdh[t, ], protocad_methyl_TCGA[cg, ],
      xlab = "RNAseqV2 counts", ylab = "methylation beta", 
      # xlim=c(0,1500), 
      ylim=c(0, 1),
      main = paste(rownames(protocad_methyl_TCGA)[cg], "vs.", rownames(RNAseqV2_pcdh)[t]),
      pch=16, bty="l", cex=0.5, cex.main=0.6, cex.lab=0.6, cex.axis=0.6
    )
  }
}
# dev.off()

#------------------------------------------Fit regression line only on the CpGs shared across all transcripts------------------------------------------
## Focus on the CpGs shared across all gene-body transcripts
load("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/20170320_meth_mRNA_dfs.RData")

library(gridExtra)
library(gridBase)
library(grid)
grobTheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex=0.5), bg_params = list(fill=NA, col=1)),
  colhead = list(fg_params=list(cex=0.5), bg_params = list(fill=NA, col=1)),
  rowhead = list(fg_params=list(cex=0.5), bg_params = list(fill=NA, col=1))
)

## Directly export from the RStudio plotting canvas:
for(cg in which(rownames(protocad_methyl_TCGA) %in% c("cg26890354", "cg15949044", "cg23445461")) ){
  par(mfrow=c(3, 4)) #reset graphics per iteration of CpG
  vps <- baseViewports()
  stat_table <- NULL #initialize table
  
  for(t in 1:nrow(RNAseqV2_pcdh)){
    ## Statistical test:
    my_fit <- lm(protocad_methyl_TCGA[cg, ] ~ RNAseqV2_pcdh[t, ])
    r.sq <- summary(my_fit)$r.squared
    p.val <- summary(my_fit)$coefficients[8]
    stat_table <- rbind(stat_table, c(
      Panel=letters[t], 
      R.squared=formatC(round(r.sq,3), format="f", digits=3), 
      P.value=formatC(round(p.val,3), format="f", digits=3)
    ))
    ## Plot w/ labels, plus a table
    plot(
      RNAseqV2_pcdh[t, ], protocad_methyl_TCGA[cg, ],
      xlab = "RNAseqV2 counts", ylab = "methylation beta", 
      # xlim=c(0,1500), 
      ylim=c(0, 1),
      main = paste(letters[t], ")", rownames(protocad_methyl_TCGA)[cg], "vs.", rownames(RNAseqV2_pcdh)[t]),
      pch=16, bty="l", cex=0.5, cex.main=0.75, cex.lab=0.75, cex.axis=0.5
    )
    abline(my_fit$coefficients)
  }
  
  ## Add summary table:
  pushViewport(vps$inner, vps$figure, vps$plot)
  grid.draw(tableGrob(stat_table, theme=grobTheme)) #with R.square AND P.value columns
}
