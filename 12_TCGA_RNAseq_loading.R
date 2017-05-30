####################################################################################################
# Loading, cleaning, and visualizing TCGA stage I/Ia/Ib ER+ breast cancer RNAseq normalized counts
# Script author: David Chen <youdinghuan.chen.gr@dartmouth.edu>
# Start date: 03/15/2017
####################################################################################################

library(matrixStats)
library(pheatmap)

load("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/RNAseq&450kSampleAnnot.Rdata")

## Select samples and cbind into data matrix:
setwd('~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/Unpacked_gene_expr_quant/')
sum(RNAseq_meta_filtered$file_name %in% list.files()) #checkpoint: 51 expected

## Binding individual normalized_results files (txt) into a data.frame:
data.RNAseqV2 <- NULL #initialize
for(f in RNAseq_meta_filtered$file_name){
  test <- as.matrix(read.table(f, header=T, sep="\t", row.names=1, stringsAsFactors=F))
  colnames(test) <- as.character(f)
  data.RNAseqV2 <- cbind(data.RNAseqV2, test)
}
dim(data.RNAseqV2)
hist(log2(data.RNAseqV2))

RNAseq_meta_filtered$Sample_Name <- substring(RNAseq_meta_filtered$cases, 1, 12)
anyDuplicated(RNAseq_meta_filtered$Sample_Name) # 0 required to proceed
all(RNAseq_meta_filtered$file_name %in% colnames(data.RNAseqV2)) #TRUE required to proceed
data.RNAseqV2 <- data.RNAseqV2[ , match(RNAseq_meta_filtered$file_name, colnames(data.RNAseqV2))]
identical(RNAseq_meta_filtered$file_name, colnames(data.RNAseqV2)) #checkpoint: TRUE required
colnames(data.RNAseqV2) <- RNAseq_meta_filtered$Sample_Name

## Subset protocadherins and visualize expression:
## Retrieve genes
geneList <- read.csv("../../Paper Drafts/Supplemental Table 1.csv", skip=1, header=F, stringsAsFactors=F)
my_genes <- geneList$V1[grep("PCDH", geneList$V1[1:25])]
paste(my_genes, collapse=", ")
## [1] "PCDHGA2, PCDHGA4, PCDHGA1, PCDHGB1, PCDHGA3, PCDHGA5, PCDHGB2, PCDHGB3, PCDHGB4, PCDHGA6, PCDHGA7"

RNAseqV2_pcdh <- data.RNAseqV2[gsub("\\|.*", "", rownames(data.RNAseqV2)) %in% my_genes, ] 
sum(is.na(RNAseqV2_pcdh))

## Variance barplot
pcdh_vars <- rowVars(RNAseqV2_pcdh); names(pcdh_vars) <- rownames(RNAseqV2_pcdh)
barplot(pcdh_vars, las=2, names.arg=gsub("\\|.*", "", names(pcdh_vars)), ylab="Sample variance in gene expression", border=F)
pheatmap(
  RNAseqV2_pcdh,
  color = colorRampPalette(c("white", "tomato2"))(1024),
  labels_row = gsub("\\|.*", "", rownames(RNAseqV2_pcdh)),
  border_color = NA, 
  main = "TCGA RNAseqV2 normalized_counts (n=51 ER+ stage I/Ia/Ib)"
)

## Save RNAseq data
# save(list=c("data.RNAseqV2", "RNAseqV2_pcdh"), file="../20170315_TCGA_lev3_RNAseq.RData", compress=TRUE)
