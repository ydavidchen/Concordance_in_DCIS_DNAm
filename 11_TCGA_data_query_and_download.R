####################################################################################################
# Downloading TCGA methylation and gene expression data from Genomic Data Commons
# Script author: David Chen <youdinghuan.chen.gr@dartmouth.edu>
# Start date: 03/15/2017
####################################################################################################

library(TCGAbiolinks)
library(DT)
library(DESeq2)

## Retrieve clinical data.frame with *complete* barcode and race:
## Hormone receptor status & distal metastasis are missing from the clinical data
clinical <- GDCquery_clinic(project="TCGA-BRCA", type="clinical")
class(clinical) <- "data.frame"
table(clinical$tumor_stage)
epi_condition <- (clinical$tumor_stage=="stage i" | clinical$tumor_stage=="stage ia" | clinical$tumor_stage == "stage ib") &
  clinical$race=="white" &
  clinical$gender=="female"
sum(epi_condition)
clin_sub <- clinical[epi_condition, ]

## Retrieve subtype data, which is NOT in the `clinical` data.frame:
subtype_df <- TCGAquery_subtype(tumor="BRCA")
table(subtype_df$AJCC.Stage) #stages
table(subtype_df$Tumor)
table(subtype_df$Metastasis)
table(subtype_df$Gender)

extra_condition <- subtype_df$ER.Status == "Positive" & 
  subtype_df$Metastasis == "M0" #no metastasis 
subtype_sub <- subtype_df[extra_condition, ]
nrow(subtype_sub)
subtype_sub <- subtype_sub[ , c("patient", "ER.Status", "Metastasis")]

# ER_pos_no_metas <- subtype_sub$patient
# sum(clin_sub$submitter_id %in% ER_pos_no_metas)
# clin_stageI_ERpos <- clin_sub[clin_sub$submitter_id %in% ER_pos_no_metas, ]

clin_stageI_ERpos <- merge(clin_sub, subtype_sub, by=1) #merge by common barcodes
identical(clin_stageI_ERpos$submitter_id, clin_stageI_ERpos$bcr_patient_barcode) #checkpoint: TRUE required
## [1] TRUE
clin_stageI_ERpos$treatments <- NULL #drop in order to export as CSV
# write.csv(clin_stageI_ERpos, file="~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/TCGA_covar_stageI_ER.csv")
my_barcodes <- clin_stageI_ERpos$submitter_id
anyNA(my_barcodes)
length(my_barcodes)

## Retrieve level 3 RNA-seq data aligned against hg19:
setwd("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/")
query.exp.hg19 <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE,
  file.type  = "normalized_results",
  platform = "Illumina HiSeq", 
  barcode = my_barcodes
)
GDCdownload(query.exp.hg19)
RNAseq_data <- GDCprepare(query.exp.hg19)
RNAseq_meta <- getResults(query.exp.hg19)
datatable(
  as.data.frame(colData(RNAseq_data)), 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 25), 
  rownames = FALSE
) #initial exploration in Viewer
# save(list=c("clin_stageI_ERpos", "my_barcodes", "RNAseq_data"), file="GDCprepared_level3_brca_RNAseq.RData", compress=TRUE)

## Retrieve level 1 DNA methylation data:
setwd("~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDCdata_HM450")
query.HM450 <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Raw microarray data",
  data.type = "Raw intensities", 
  experimental.strategy = "Methylation array", 
  legacy = TRUE,
  file.type = ".idat",
  platform = "Illumina Human Methylation 450", 
  barcode = my_barcodes
)
GDCdownload(query.HM450)
HM450_meta <- getResults(query.HM450) #clin annotation for downloaded data

## Confirm only 2 types of tissue specimens exist: tumor vs. normal
table(droplevels(HM450_meta$tissue.definition))
table(droplevels(RNAseq_meta$tissue.definition))

mean(substr(HM450_meta$cases, 1, 16) %in% substr(RNAseq_meta$cases, 1, 16)) 
mean(substr(RNAseq_meta$cases, 1, 16) %in% substr(HM450_meta$cases, 1, 16)) 

## Filtering: 
## Drop normals
## And keep rows with both methylation and RNAseq data files
file_cond <- substr(RNAseq_meta$cases,1,16) %in% substr(HM450_meta$cases,1,16) &
  RNAseq_meta$tissue.definition == "Primary solid Tumor"
sum(file_cond)

RNAseq_meta_filtered <- RNAseq_meta[file_cond, ] #51 total
HM450_meta_filtered <- HM450_meta[substr(HM450_meta$cases,1,16) %in% substr(RNAseq_meta_filtered$cases,1,16), ] #51 subjects * 2 idats per subject = 102

# save(list = c("clin_stageI_ERpos", "HM450_meta", "RNAseq_meta", "HM450_meta_filtered", "RNAseq_meta_filtered"),
#      file = "~/Dropbox (Christensen Lab)/Christensen Lab - Fall 2016/F2016 DCIS pathology paper/GDC_ERpos_datasets/RNAseq&450kSampleAnnot.Rdata", compress=T)
