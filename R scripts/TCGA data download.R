### Download data from GDC portal ###


#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/")


# Get a list of projects #
gdcprojects = getGDCprojects()
getProjectSummary('TCGA-GBM')

# build a query to retrieve data #
query_TCGA_cnv <- GDCquery(project = 'TCGA-GBM',
        data.category = 'Copy Number Variation',
        #sample.type = "Primary Blood Derived Cancer - Peripheral Blood",
        data.type = "Gene Level Copy Number",
        workflow.type = 'ASCAT3')
        #barcode = luad_cnv_list)

output_query_TCGA <- getResults(query_TCGA_cnv)


# build a query to retrieve Tumor gene expression data #
query_TCGA_rna <- GDCquery(project = 'TCGA-GBM',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor",
                       access = 'open')
                       #barcode = luad_rna_tumor)

# build a query to retrieve Normal gene expression data #
query_TCGA_rna_normal <- GDCquery(project = 'TCGA-UCS',
                           data.category = 'Transcriptome Profiling',
                           experimental.strategy = 'RNA-Seq',
                           workflow.type = 'STAR - Counts',
                           data.type = "Gene Expression Quantification",
                           sample.type = "Solid Tissue Normal",
                           access = 'open')
                           #barcode = luad_rna_normal)
#download data
GDCdownload(query_TCGA_rna)

# clinical data download #
clinical_luad <- GDCquery(project = "TCGA-LAML", data.category = "Clinical", 
                          data.format = "bcr xml", barcode = barcode_tumor)
GDCdownload(clinical_luad)

clinical_luad <- GDCprepare_clinic(clinical_luad, clinical.info = "patient")
clinical_luad <- clinical_luad %>% select(1,6,7,69)
save(gbml_cnv_tumor, file = "gbm_cnv_tumor.Rdata")


### Prepare data ###
gbm_cnv <- GDCprepare(query_TCGA_cnv, summarizedExperiment = TRUE)
gbm_cnv_tumor <- assay(gbm_cnv, 'copy_number', rownames = TRUE) %>% as.data.frame()

gbm_rna <- GDCprepare(query_TCGA_rna, summarizedExperiment = TRUE)
gbm_rna_tumor <- assay(gbm_rna, 'unstranded', rownames = TRUE)

gene_name <- as.data.frame(gbm_rna@rowRanges@elementMetadata@listData[["gene_name"]]) 
colnames(gene_name)[1] <- "GeneID"
gbm_rna_tumor <- cbind(gene_name, gbm_rna_tumor)
gbm_rna_tumor <- gbm_rna_tumor[!duplicated(gbm_rna_tumor$GeneID), ] %>% remove_rownames %>% column_to_rownames(var="GeneID") %>%
  na.omit()

#substring columns
colnames(gbm_rna_tumor) <- substr(colnames(gbm_rna_tumor), 1, 12)


save(laml_cnv_tumor, file = '~/Documents/PhD_AI/TCGA/aml_cancer/laml_cnv_tumor.Rdata')

#clinical data
clinical_brca <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", 
                          data.format = "bcr xml", barcode = barcode_tumor)
GDCdownload(clinical_brca)

clinical_brca <- GDCprepare_clinic(clinical_brca, clinical.info = "patient")
clinical_brca <- clinical_brca %>% select(1,6,22,107)
save(clinical_brca, file = "~/model_data/TCGA/breast_cancer/data/clinical_brca.Rdata")

gbm_rna_tumor <- gbm_rna_tumor[(rownames(gbm_rna_tumor) %in% rownames(gbml_cnv_tumor)),]
gbm_rna_tumor <- gbm_rna_tumor[,(colnames(gbm_rna_tumor) %in% colnames(gbml_cnv_tumor))]

#rownames(countdata) <- seqdata[,1]
#substr("ThisIsAString", start=1, stop=5) #shorten sample names
#colnames(countdata) <- substr(colnames(countdata), 1, 7)
#colnames(countdata)==sampleinfo$SampleName
#all(colnames(countdata)==sampleinfo$SampleName)


#Filtering low counts genes
gbm_rna_tumor <- gbm_rna_tumor[which(rowSums(gbm_rna_tumor)>100),]
gbm_rna_tumor <- laml_cnv_tumor[(rownames(laml_cnv_tumor) %in% rownames(laml_rna_tumor)),] #delete rows by nam



