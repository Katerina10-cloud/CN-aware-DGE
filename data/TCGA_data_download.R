### Download data from GDC portal ###


#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/breast_cancer/data/")


# Get a list of projects #
gdcprojects = getGDCprojects()
getProjectSummary('TCGA-BRCA')

# build a query to retrieve data #
query_TCGA_cnv <- GDCquery(project = 'TCGA-BRCA',
        data.category = 'Copy Number Variation',
        #sample.type = "Primary Blood Derived Cancer - Peripheral Blood",
        data.type = "Gene Level Copy Number",
        workflow.type = 'ASCAT3')
        #barcode = luad_cnv_list)

output_query_TCGA <- getResults(query_TCGA_cnv)


# build a query to retrieve Tumor gene expression data #
query_TCGA_rna_tumor <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor",
                       access = 'open')
                       #barcode = luad_rna_tumor)

# build a query to retrieve Normal gene expression data #
query_TCGA_rna_normal <- GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Transcriptome Profiling',
                           experimental.strategy = 'RNA-Seq',
                           workflow.type = 'STAR - Counts',
                           data.type = "Gene Expression Quantification",
                           sample.type = "Solid Tissue Normal",
                           access = 'open')
                           #barcode = luad_rna_normal)
#download data
GDCdownload(query_TCGA_rna_normal)

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

brca_rna <- GDCprepare(query_TCGA_rna_normal, summarizedExperiment = TRUE)
brca_rna_normal <- assay(brca_rna, 'unstranded', rownames = TRUE)

gene_name <- as.data.frame(brca_rna@rowRanges@elementMetadata@listData[["gene_name"]]) 
colnames(gene_name)[1] <- "GeneID"
brca_rna_normal <- cbind(gene_name, brca_rna_normal)
brca_rna_normal <- brca_rna_normal[!duplicated(brca_rna_normal$GeneID), ] %>% remove_rownames %>% column_to_rownames(var="GeneID") %>%
  na.omit()

#substring columns
colnames <- substr(colnames(cnv_tumor), 1, 12)
colnames(cnv_tumor) <- colnames


#clinical data
clinical_brca <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", 
                          data.format = "bcr xml", barcode = barcode_tumor)
GDCdownload(clinical_brca)

clinical_brca <- GDCprepare_clinic(clinical_brca, clinical.info = "patient")
clinical_brca <- clinical_brca %>% select(1,6,22,107)
save(clinical_brca, file = "~/model_data/TCGA/breast_cancer/data/clinical_brca.Rdata")

gbm_rna_tumor <- gbm_rna_tumor[(rownames(gbm_rna_tumor) %in% rownames(gbml_cnv_tumor)),]
gbm_rna_tumor <- gbm_rna_tumor[,(colnames(gbm_rna_tumor) %in% colnames(gbml_cnv_tumor))]


#Filtering low counts genes
gbm_rna_tumor <- gbm_rna_tumor[which(rowSums(gbm_rna_tumor)>100),]
brca_rna_normal <- brca_rna_normal[(rownames(brca_rna_normal) %in% rownames(brca_rna_tum)),] #delete rows by name

# Reordering to match datasets #
colnames_idx <- match(colnames(cnv_tumor), colnames(brca_rna_tum))
brca_rna_tum <- brca_rna_tum[,colnames_idx]

brca_rna_tum <- brca_rna_tum[order(match(rownames(brca_rna_tum), rownames(cnv_tumor))),]

x <- colnames(brca_rna_tum)
names(brca_rna_tum) <- paste(x , "-01A")

save(metadata, file = 'metadata.Rdata')
