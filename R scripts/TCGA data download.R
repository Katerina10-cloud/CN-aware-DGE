#loading libraries
#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/")


# Get a list of projects #
gdcprojects = getGDCprojects()
getProjectSummary('TCGA-OV')

# build a query to retrieve data #
query_TCGA_cnv <- GDCquery(project = 'TCGA-LAML',
        data.category = 'Copy Number Variation',
        #sample.type = "Primary Blood Derived Cancer - Peripheral Blood",
        data.type = "Gene Level Copy Number",
        workflow.type = 'ASCAT3')
        #barcode = luad_cnv_list)

output_query_TCGA <- getResults(query_TCGA_cnv)


# build a query to retrieve Tumor gene expression data #
query_TCGA_rna <- GDCquery(project = 'TCGA-LAML',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       #sample.type = "Primary Tumor",
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
save(clinical_luad, file = "~/model_data/TCGA/lung_cancer/LUAD/data/clinical_luad.Rdata")


### Prepare data ###
laml_cnv <- GDCprepare(query_TCGA_cnv, summarizedExperiment = TRUE)
laml_cnv_tumor <- assay(laml_cnv, 'copy_number', rownames = TRUE)

laml_rna <- GDCprepare(query_TCGA_rna, summarizedExperiment = TRUE)
laml_rna_tumor <- assay(laml_rna, 'unstranded', rownames = TRUE)

gene_name <- as.data.frame(laml_rna@rowRanges@elementMetadata@listData[["gene_name"]]) 
colnames(gene_name)[1] <- "GeneID"
laml_rna_tumor <- as.data.frame(laml_rna_tumor)
laml_rna_tumor <- cbind(gene_name, laml_rna_tumor)
laml_rna_tumor <- laml_rna_tumor[!duplicated(laml_rna_tumor$GeneID), ] %>% remove_rownames %>% column_to_rownames(var="GeneID")
laml_rna_tumor <- na.omit(laml_rna_tumor)

#substring columns
colnames(laml_rna_tumor) <- substr(colnames(laml_rna_tumor), 1, 12)


save(laml_cnv_tumor, file = '~/Documents/PhD_AI/TCGA/aml_cancer/laml_cnv_tumor.Rdata')

#clinical data
clinical_brca <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", 
                          data.format = "bcr xml", barcode = barcode_tumor)
GDCdownload(clinical_brca)

clinical_brca <- GDCprepare_clinic(clinical_brca, clinical.info = "patient")
clinical_brca <- clinical_brca %>% select(1,6,22,107)
save(clinical_brca, file = "~/model_data/TCGA/breast_cancer/data/clinical_brca.Rdata")

laml_rna_tumor <- laml_rna_tumor[(rownames(laml_rna_tumor) %in% rownames(res_allGenes)),]
laml_rna_tumor <- laml_rna_tumor[,(colnames(laml_rna_tumor) %in% colnames(laml_cnv_tumor))]

#rownames(countdata) <- seqdata[,1]
#substr("ThisIsAString", start=1, stop=5) #shorten sample names
#colnames(countdata) <- substr(colnames(countdata), 1, 7)
#colnames(countdata)==sampleinfo$SampleName
#all(colnames(countdata)==sampleinfo$SampleName)


#Filtering low counts genes
laml_rna_tumor <- laml_rna_tumor[which(rowSums(laml_rna_tumor)>100),]
laml_cnv_tumor <- laml_cnv_tumor[(rownames(laml_cnv_tumor) %in% rownames(laml_rna_tumor)),] #delete rows by nam



