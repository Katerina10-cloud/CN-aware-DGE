BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

#Get a list of projects
gdcprojects = getGDCprojects()
getProjectSummary('TCGA-LUAD')

#Building a query
query_TCGA <- GDCquery(project = 'TCGA-LUAD',
        data.category = 'Transcriptome Profiling')

output_query_TCGA <- getResults(query_TCGA)

#build a query to retrieve gene expression data
query_TCGA <- GDCquery(project = 'TCGA-LUAD',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = c())
getResults(query_TCGA)

#download
GDCdownload(query_TCGA)

#prepare data
tcga_luad_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
luad_rna <- assay(tcga_luad_data, 'unstranded')
