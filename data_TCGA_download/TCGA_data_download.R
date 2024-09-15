### Download data from GDC portal ###

#rm(list=ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/")
pkgs <- c("tidyverse", "TCGAbiolinks", "SummarizedExperiment")
sapply(pkgs, require, character.only = TRUE)


# Get a list of projects #
gdcprojects = getGDCprojects()
#getProjectSummary('TCGA-LGG')

# build a query to retrieve data #
query_TCGA_cnv <- TCGAbiolinks::GDCquery(project = 'TCGA-KIRC',
        data.category = 'Copy Number Variation',
        #sample.type = c("Primary Tumor"),
        data.type = "Gene Level Copy Number",
        workflow.type = 'ASCAT3')
    

#output_query_TCGA <- TCGAbiolinks::getResults(query_TCGA_cnv)

#output_query_TCGA <- output_query_TCGA %>% 
  #dplyr::filter(sample_type=="Primary Tumor;Blood Derived Normal")
query_TCGA_cnv[[1]][[1]] <- query_TCGA_cnv[[1]][[1]] %>% 
  dplyr::filter(sample_type=="Primary Tumor;Solid Tissue Normal")
  


# build a query to retrieve Tumor gene expression data #
query_TCGA_rna_tumor <- TCGAbiolinks::GDCquery(project = 'TCGA-KIRC',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor",
                       access = 'open')
              

# build a query to retrieve Normal gene expression data #
query_TCGA_rna_normal <- TCGAbiolinks::GDCquery(project = 'TCGA-KIRC',
                           data.category = 'Transcriptome Profiling',
                           experimental.strategy = 'RNA-Seq',
                           workflow.type = 'STAR - Counts',
                           data.type = "Gene Expression Quantification",
                           sample.type = "Solid Tissue Normal",
                           access = 'open')
                      
# download data
GDCdownload(query_TCGA_cnv)
GDCdownload(query_TCGA_rna_tumor)
GDCdownload(query_TCGA_rna_normal)


# Prepare data #
cnv_tumor <- TCGAbiolinks::GDCprepare(query_TCGA_cnv, summarizedExperiment = TRUE)
cnv_tum <- assay(cnv_tumor, 'copy_number', rownames = TRUE) %>% as.data.frame()
gene_name <- as.data.frame(cnv_tumor@rowRanges@elementMetadata@listData[["gene_name"]])
colnames(gene_name)[1] <- "GeneID"
cnv_tum <- cbind(gene_name, cnv_tum)

cnv_tum <- cnv_tum[!duplicated(cnv_tum$GeneID), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="GeneID") %>%
  na.omit()

rna_normal <- TCGAbiolinks::GDCprepare(query_TCGA_rna_normal, summarizedExperiment = TRUE)
rna_norm <- assay(rna_normal, 'unstranded', rownames = TRUE)
gene_name <- as.data.frame(rna_normal@rowRanges@elementMetadata@listData[["gene_name"]]) 
colnames(gene_name)[1] <- "GeneID"
rna_norm <- cbind(gene_name, rna_norm)

rna_norm <- rna_norm[!duplicated(rna_norm$GeneID), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="GeneID") %>%
  na.omit()


rna_tumor <- TCGAbiolinks::GDCprepare(query_TCGA_rna_tumor, summarizedExperiment = TRUE)
rna_tum <- assay(rna_tumor, 'unstranded', rownames = TRUE) %>% as.data.frame()
gene_name <- as.data.frame(rna_tumor@rowRanges@elementMetadata@listData[["gene_name"]]) 
colnames(gene_name)[1] <- "GeneID"
rna_tum <- cbind(gene_name, rna_tum)

rna_tum <- rna_tum[!duplicated(rna_tum$GeneID), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="GeneID") %>%
  na.omit()

#substring columns
colnames(rna_norm) <- substr(colnames(rna_norm), 1, 12)
colnames(rna_tum) <- substr(colnames(rna_tum), 1, 12)
colnames(cnv_tum) <- substr(colnames(cnv_tum), 1, 12)

# clinical data
clinical <- TCGAbiolinks::GDCquery(project = "TCGA-COAD", data.category = "Clinical", 
                          data.format = "bcr xml")
GDCdownload(clinical)

clinical <- TCGAbiolinks::GDCprepare_clinic(clinical, clinical.info = "patient")

saveRDS(clinical, file = "colon/clinical_coad.RDS")


# Data preprocessing #
rna_tum <- rna_tum[,(colnames(rna_tum) %in% colnames(cnv_tum))]
rna_norm <- rna_norm[,(colnames(rna_norm) %in% colnames(cnv_tum))]
cnv_tum <- cnv_tum[,(colnames(cnv_tum) %in% colnames(rna_norm))]

#Filtering low counts genes
#lihc_rna_tum <- lihc_rna_tum[which(rowSums(lihc_rna_tum)>50),]

# Reordering to match datasets #
colnames_idx <- match(colnames(rna_tum), colnames(cnv_tum))
cnv_tum <- cnv_tum[,colnames_idx]
colnames_idx <- match(colnames(cnv_tum), colnames(rna_norm))
rna_norm <- rna_norm[,colnames_idx]

  
rownames_idx <- match(rownames(rna_tum), rownames(cnv_tum))
cnv_tum <- cnv_tum[rownames_idx,]

x <- colnames(rna_norm)
names(rna_norm) <- paste(x,"-11A")

x <- colnames(rna_tum)
names(rna_tum) <- paste(x,"-01A")

x <- colnames(cnv_tum)
names(cnv_tum) <- paste(x,"-11A")

saveRDS(cnv_tum, file = "kidney/cnv_tumor.RDS")
saveRDS(rna_norm, file = "kidney/rna_norm.RDS")
saveRDS(rna_tum, file = "kidney/rna_tum.RDS")

# Formatting strings | Select samples
lihc_cnv_seg$Sample <- stringr::str_sub(lihc_cnv_seg$Sample,1,12)
colnames(lihc_cnv_tum) <- stringr::str_sub(colnames(lihc_cnv_tum),1,12)
lihc_cnv_seg <- lihc_cnv_seg[lihc_cnv_seg$Sample %in% colnames(lihc_cnv_tum),]

lihc_cnv_seg <- lihc_cnv_seg %>%
  dplyr::rename(
    chr = Chromosome,
    start = Start,
    end = End,
    seg_mean = Segment_Mean,
    sample = Sample) %>%
  dplyr::select(chr, start, end, seg_mean, sample) %>%
  as.data.frame()

save(lihc_cnv_seg, file = "liver_LIHC/lihc_cnv_seg.Rdata")
