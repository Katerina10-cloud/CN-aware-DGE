### Download data from GDC portal ###

#rm(list=ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/")
pkgs <- c("tidyverse", "TCGAbiolinks", "SummarizedExperiment")
sapply(pkgs, require, character.only = TRUE)


# Get a list of projects #
gdcprojects = getGDCprojects()
getProjectSummary('TCGA-BRCA')

# build a query to retrieve data #
query_TCGA_cnv <- TCGAbiolinks::GDCquery(project = 'TCGA-LUSC',
                                         data.category = "Copy Number Variation",
                                         #sample.type = c("Recurrent Tumor", "Recurrent Tumor;Blood Derived Normal",
                                                         #"Blood Derived Normal;Recurrent Tumor"),
                                         #data.type = "Gene Level Copy Number")
                                         data.type = "Copy Number Segment")
                                         #workflow.type = 'ASCAT3'
                                         #workflow.type = "DNAcopy")
    

#output_query_TCGA <- TCGAbiolinks::getResults(query_TCGA_cnv)

#output_query_TCGA <- output_query_TCGA %>% 
  #dplyr::filter(sample_type=="Primary Tumor;Blood Derived Normal")

query_TCGA_cnv_rec[[1]][[1]] <- query_TCGA_cnv_rec[[1]][[1]] %>% 
  #dplyr::filter(sample_type=="Primary Tumor;Solid Tissue Normal")
  dplyr::filter(sample_type==c("Recurrent Tumor", "Recurrent Tumor;Blood Derived Normal",
                               "Blood Derived Normal;Recurrent Tumor"))


# build a query to retrieve Tumor gene expression data #
query_TCGA_rna_tumor <- TCGAbiolinks::GDCquery(project = 'TCGA-LGG',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor",
                       access = 'open')
              

# build a query to retrieve Nor2mal gene expression data #
query_TCGA_rna_rec <- TCGAbiolinks::GDCquery(project = 'TCGA-LGG',
                           data.category = 'Transcriptome Profiling',
                           experimental.strategy = 'RNA-Seq',
                           workflow.type = 'STAR - Counts',
                           data.type = "Gene Expression Quantification",
                           sample.type = "Recurrent Tumor",
                           access = 'open')
                      
# download data
GDCdownload(query_TCGA_cnv)
GDCdownload(query_TCGA_cnv_rec)
GDCdownload(query_TCGA_rna_tumor)
GDCdownload(query_TCGA_rna_rec)


# Prepare data #
#query_TCGA_cnv[[1]][[1]]$sample.submitter_id <- substr(query_TCGA_cnv[[1]][[1]]$sample.submitter_id, 1, 20)

cnv_tumor <- TCGAbiolinks::GDCprepare(query_TCGA_cnv, summarizedExperiment = T)
cnv_tum <- assay(cnv_tumor, 'copy_number', rownames = TRUE) %>% as.data.frame()
gene_name <- as.data.frame(cnv_tumor_rec@rowRanges@elementMetadata@listData[["gene_name"]])
colnames(gene_name)[1] <- "GeneID"
cnv_tum_rec <- cbind(gene_name, cnv_tum_rec)

#cnv_tum <- cnv_tumor %>% select(-c(gene_id, chromosome, start, end)) %>% 
  #na.omit()

cnv_tum_rec <- cnv_tum_rec[!duplicated(cnv_tum_rec$GeneID), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="GeneID") %>%
  na.omit()

colnames(cnv_tum_rec) <- substr(colnames(cnv_tum_rec), 1, 12)
#cnv_tum <- cnv_tum[, !duplicated(colnames(cnv_tum))]

rna_normal <- TCGAbiolinks::GDCprepare(query_TCGA_rna_rec, summarizedExperiment = TRUE)
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


## Clinical data ##
colnames(rna_tum) <- substr(colnames(rna_tum), 1, 12)
patient_barcodes <- c(colnames(rna_tum))

clinical_query <- TCGAbiolinks::GDCquery(project = "TCGA-LUSC", data.category = "Clinical",
                                          data.format = "bcr xml", barcode = patient_barcodes)
GDCdownload(clinical_query)

clinical_data <- TCGAbiolinks::GDCprepare_clinic(clinical_query, clinical.info = "patient")

saveRDS(clinical_data, file = "brca/clinical_full.RDS")


# Data preprocessing #
rna_tum <- rna_tum[,(colnames(rna_tum) %in% colnames(cnv_tumor))]
rna_normal<- rna_normal[,(colnames(rna_norm) %in% colnames(cnv_tumor))]
cnv_tumor <- cnv_tumor[,(colnames(cnv_tumor) %in% colnames(rna_tum))]

#Filtering low counts genes
#lihc_rna_tum <- lihc_rna_tum[which(rowSums(lihc_rna_tum)>50),]

# Reordering to match datasets #
colnames_idx <- match(colnames(rna_tum), colnames(rna_norm))
rna_norm <- rna_norm[,colnames_idx]
colnames_idx <- match(colnames(cnv_tum), colnames(rna_norm))
rna_norm <- rna_norm[,colnames_idx]

rownames_idx <- match(rownames(cnv_tum), rownames(rna_tum))
rna_tum <- rna_tum[rownames_idx,] %>% na.omit()

x <- colnames(rna_norm)
names(rna_norm) <- paste(x,"-11A")

x <- colnames(rna_tum)
names(rna_tum) <- paste(x,"-01A")

x <- colnames(cnv_tumor)
names(cnv_tumor) <- paste(x,"-11A")

saveRDS(cnv_tum_rec, file = "glioma/recurrent/cnv_tumor_rec.RDS")
saveRDS(rna_norm, file = "glioma/recurrent/rna_tumor_rec.RDS")
saveRDS(rna_tum, file = "glioma/primary/rna_tumor.RDS")

# Formatting strings | Select samples | CN segment data
cnv_seg$Sample <- stringr::str_sub(cnv_seg$Sample,1,12)
colnames(cnv_tumor) <- stringr::str_sub(colnames(cnv_tumor),1,12)
cnv_seg <- cnv_seg[cnv_seg$Sample %in% colnames(cnv_tumor),]

cnv_seg <- cnv_seg %>%
  dplyr::rename(
    chr = Chromosome,
    start = Start,
    end = End,
    seg_mean = Segment_Mean,
    sample = Sample) %>%
  dplyr::select(chr, start, end, seg_mean, sample) %>%
  as.data.frame()

save(lihc_cnv_seg, file = "liver_LIHC/lihc_cnv_seg.Rdata")
