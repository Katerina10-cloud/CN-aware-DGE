### Download data from GDC portal ###

rm(list=ls())
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "TCGAbiolinks", "SummarizedExperiment")
sapply(pkgs, require, character.only = TRUE)


# Get a list of projects #
gdcprojects = getGDCprojects()
getProjectSummary('TCGA-BRCA')

# build a query to retrieve data #
query_TCGA_cnv <- TCGAbiolinks::GDCquery(project = 'TCGA-LIHC',
        data.category = 'Copy Number Variation',
        sample.type = "Primary Tumor",
        data.type = "Copy Number Segment",
        workflow.type = 'DNAcopy')
        #barcode = luad_cnv_list)

output_query_TCGA <- TCGAbiolins::getResults(query_TCGA_cnv)


# build a query to retrieve Tumor gene expression data #
query_TCGA_rna_tumor <- TCGAbiolinks::GDCquery(project = 'TCGA-LIHC',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor",
                       access = 'open')
              

# build a query to retrieve Normal gene expression data #
query_TCGA_rna_normal <- TCGAbiolinks::GDCquery(project = 'TCGA-LIHC',
                           data.category = 'Transcriptome Profiling',
                           experimental.strategy = 'RNA-Seq',
                           workflow.type = 'STAR - Counts',
                           data.type = "Gene Expression Quantification",
                           sample.type = "Solid Tissue Normal",
                           access = 'open')
                      
# download data
GDCdownload(query_TCGA_cnv)

# clinical data download #
clinical_luad <- TCGAbiolinks::GDCquery(project = "TCGA-LAML", data.category = "Clinical", 
                          data.format = "bcr xml", barcode = barcode_tumor)
GDCdownload(clinical_luad)

clinical_luad <- TCGAbiolinks::GDCprepare_clinic(clinical_luad, clinical.info = "patient")
clinical_luad <- clinical_luad %>% select(1,6,7,69)

# Prepare data #
lihc_cnv_seg <- TCGAbiolinks::GDCprepare(query_TCGA_cnv)

lihc_cnv_tumor <- TCGAbiolinks::GDCprepare(query_TCGA_cnv, summarizedExperiment = TRUE)
lihc_cnv_tum <- assay(lihc_cnv_tumor, 'copy_number', rownames = TRUE) %>% as.data.frame()
gene_name <- as.data.frame(lihc_cnv_tumor@rowRanges@elementMetadata@listData[["gene_name"]])
colnames(gene_name)[1] <- "GeneID"
lihc_cnv_tum <- cbind(gene_name, lihc_cnv_tum)

lihc_cnv_tum <- lihc_cnv_tum[!duplicated(lihc_cnv_tum$GeneID), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="GeneID") %>%
  na.omit()

lihc_rna_normal <- TCGAbiolinks::GDCprepare(query_TCGA_rna_normal, summarizedExperiment = TRUE)
lihc_rna_norm <- assay(lihc_rna_normal, 'unstranded', rownames = TRUE)
gene_name <- as.data.frame(lihc_rna_normal@rowRanges@elementMetadata@listData[["gene_name"]]) 
colnames(gene_name)[1] <- "GeneID"
lihc_rna_norm <- cbind(gene_name, lihc_rna_norm)

lihc_rna_norm <- lihc_rna_norm[!duplicated(lihc_rna_norm$GeneID), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="GeneID") %>%
  na.omit()

#substring columns
colnames(lihc_rna_tum) <- substr(colnames(lihc_rna_tum), 1, 12)
colnames(lihc_cnv_tum) <- substr(colnames(lihc_cnv_tum), 1, 12)

# clinical data
clinical_brca <- TCGAbiolinks::GDCquery(project = "TCGA-BRCA", data.category = "Clinical", 
                          data.format = "bcr xml", barcode = barcode_tumor)
GDCdownload(clinical_brca)

clinical_brca <- TCGAbiolinks::GDCprepare_clinic(clinical_brca, clinical.info = "patient")
clinical_brca <- clinical_brca %>% select(1,6,22,107)
save(clinical_brca, file = "~/model_data/TCGA/breast_cancer/data/clinical_brca.Rdata")

lihc_rna_norm <- lihc_rna_norm[(rownames(lihc_rna_norm) %in% rownames(lihc_cnv_tum)),]
lihc_cnv_tum <- lihc_cnv_tum[,(colnames(lihc_cnv_tum) %in% colnames(lihc_rna_tum))]


#Filtering low counts genes
lihc_rna_tum <- lihc_rna_tum[which(rowSums(lihc_rna_tum)>50),]

# Reordering to match datasets #
colnames_idx <- match(colnames(lihc_rna_tum), colnames(lihc_cnv_tum))
lihc_cnv_tum <- lihc_cnv_tum[,colnames_idx]

rownames_idx <- match(rownames(lihc_rna_tum), rownames(lihc_cnv_tum))
lihc_cnv_tum <- lihc_cnv_tum[rownames_idx,]

#brca_rna_tum <- brca_rna_tum[order(match(rownames(brca_rna_tum), rownames(cnv_tumor))),]

x <- colnames(lihc_rna_norm)
names(lihc_rna_norm) <- paste(x,"-11A")

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
