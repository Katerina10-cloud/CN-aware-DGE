BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

#Get a list of projects
gdcprojects = getGDCprojects()
getProjectSummary('TCGA-LUAD')

##build a query to retrieve Copy Number data
query_TCGA <- GDCquery(project = 'TCGA-LUAD',
        data.category = 'Copy Number Variation',
        #sample.type = "Primary Tumor",
        data.type = "Gene Level Copy Number",
        workflow.type = 'ASCAT3',
        barcode = c("TCGA-50-5932-01A-11D-1752-01", "TCGA-49-6742-01A-11D-1854-01",
                    "TCGA-44-6147-01A-11D-A273-01", "TCGA-55-6979-01A-11D-1943-01",
                    "TCGA-50-5931-01A-11D-1752-01", "TCGA-38-4626-01A-01D-1549-01",
                    "TCGA-91-6835-01A-11D-1854-01", "TCGA-49-4490-01A-21D-1854-01",
                    "TCGA-55-6970-01A-11D-1943-01", "TCGA-55-6969-01A-11D-1943-01",
                    "TCGA-91-6831-01A-11D-1854-01", "TCGA-50-5936-01A-11D-1624-01",
                    "TCGA-44-2665-01A-01D-0944-01", "TCGA-44-6776-01A-11D-1854-01",
                    "TCGA-55-6978-01A-11D-1943-01", "TCGA-49-6744-01A-11D-1854-01",
                    "TCGA-55-6982-01A-11D-1943-01", "TCGA-44-3396-01A-01D-1204-01",
                    "TCGA-49-6761-01A-31D-1943-01", "TCGA-44-6778-01A-11D-1854-01",
                    "TCGA-44-6145-01A-11D-1752-01", "TCGA-50-6595-01A-12D-1854-01",
                    "TCGA-44-6148-01A-11D-1752-01", "TCGA-91-6828-01A-11D-1854-01",
                    "TCGA-91-6849-01A-11D-1943-01", "TCGA-91-6847-01A-11D-1943-01",
                    "TCGA-44-2668-01A-01D-0944-01", "TCGA-44-2662-01A-01D-0944-01",
                    "TCGA-44-5645-01A-01D-1624-01", "TCGA-55-6985-01A-11D-1943-01",
                    "TCGA-44-2657-01A-01D-1549-01", "TCGA-73-4676-01A-01D-1752-01",
                    "TCGA-50-5939-01A-11D-1624-01", "TCGA-50-5935-01A-11D-1752-01",
                    "TCGA-55-6981-01A-11D-1943-01", "TCGA-91-6829-01A-21D-1854-01",
                    "TCGA-44-6777-01A-11D-1854-01", "TCGA-55-6983-01A-11D-1943-01",
                    "TCGA-91-6836-01A-21D-1854-01", "TCGA-55-6986-01A-11D-1943-01",
                    "TCGA-38-4625-01A-01D-1204-01", "TCGA-50-5933-01A-11D-1752-01",
                    "TCGA-38-4632-01A-01D-1752-01", "TCGA-55-6971-01A-11D-1943-01",
                    "TCGA-55-6975-01A-11D-1943-01", "TCGA-49-6743-01A-11D-1854-01",
                    "TCGA-49-4512-01A-21D-1854-01", "TCGA-55-6972-01A-11D-1943-01",
                    "TCGA-38-4627-01A-01D-1204-01", "TCGA-44-3398-01A-01D-1549-01",
                    "TCGA-44-2655-01A-01D-0944-01", "TCGA-44-6146-01A-11D-1752-01",
                    "TCGA-44-2661-01A-01D-1549-01", "TCGA-49-6745-01A-11D-1854-01"))

output_query_TCGA <- getResults(query_TCGA)

#build a query to retrieve gene expression data
query_TCGA <- GDCquery(project = 'TCGA-LUAD',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor",
                       access = 'open',
                       barcode = c("TCGA-50-5932-01A-11R-1755-07", "TCGA-49-6742-01A-11R-1858-07",
                                   "TCGA-44-6147-01A-11R-1755-07", "TCGA-55-6979-01A-11R-1949-07",
                                   "TCGA-50-5931-01A-11R-1755-07", "TCGA-38-4626-01A-01R-1206-07",
                                   "TCGA-91-6835-01A-11R-1858-07", "TCGA-50-5930-01A-11R-1755-07",
                                   "TCGA-49-4490-01A-21R-1858-07", "TCGA-55-6970-01A-11R-1949-07",
                                   "TCGA-55-6969-01A-11R-1949-07", "TCGA-91-6831-01A-11R-1858-07",
                                   "TCGA-50-5936-01A-11R-1628-07", "TCGA-44-2665-01A-01R-0946-07",
                                   "TCGA-44-6776-01A-11R-1858-07", "TCGA-55-6978-01A-11R-1949-07",
                                   "TCGA-49-6744-01A-11R-1858-07", "TCGA-55-6982-01A-11R-1949-07",
                                   "TCGA-44-3396-01A-01R-1206-07", "TCGA-49-6761-01A-31R-1949-07",
                                   "TCGA-44-6778-01A-11R-1858-07", "TCGA-44-6145-01A-11R-1755-07",
                                   "TCGA-50-6595-01A-12R-1858-07", "TCGA-44-6148-01A-11R-1755-07",
                                   "TCGA-91-6828-01A-11R-1858-07", "TCGA-91-6849-01A-11R-1949-07",
                                   "TCGA-91-6847-01A-11R-1949-07", "TCGA-44-2668-01A-01R-A278-07",
                                   "TCGA-44-2662-01A-01R-0946-07", "TCGA-44-5645-01A-01R-1628-07",
                                   "TCGA-55-6985-01A-11R-1949-07", "TCGA-55-6984-01A-11R-1949-07",
                                   "TCGA-44-2657-01A-01R-1107-07", "TCGA-73-4676-01A-01R-1755-07",
                                   "TCGA-50-5939-01A-11R-1628-07", "TCGA-50-5935-01A-11R-1755-07",
                                   "TCGA-55-6981-01A-11R-1949-07", "TCGA-91-6829-01A-21R-1858-07",
                                   "TCGA-55-6980-01A-11R-1949-07", "TCGA-44-6777-01A-11R-1858-07",
                                   "TCGA-55-6983-01A-11R-1949-07", "TCGA-91-6836-01A-21R-1858-07",
                                   "TCGA-55-6986-01A-11R-1949-07", "TCGA-38-4625-01A-01R-1206-07",
                                   "TCGA-50-5933-01A-11R-1755-07", "TCGA-38-4632-01A-01R-1755-07",
                                   "TCGA-55-6971-01A-11R-1949-07", "TCGA-55-6975-01A-11R-1949-07",
                                   "TCGA-49-6743-01A-11R-1858-07", "TCGA-49-4512-01A-21R-1858-07",
                                   "TCGA-55-6972-01A-11R-1949-07", "TCGA-38-4627-01A-01R-1206-07",
                                   "TCGA-44-3398-01A-01R-1107-07", "TCGA-44-2655-01A-01R-0946-07",
                                   "TCGA-44-6146-01A-11R-A278-07", "TCGA-44-2661-01A-01R-1107-07",
                                   "TCGA-49-6745-01A-11R-1858-07"))

##build a query to retrieve Copy Number data
getResults(query_TCGA)


#download
GDCdownload(query_TCGA)

#prepare data
tcga_luad_cnv <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
luad_cnv_tumor <- assay(tcga_luad_cnv, 'copy_number', rownames = TRUE)

gene_name <- as.data.frame(tcga_luad_cnv@rowRanges@elementMetadata@listData[["gene_name"]]) %>% 
  colnames(gene_name)[1] <- "GeneID"
luad_cnv_tumor <- as.data.frame(luad_cnv_tumor)
luad_cnv_tumor <- cbind(gene_name, luad_cnv_tumor)
luad_cnv_tumor <- luad_cnv_tumor %>%  luad_cnv_tumor[!duplicated(luad_cnv_tumor$GeneID), ] %>% 
  remove_rownames %>% column_to_rownames(var="GeneID") %>% na.omit(luad_cnv_tumor)

#common_cols <- intersect(colnames(luad_rna_normal), colnames(luad_rna_tumor))
luad_rna_normal <- luad_rna_normal[, c(-8,-14,-28,-30,-32,-33,-39,-40,-42,-47,-55,-57,-58)]

colnames(luad_rna_normal) <- c('s1_normal','s2_normal','s3_normal','s4_normal','s5_normal', 's6_normal',
                              's7_normal','s8_normal','s9_normal','s10_normal','s11_normal', 's12_normal',
                              's13_normal','s14_normal','s15_normal','s16_normal','s17_normal', 's18_normal',
                              's19_normal','s20_normal','s21_normal','s22_normal','s23_normal', 's24_normal',
                              's25_normal','s26_normal','s27_normal','s28_normal','s29_normal', 's30_normal',
                              's31_normal','s32_normal','s33_normal','s34_normal','s35_normal', 's36_normal',
                              's37_normal','s38_normal','s39_normal','s40_normal','s41_normal', 's42_normal',
                              's43_normal','s44_normal','s45_normal','s46_normal')
#rownames(countdata) <- seqdata[,1]
#substr("ThisIsAString", start=1, stop=5) #shorten sample names
#colnames(countdata) <- substr(colnames(countdata), 1, 7)
#colnames(countdata)==sampleinfo$SampleName
#all(colnames(countdata)==sampleinfo$SampleName)

#Filtering low counts genes
rna_normal_tumor_3 <- rna_normal_tumor_3[which(rowSums(rna_normal_tumor_3)>10),]
cnv <- cnv[(rownames(cnv) %in% rownames(rna_normal_tumor_3)),] #delete rows by name

save(luad_cnv_tumor, file = "~/model_data/TCGA/lung_cancer/cnv_luad.Rdata")
write.csv(rna_cnv, "model_data/TCGA/lung_cancer/last_test/rna_cnv.csv")

luad_cnv_tumor <- luad_cnv_tumor %>% remove_rownames %>% column_to_rownames(var="GeneID")
rna_normal <- rna_normal_tumor %>% select(1:46)
rna_tumor <- rna_normal_tumor %>% select(47:92)
cnv_tumor_3 <- luad_cnv_tumor %>% select(8,16,27)
rna_normal_tumor_3 <- rna_normal_tumor %>% as.data.frame() %>% select(8,16,27,54,62,73)
rna_normal_3 <- rna_normalized_normal %>% as.data.frame() %>% select(8,16,27)
rna_3 <- cbind(rna_tumor_3, rna_normal_3)

