# RNA vs CNV analysis

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(DESeq2)

#Exploration of CNV and RNAseq data 
statRes_map_CNV = read.csv('~/model_fit_Python/model_results/results_luad/results_3/statRes_map_CNV.csv',header=TRUE)
statRes_map_noCNV = read.csv('~/model_fit_Python/model_results/results_brca/statRes_map_noCNV.csv',header=TRUE)
metadata = read.csv('model_fit_Python/model_data/metadata.csv',header=TRUE)

statRes_map_noCNV <- statRes_map_noCNV %>% remove_rownames %>% column_to_rownames(var="GeneID")
rna_counts <- rna_counts[(rownames(rna_counts) %in% rownames(res_allGenes)),]
res_rna_cnv <- res_rna_cnv %>% relocate(padj_2, .after = padj_1)

#res_allGenes$GeneID <- rownames(res_allGenes)
#res_allGenes <- res_allGenes %>% mutate(difference = B1_2 - B1_1)
#cnv <- cnv[(rownames(cnv) %in% rownames(res_allGenes)),]
#cnv <- cnv %>% mutate(cnv_mean = rowMeans(cnv))
#rna_tumor = read.delim("model_data/TCGA/lung_cancer/LUSC/s1/s1_rna_tumor.tsv", header=TRUE, sep="\t")
#rna <- read.delim("~/model_data/TCGA/lung_cancer/LUSC/s10/s10_rna_tumor.tsv", header=TRUE, sep="\t")
#cnv_brca <- brca_cnv_tumor %>% select(7:10, 13, 14, 18, 19, 22, 23, 25:27, 29, 34, 35, 37, 41, 45, 47, 49, 53, 54, 56:58, 61, 62, 64, 69, 72, 73, 76, 83, 90, 92:96, 102, 107, 109)
#cnv_3 <- cnv_brca %>% select(23,24,36)
#rna_norm_3 <- brca_rna_norm %>% select(55,57,93)
#rna_tum_3 <- brca_rna_tum %>% select(55,57,92)
#save(metadata, file = "~/model_data/TCGA/lung_cancer/LUAD/metadata.Rdata")
#write.csv(rna, file = "~/model_data/TCGA/lung_cancer/LUAD/rna_cnv_10.csv")

#Data manipulation and cleaning
stat_res_luad <- stat_res_luad %>% mutate(difference = B1_2 - B1_1)

resFit_merged <- cbind(statRes_map_noCNV, statRes_map_CNV) %>% 
  select() %>% 
  mutate(Difference = B1_2-B1_1) %>% 
  relocate(B0_2, .before = B1_2) %>% 
  remove_rownames %>% column_to_rownames(var="GeneID")

rna_tumor <- rna_tumor %>% select(2,4) %>% na.omit() %>% 
  colnames(rna_tumor)[2] <- "s1_tumor" %>% 
  rna_tumor[!duplicated(rna_tumor$GeneID), ] %>% #remove dublicates 
  remove_rownames %>% column_to_rownames(var="GeneID") 

rna_normal <- rna_normal %>% select(2,4) %>% na.omit() %>% 
  colnames(rna_normal)[2] <- "s1_tumor" %>% 
  rna_normal[!duplicated(rna_normal$GeneID), ] %>% #remove dublicates
  remove_rownames %>% column_to_rownames(var="GeneID")

cnv_tumor <- luad_cnv_tumor %>% replace(cnv_tumor, cnv_tumor>5, 6) %>% 
  mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  subset(cnv_tumor, cnv_mean <= 0.9 | cnv_mean >= 1.5) %>% 
  select(1:10) %>%
  remove_rownames %>% column_to_rownames(var="Row.names")

cnv_normal <- data.frame("TCGA-50-5932-11A" = rep(1, nrow(luad_cnv)), "TCGA-44-6147-11A" = rep(1, nrow(luad_cnv)),
                         "TCGA-55-6979-11A" = rep(1, nrow(luad_cnv)), "TCGA-50-5931-11A" = rep(1, nrow(luad_cnv)),
                         "TCGA-91-6835-11A" = rep(1, nrow(luad_cnv)), "TCGA-44-6776-11A" = rep(1, nrow(luad_cnv)),
                         "TCGA-44-6778-11A" = rep(1, nrow(luad_cnv)), "TCGA-44-6145-11A" = rep(1, nrow(luad_cnv)),
                         "TCGA-50-5939-11A" = rep(1, nrow(luad_cnv)), "TCGA-50-5935-11A" = rep(1, nrow(luad_cnv)))

#Metadata generation
metadata <- data.frame(patID = colnames(luad_rna), 
                       condition = rep(c("A", "B"), each = 10)) 
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="patID")  
metadata$condition <- as.factor(metadata$condition)

#Correction of normal RNAseq counts for CNV
cnv <- cnv/2
cnv <- cnv + 10e-9
rna <- luad_rna * cnv

#Counts normalization
load("~/model_data/TCGA/lung_cancer/LUSC/cnv_lusc.Rdata")
load("~/model_data/TCGA/lung_cancer/LUSC/rna_lusc.Rdata")

#Counts normalization
rna_normalized <- luad_rna %>%  as.matrix()
#rna_normalized <- DESeq2::varianceStabilizingTransformation(rna_normalized)
rna_log_normalized <- DESeq2::rlog(rna_normalized)

rna_zscore <- t(scale(t(rna_log_normalized)))
rna_zscore_normal <- rna_zscore %>% as.data.frame() %>% select(1:10)
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% select(11:20) 
rna_zscore_normal <- rna_zscore_normal %>% as.data.frame %>% mutate(rna_mean = rowMeans(rna_zscore_normal)) %>% select(11)
rna_zscore_tumor <- rna_zscore_tumor %>% as.data.frame %>% mutate(rna_mean = rowMeans(rna_zscore_tumor)) %>% select(11)

luad_cnv <- luad_cnv %>%
  mutate(cnv_mean = rowMeans(luad_cnv))

#CNV factorization
cnv <- cnv %>% 
  mutate(cnv = case_when(
    cnv_mean <= 0.5 ~ "0",
    cnv_mean > 0.5 & cnv_mean <= 1.6 ~ "1",
    cnv_mean > 1.6 & cnv_mean <= 2.4 ~ "2",
    cnv_mean > 2.4 & cnv_mean <= 3.5 ~ "3",
    cnv_mean > 3.5 & cnv_mean <= 4.5 ~ "4",
    cnv_mean > 4.5 ~ "5")) 

cnv <- cnv %>% 
  mutate(cn_group = case_when(
    cnv == "2" ~ "diploid",
    cnv == "0" ~ "cn_loss",
    cnv == "1" ~ "cn_loss",
    cnv == "3" ~ "cn_gain",
    cnv == "4" ~ "cn_gain",
    cnv == "5" ~ "cn_amplification"))

#gene group factorization
res_noCNV <- res_noCNV %>% 
  mutate(gene_group = case_when(
    B1_1 < -0.6 & padj < 0.05 ~ "DEG",
    B1_1 > 0.6 & padj < 0.05 ~ "DEG",
    B1_1 <= 0.6 & B1_1 >= -0.6 ~ "noDEG",
    B1_1 >= 0.6 & padj > 0.05 ~ "not significant",
    B1_1 <= -0.6 & padj > 0.05 ~ "not significant"))

#Gene group facrorization based on Effect size difference
deg <- deg %>% 
  mutate(gene_group = case_when(
    Difference <= -1.0 ~ "super-dosage",
    Difference >= 1.0 ~ "super-dosage",
    Difference > -1.0 & Difference < -0.4 ~ "dosage-sensitive",
    Difference < 1.0 & Difference > 0.4 ~ "dosage-sensitive",
    Difference >= -0.4 & Difference <= 0.4 ~ "dosage_insensitive"
  ))

#getting DE genes
sum(statRes_map_noCNV$padj < 0.05 & statRes_map_noCNV$B1_1 >= 0.6, na.rm=TRUE) #up_regulated
sum(statRes_map_noCNV$padj < 0.05 & statRes_map_noCNV$B1_1 <= -0.6, na.rm=TRUE) #down-reg

sum(res_allG$B1_1 < 1.0 & res_allG$B1_1 > -1 & res_allG$padj > 0.05, na.rm=TRUE)
sum(res_noCNV$gene_group == "not significant" & res_noCNV$cn_group == "cn_amplification", na.rm=TRUE)
sum(deg$gene_group == "dosage-insensitive" & deg$cnv == "5", na.rm=TRUE)

deg_up <- subset(stat_res_luad, stat_res_luad$padj_1 < 0.05 & stat_res_luad$B1_1 > 0.6)
deg_down <- subset(stat_res_luad, stat_res_luad$padj_1 < 0.05 & stat_res_luad$B1_1 < -0.6)
deg <- rbind(deg_up, deg_down)

#Calculate SD across rows
cnv_sd <- transform(cnv_filtered, sd=apply(cnv_filtered, 1, sd, na.rm=TRUE)) %>% 
  subset(cnv_sd, sd < 0.8) %>% 
  select(1,3,5:10)





