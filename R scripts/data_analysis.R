# RNA vs CNV analysis

library(ggplot2)
library(cowplot)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(DESeq2)

#Exploration of CNV and RNAseq data 
#Load datasets
#load("../data/patMeta_enc.RData")
#source("../code/utils.R")
statRes_map_CNV = read.csv('~/model_fit_Python/model_results/results_luad/results_3/statRes_map_CNV.csv',header=TRUE)
statRes_map_noCNV = read.csv('~/model_fit_Python/model_results/results_brca/statRes_map_noCNV.csv',header=TRUE)
res_nocnv2 = read.csv('~/model_fit_Python/model_results/results_sinthetic_brca/heterog_cnv/res_sint_nocnv.csv',header=TRUE)
rna_normal = read.csv('~/model_data/TCGA/breast_cancer/sinthetic_data/rna_normal.csv',header=TRUE)

#Data preprocessing
cnv_tumor <- cnv_tumor %>% remove_rownames %>% column_to_rownames(var="Row.names")
rna_norm <- rna_norm[(rownames(rna_norm) %in% rownames(cnv_tumor)),]
res_rna_cnv <- res_rna_cnv %>% relocate(padj_2, .after = padj_1)

cnv_tumor <- replace(cnv, cnv>5, 5)

#Data save
save(p6, file = "~/model_fit_Python/model_results/results_sinthetic_brca/p6.Rdata")
write.csv(rna_cnv, file = "~/model_fit_Python/model_data/test_tum_vs_norm_46/rna_cnv.csv")

#Data manipulation and cleaning
stat_res_luad <- stat_res_luad %>% mutate(difference = B1_2 - B1_1)

resFit_merged <- cbind(statRes_map_noCNV, statRes_map_CNV) %>% 
  select() %>% 
  mutate(Difference = B1_2-B1_1) %>% 
  relocate(B0_2, .before = B1_2) %>% 
  column_to_rownames(var="GeneID")

rna_tumor <- rna_tumor %>% select(2,4) %>% na.omit() %>% 
  colnames(rna_tumor)[2] <- "s1_tumor" %>% 
  rna_tumor[!duplicated(rna_tumor$GeneID), ] %>% #remove dublicates 
  column_to_rownames(var="GeneID") 

rna_normal <- rna_normal %>% select(2,4) %>% na.omit() %>% 
  colnames(rna_normal)[2] <- "s1_tumor" %>% 
  rna_normal[!duplicated(rna_normal$GeneID), ] %>% #remove dublicates
  column_to_rownames(var="GeneID")

cnv_tumor <- luad_cnv_tumor %>% replace(cnv_tumor, cnv_tumor>5, 6) %>% 
  mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  subset(cnv_tumor, cnv_mean <= 0.9 | cnv_mean >= 1.5) %>% 
  select(1:10) %>%
  column_to_rownames(var="Row.names")

#Metadata generation
rna_nor <- rna_normal %>% gsub("-", ".", rna_normal)
metadata <- data.frame(patID = colnames(rna_norm_tum), 
                       condition = rep(c("A", "B"), each = 45)) 
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="patID")  
metadata$condition <- as.factor(metadata$condition)
all.equal(colnames(),rownames())

#Generating simulated data
#starting from normal mRNA counts
rna_norm_tum <- cbind(brca_rna_norm, brca_cnv_tumor)
# remove genes with low counts
rna_norm_tum <- rna_norm_tum[rowSums(rna_norm_tum) > 500,] 

rna_norm <- rna %>% select(1:50)
rna_tum <- rna %>%  select(51:100)

group1 <- rna_normal[1:2000,]
group2 <- rna_normal[2001:4000,]
group3 <- rna_normal[4001:6000,]
group4 <- rna_normal[6001:8000,]
group5 <- rna_normal[8001:10000,]

cnv1 <- matrix(1, nrow(group1), 50)
cnv2 <- matrix(2, nrow(group1), 50)
cnv3 <- matrix(3, nrow(group1), 50)
cnv4 <- matrix(4, nrow(group1), 50)
cnv5 <- matrix(5, nrow(group1), 50)

cnv <- rbind(cnv1, cnv2, cnv3, cnv4, cnv5)
rna_normal <- rbind(group1, group2, group3, group4, group5)
rna_tumor <- rna_tumor[(rownames(rna_tumor) %in% rownames(rna_normal)),]
rna_norm_tum <- cbind(rna_norm, rna_tum)

rna_norm <- rna_norm_tum %>% select(1:110)
colnames(cnv_norm) <- colnames(rna_norm)
rownames(cnv_norm) <- rownames(rna_norm)

cnv_norm <- matrix(1, nrow(rna_norm), 45)
cnv <- cbind(cnv_tumor, cnv_norm)


#Correction of normal RNAseq counts for CNV
cnv_tumor <- cnv_tumor/2   
cnv <- cnv + 10e-9
rna_tumor <- rna_tumor * cnv_tum
rna_cnv <- rna_norm_tum * cnv


#Counts normalization
rna_normalized <- luad_rna %>%  as.matrix()
#rna_normalized <- DESeq2::varianceStabilizingTransformation(rna_normalized)
rna_log_normalized <- DESeq2::rlog(rna_normalized)

rna_zscore <- t(scale(t(rna_log_normalized)))
rna_zscore_normal <- rna_zscore %>% as.data.frame() %>% select(1:10)
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% select(11:20) 
rna_zscore_normal <- rna_zscore_normal %>% as.data.frame() %>%
  mutate(rna_mean = rowMeans(rna_zscore_normal)) %>% select(11)
rna_zscore_tumor <- rna_zscore_tumor %>% as.data.frame() %>% 
  mutate(rna_mean = rowMeans(rna_zscore_tumor)) %>% select(11)

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





