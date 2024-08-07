###----------------------------------------------------------###
### RNA vs CNV analysis ###
###----------------------------------------------------------###

library(ggplot2)
library(cowplot)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(DESeq2)

#Exploration of CNV and RNAseq data 
#source("../code/utils.R")

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/relationship_plot")

res2 = read.csv('res_sim_cnv.csv',header=TRUE)


#Data preprocessing
cnv_tumor <- cnv_tumor %>% remove_rownames %>% column_to_rownames(var="Row.names")
rna_norm <- rna_norm[(rownames(rna_norm) %in% rownames(cnv_tumor)),]
res_rna_cnv <- res_rna_cnv %>% relocate(padj_2, .after = padj_1)

#cnv_tumor <- replace(cnv, cnv>5, 5)

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

#Formatter of column names#
colnames(luad_rna_sd) <- gsub(pattern = "\\.", replacement = "-", colnames(luad_rna_sd))

luad_tum <- luad_rna_sd %>% select(46:90)
luad_norm <- luad_rna_sd %>% select(1:45)

colnames(lihc_cnv) <- stringr::str_sub(colnames(lihc_cnv),1,12)
lihc_tum <- lihc_rna_tum[,colnames(lihc_rna_tum) %in% colnames(lihc_cnv)]
lihc_norm <- lihc_rna_norm[,colnames(lihc_rna_norm) %in% colnames(lihc_cnv)]

x <- colnames(lihc_norm)
names(lihc_norm) <- paste(x,"-11A")
lihc_rna <- cbind(lihc_norm, lihc_tum)

save(lihc_rna, file = "lihc_rna.Rdata")

#Metadata generation
rna_nor <- rna_normal %>% gsub("-", ".", rna_normal)
rna_counts <- cbind(brca_rna_normal, brca_rna_tum)
metadata <- data.frame(patID = colnames(rna_nocnv), 
                       condition = rep(c("A", "B"), each = 96)) 
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="patID")  
metadata$condition <- as.factor(metadata$condition)
all.equal(colnames(),rownames())

#Generating simulated data
#starting from normal mRNA counts
rna_norm_tum <- cbind(brca_rna_norm, brca_cnv_tumor)

# remove genes with low counts
luad_rna_norm <- luad_rna_norm[rowSums(luad_rna_norm) > 200,] 


#Correction of normal RNAseq counts for CNV
cnv_tumor <- cnv_tumor/2   
cnv <- cnv + 10e-9
rna_tumor <- rna_tumor * cnv_tum
rna_cnv <- rna_counts * cnv

#Counts normalization
brca_cnv <- apply(brca_cnv, 2, function(x) ifelse(x > 10, 10, x)) 
rna_normalized <- lihc_rna %>%  as.matrix()
rna_log_normalized <- DESeq2::rlog(rna_normalized)

# Z-score trasformation
rna_zscore <- t(scale(t(rna_log_normalized)))
rna_zscore_normal <- rna_zscore %>% as.data.frame() %>% select(1:15)
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% select(16:30) 
rna_zscore_normal <- rna_zscore_normal %>% as.data.frame() %>%
  mutate(rna_mean = rowMeans(rna_zscore_normal)) %>% select(rna_mean)
rna_zscore_tumor <- rna_zscore_tumor %>% as.data.frame() %>% 
  mutate(rna_mean = rowMeans(rna_zscore_tumor)) %>% select(rna_mean)

lihc_cnv <- lihc_cnv %>% as.data.frame() %>% 
  mutate(cnv_mean = rowMeans(lihc_cnv))


#CNV factorization
cnv <- lihc_cnv %>% 
  mutate(cnv = case_when(
    cnv_mean <= 0.5 ~ "0",
    cnv_mean > 0.5 & cnv_mean <= 1.8 ~ "1",
    cnv_mean > 1.8 & cnv_mean <= 2.4 ~ "2",
    cnv_mean > 2.4 & cnv_mean <= 3.5 ~ "3",
    cnv_mean > 3.5 & cnv_mean <= 4.4 ~ "4",
    cnv_mean > 4.4 ~ "5")) %>% 
  select(cnv)

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
sum(res1_nocnv$padj < 0.05 & res1_nocnv$log2FoldChange >= 0.5, na.rm=TRUE) #up_regulated
sum(res1_nocnv$padj < 0.05 & res1_nocnv$log2FoldChange <= -0.5, na.rm=TRUE) #down-reg

sum(res_allG$B1_1 < 1.0 & res_allG$B1_1 > -1 & res_allG$padj > 0.05, na.rm=TRUE)
sum(res_noCNV$gene_group == "not significant" & res_noCNV$cn_group == "cn_amplification", na.rm=TRUE)
sum(deg$gene_group == "dosage-insensitive" & deg$cnv == "5", na.rm=TRUE)

deg_up <- subset(res1_nocnv, res1_nocnv$padj < 0.05 & res1_nocnv$log2FoldChange > 0.5)
deg_down <- subset(res1_nocnv, res1_nocnv$padj < 0.05 & res1_nocnv$log2FoldChange < -0.5)
deg <- rbind(deg_up, deg_down)
deg_cnv <- rownames(deg)

#Separate CNV regulated genes
genes_cnvReg_1 <- subset(deg, deg$padj_2 > 0.05 & deg$difference > 0.2)
genes_cnvReg_2 <- subset(deg, deg$padj_2 > 0.05 & deg$difference < -0.2)
genes_cnvReg <- rbind(genes_cnvReg_1, genes_cnvReg_2)


#Calculate SD across rows and columns 
library(matrixStats)

rna_counts <- transform(rna_counts, sd=apply(rna_counts, 1, sd, na.rm=TRUE)) %>% 
  subset(sd > 25.0) %>% 
  select(1:192) %>% 
  as.matrix()

cnv_tumor <- as.matrix(lihc_cnv_tum)
cnv_tumor_sd <- matrixStats::colSds(cnv_tumor) %>% 
  as.data.frame() %>% 
  setNames("sd")  %>% 
  subset(sd > 1)

brca_rna_normal <- rna_counts[,1:96]
brca_rna_tum <- rna_counts[,97:192]

cnv_tumor <- cnv_tumor[,colnames(cnv_tumor) %in% rownames(cnv_tumor_sd)]
colnames(cnv_tumor) <- colnames(brca_rna_tum)


###----------------------------------------------------###
### EdgeR DE test ###
###----------------------------------------------------###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/de_fit_Python/data_simulation/")
library(edgeR)

rna_nocnv = read.csv('sim3_real_data/brca/rna_nocnv.csv', header=TRUE)

metadata <- read.csv('metadata.csv', header=TRUE)
rna_nocnv <- rna_nocnv %>% remove_rownames %>% column_to_rownames(var="X")

# Round counts #
rna_cnv <- ceiling(rna_cnv)
rna_cnv <- ceiling(rna_cnv)

#Creating a DGEList object
gList <- DGEList(counts=rna_cnv, genes=rownames(rna_cnv))

#Normalization
gList <- calcNormFactors(gList, method="TMM")

# Design matrix
designMat <- model.matrix(~ condition, metadata)

# Estimating dispersions
gList <- estimateGLMCommonDisp(gList, design=designMat)
gList <- estimateGLMTrendedDisp(gList, design=designMat)
gList <- estimateGLMTagwiseDisp(gList, design=designMat)

# Model fit
fit <- glmFit(gList, designMat)
lrt <- glmLRT(fit, coef=ncol(fit$design), contrast = NULL)
edgeR_result <- topTags(lrt)
res_edger <- topTags(lrt, n=15000)$table                       

save(res_edger, file = "results/sim_3/brca/res1_edger.Rdata")
