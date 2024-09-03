### Identify homogeneous groups of patients based on their CNV profiles ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "DESeq2")
sapply(pkgs, require, character.only = TRUE)

source("CN-aware-DGE/R/utils.R")
set.seed(12345)

# Input data
data_path <- "TCGA/lung_cancer/LUAD/cnv_tumor.RDS"
dataset_name <- "LUAD_cnv"

# Clustering patients 
clustering_res <- clustering_patients_cnv(dataset_name, data_path)
clustering_res[[2]]
cnv_tumor <- clustering_res[[1]]

# Cluster selection and analysis
cnv_filt <- cnv_tumor[cnv_tumor$Cluster %in% c(2),]
cnv_filt <- as.matrix(t(cnv_filt))
cnv_filt <- apply(cnv_filt, 2, function(x) ifelse(x > 10, 10, x)) 

hist(rowMeans(cnv_filt),
     main = "LUAD - cluster 3", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 8)

cnv_mean <- cnv_filt %>% 
  as.data.frame() %>% 
  mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  select(cnv_mean)


## RNA data processing ##

data_path <- "TCGA/lung_cancer/LUAD/rna_counts.RDS"
dataset_name <- "LUAD_rna"

rna <- rna_processing(dataset_name, data_path, cnv_filt)
rna_norm <- rna[[1]]
rna_tum <- rna[[2]]

# Gene selection #

# Exclude genes with low expression in normal tissue #
low_expression_threshold <- 15
expression_summary <- data.frame(
  Gene = rownames(rna_norm),
  MeanExpression = rowMeans(rna_norm)
)

summary(expression_summary$MeanExpression)

filtered_genes <- expression_summary %>%
  filter(MeanExpression > low_expression_threshold)
rna_norm <- rna_norm[filtered_genes$Gene, ]
rna_tum <- rna_tum[filtered_genes$Gene, ]
rna <- cbind(rna_norm, rna_tum)

# Expression variability check #
gene_iqr <- apply(rna, 1, IQR)
gene_sd <- apply(rna, 1, sd)

variability_summary <- data.frame(
  Gene = rownames(rna),
  IQR = gene_iqr,
  SD = gene_sd
)

iqr_threshold <- quantile(variability_summary$IQR, 0.75)
#sd_threshold <- quantile(variability_summary$SD, 0.75)

filtered_genes_iqr <- variability_summary %>%
  filter(IQR > iqr_threshold)

#filtered_genes_sd <- variability_summary %>%
#filter(SD > sd_threshold)

rna_filt <- rna[filtered_genes_iqr$Gene, ]


## Correlation analysis ##

# Calculate Pearson correlation between CNV and expression for a gene
rna <- rna %>%  as.matrix()
rna_log_normalized <- DESeq2::vst(rna)
rna_tumor <- rna_log_normalized[,21:40]

common_genes <- intersect(rownames(cnv_filt), rownames(rna_tumor))
cnv_filt <- cnv_filt[common_genes, ]
rna_tumor <- rna_tumor[common_genes, ]

# Identify significant dosage dependent genes
gene_correlations <- sapply(rownames(cnv_filt), function(gene) {
  calculate_correlation(cnv_filt[gene, ], rna_tumor[gene, ])
})

cor_results <- data.frame(
  Gene = rownames(cnv_filt),
  Correlation = gene_correlations
)
rownames(cor_results) <- cor_results$Gene

gene_p_values<- sapply(rownames(cnv_filt), function(gene) {
  calculate_p_value(cnv_filt[gene, ], rna_tumor[gene, ])
})

cor_results$p_value <- gene_p_values
cor_results$adjusted_p_value <- p.adjust(cor_results$p_value, method = "BH")

pvalue_threshold <- 0.05
correlation_threshold <- 0.3

significant_genes <- cor_results%>%
  filter(adjusted_p_value < pvalue_threshold & Correlation > correlation_threshold)


## Identify genes with linear dosage dependance (linear model fit) ##
rna <- rna %>%  as.matrix()
rna_log_normalized <- DESeq2::vst(rna)
rna_tumor <- rna_log_normalized[,21:40]

common_genes <- intersect(rownames(cnv_filt), rownames(rna_tumor))
cnv_filt <- cnv_filt[common_genes, ]
rna_tumor <- rna_tumor[common_genes, ]

model_results <- t(sapply(rownames(cnv_filt), function(gene) {
  fit_linear_model(cnv_filt[gene, ], rna_tumor[gene, ])
}))

results_df <- data.frame(Gene = rownames(model_results), Slope = model_results[, 1], p.value = model_results[, 2])
results_df$adj.p.value <- p.adjust(results_df$p.value, method = "BH")

# Filter for Linear Dosage Dependence #
pvalue_threshold <- 0.05
slope_threshold <- 0.2

dosage_dependent_genes <- results_df %>%
  filter(adj.p.value < pvalue_threshold & abs(Slope) > slope_threshold)

dosage_independent_genes <- results_df %>%
  filter(adj.p.value >= pvalue_threshold | abs(Slope) <= slope_threshold)

# Subset both datasets to common genes #
common_genes <- intersect(rownames(cnv_filt), rownames(dosage_dependent_genes))
cnv_filt <- cnv_filt[common_genes, ]
luad_rna <- luad_rna[common_genes, ]

cnv_normal <- matrix(2, nrow(luad_rna), 20)
rownames(cnv_normal) <- rownames(cnv_filt)
cnv <- cbind(cnv_normal, cnv_filt)
colnames(cnv) <- colnames(rna)
cnv <- cnv/2
colnames(rna) <- paste0("sample", 1:(ncol(rna)))
colnames(cnv) <- colnames(rna)

#Generate metadata#
metadata <- data.frame(patID = colnames(rna),
                       condition = rep(c("A", "B"), each = 20))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 
metadata$condition <- as.factor(metadata$condition)

saveRDS(cnv_filt, file = "TCGA/lung_cancer/LUAD/cnv_filt.RDS")
write.csv(cnv, file = "TCGA/lung_cancer/LUAD/cnv_test_4.csv", row.names = T)
write.csv(rna, file = "TCGA/lung_cancer/LUAD/rna_test_4.csv", row.names = T)
write.csv(metadata, file = "TCGA/lung_cancer/LUAD/metadata_4.csv", row.names = T)


### Explortive analysis RNA z-score vs CNV ###

# Normalization 
rna_log_normalized <- rna %>% as.matrix() %>% DESeq2::rlog()

# Z-score trasformation
rna_zscore <- t(scale(t(rna_log_normalized)))
rna_zscore_normal <- rna_zscore %>% as.data.frame() %>% dplyr::select(1:32)
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% dplyr::select(33:64) 

rna_zscore_normal <- rna_zscore_normal %>% 
  as.data.frame() %>%
  dplyr::mutate(rna_mean = rowMeans(rna_zscore_normal)) %>% 
  dplyr::select(rna_mean)

rna_zscore_tumor <- rna_zscore_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(rna_mean = rowMeans(rna_zscore_tumor)) %>% 
  dplyr::select(rna_mean)

#CNV factorization
cnv <- cnv_mean %>% 
  dplyr::mutate(cnv = case_when(
    cnv_mean <= 0.5 ~ "0",
    cnv_mean > 0.5 & cnv_mean <= 1.8 ~ "1",
    cnv_mean > 1.8 & cnv_mean <= 2.4 ~ "2",
    cnv_mean > 2.4 & cnv_mean <= 3.5 ~ "3",
    cnv_mean > 3.5 & cnv_mean <= 4.4 ~ "4",
    cnv_mean > 4.4 ~ "5")) %>% 
  dplyr::select(cnv)

cnv <- cnv[rownames(cnv) %in% rownames(rna_zscore_normal),]
cnv <- cnv %>% dplyr::select(cnv)

#cnv <- cnv %>% 
  #mutate(cn_group = case_when(
    #cnv == "2" ~ "diploid",
    #cnv == "0" ~ "cn_loss",
    #cnv == "1" ~ "cn_loss",
    #cnv == "3" ~ "cn_gain",
    #cnv == "4" ~ "cn_gain",
    #cnv == "5" ~ "cn_amplification"))

#gene group factorization
#res_noCNV <- res_noCNV %>% 
  #mutate(gene_group = case_when(
    #B1_1 < -0.6 & padj < 0.05 ~ "DEG",
    #B1_1 > 0.6 & padj < 0.05 ~ "DEG",
    #B1_1 <= 0.6 & B1_1 >= -0.6 ~ "noDEG",
    #B1_1 >= 0.6 & padj > 0.05 ~ "not significant",
    #B1_1 <= -0.6 & padj > 0.05 ~ "not significant"))

#Gene group facrorization based on Effect size difference
#deg <- deg %>% 
  #mutate(gene_group = case_when(
    #Difference <= -1.0 ~ "super-dosage",
    #Difference >= 1.0 ~ "super-dosage",
    #Difference > -1.0 & Difference < -0.4 ~ "dosage-sensitive",
    #Difference < 1.0 & Difference > 0.4 ~ "dosage-sensitive",
    #Difference >= -0.4 & Difference <= 0.4 ~ "dosage_insensitive"
  #))

#Calculate SD across rows and columns 
#library(matrixStats)

#rna <- transform(rna, sd=apply(rna, 1, sd, na.rm=TRUE)) %>% 
  #subset(sd > 25.0) %>% 
  #select(1:192) %>% 
  #as.matrix()

#cnv_tumor <- as.matrix(cnv_tum)
#cnv_tumor_sd <- matrixStats::colSds(cnv_tumor) %>% 
  #as.data.frame() %>% 
  #setNames("sd")  %>% 
  #subset(sd > 1)