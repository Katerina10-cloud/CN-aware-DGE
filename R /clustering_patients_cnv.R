# Identify homogeneous groups of patients based on their CNV profiles #

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/")
pkgs <- c("dplyr", "ggplot2", "cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "DESeq2")
sapply(pkgs, require, character.only = TRUE)

load("~/Documents/PhD_AI/TCGA/lung_cancer/LUAD/cnv_tumor.Rdata")

cnv_tumor <- as.matrix(t(luad_cnv_tumor)) 

# Normalization
cnv_data_normalized <- scale(cnv_tumor)

# Using PCA for dimensionality reduction
pca_result <- prcomp(cnv_data_normalized, scale. = TRUE)
pca_data <- as.data.frame(pca_result$x[, 1:10])

# Clustering
set.seed(123)  
kmeans_result <- kmeans(pca_data, centers = 3)

# Add cluster assignments to the original data
cnv_tumor <- as.data.frame(cnv_tumor)
cnv_tumor$Cluster <- kmeans_result$cluster

#silhouette_result <- silhouette(kmeans_result$cluster, dist(pca_data))
#avg_silhouette_width <- mean(silhouette_result[, "sil_width"])
#print(paste("Average Silhouette Width:", avg_silhouette_width))

# Visualization
fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
  ggtitle("PCA of CNV Profiles")+
  theme_classic()

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

# RNA data processing #
load("~/Documents/PhD_AI/TCGA/lung_cancer/LUAD/rna_counts.Rdata")

colnames(luad_rna_sd) <- gsub(pattern = "\\.", replacement = "-", colnames(luad_rna_sd))
luad_tum <- luad_rna_sd %>% select(46:90)
luad_norm <- luad_rna_sd %>% select(1:45)

colnames(luad_norm) <- stringr::str_sub(colnames(luad_norm),1,12)
colnames(luad_tum) <- stringr::str_sub(colnames(luad_tum),1,12)

cnv_filt <- as.data.frame(cnv_filt)
luad_tum <- luad_tum[,colnames(luad_tum) %in% colnames(cnv_filt)]
luad_norm <- luad_norm[,colnames(luad_norm) %in% colnames(cnv_filt)]

x <- colnames(luad_norm)
names(luad_norm) <- paste(x,"-11A")

x <- colnames(luad_tum)
names(luad_tum) <- paste(x,"-01A")

# Gene selection #

# Exclude genes with low expression in normal tissue #
low_expression_threshold <- 15
expression_summary <- data.frame(
  Gene = rownames(luad_norm),
  MeanExpression = rowMeans(luad_norm)
)

summary(expression_summary$MeanExpression)

filtered_genes <- expression_summary %>%
  filter(MeanExpression > low_expression_threshold)
luad_norm <- luad_norm[filtered_genes$Gene, ]
luad_tum <- luad_tum[filtered_genes$Gene, ]
luad_rna <- cbind(luad_norm, luad_tum)

# Expression variability check #
gene_iqr <- apply(luad_rna, 1, IQR)
gene_sd <- apply(luad_rna, 1, sd)

variability_summary <- data.frame(
  Gene = rownames(luad_rna),
  IQR = gene_iqr,
  SD = gene_sd
)

iqr_threshold <- quantile(variability_summary$IQR, 0.75)
#sd_threshold <- quantile(variability_summary$SD, 0.75)

filtered_genes_iqr <- variability_summary %>%
  filter(IQR > iqr_threshold)

#filtered_genes_sd <- variability_summary %>%
  #filter(SD > sd_threshold)

luad_rna_filt <- luad_rna[filtered_genes_iqr$Gene, ]


# Correlation analysis #

# Calculate Pearson correlation between CNV and expression for a gene
luad_rna <- luad_rna %>%  as.matrix()
rna_log_normalized <- DESeq2::vst(luad_rna)
rna_tumor <- rna_log_normalized[,21:40]

common_genes <- intersect(rownames(cnv_filt), rownames(rna_tumor))
cnv_filt <- cnv_filt[common_genes, ]
rna_tumor <- rna_tumor[common_genes, ]

calculate_correlation <- function(cnv_values, expression_values) {
  cnv_values <- as.numeric(as.character(cnv_values))
  expression_values <- as.numeric(as.character(expression_values))
  valid_indices <- !is.na(cnv_values) & !is.na(expression_values)
  # Check if there are enough valid values to calculate correlation
  if (sum(valid_indices) > 1) {
    return(cor.test(cnv_values[valid_indices], expression_values[valid_indices], method = "pearson")$estimate)
  } else {
    return(NA)  # Return NA if not enough valid data
  }
}

calculate_p_value <- function(cnv_values, expression_values) {
  cnv_values <- as.numeric(as.character(cnv_values))
  expression_values <- as.numeric(as.character(expression_values))
  valid_indices <- !is.na(cnv_values) & !is.na(expression_values)
  # Check if there are enough valid values to calculate correlation
  if (sum(valid_indices) > 1) {
    return(cor.test(cnv_values[valid_indices], expression_values[valid_indices], method = "pearson")$p.value)
  } else {
    return(NA)  # Return NA if not enough valid data
  }
  
}

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



# Data generation and processing #

# Subset both datasets to common genes #
common_genes <- intersect(rownames(cnv_filt), rownames(significant_genes))
cnv_filt <- cnv_filt[common_genes, ]
luad_rna <- luad_rna[common_genes, ]

cnv_normal <- matrix(2, nrow(luad_rna), 20)
rownames(cnv_normal) <- rownames(cnv_filt)
cnv <- cbind(cnv_normal, cnv_filt)
colnames(cnv) <- colnames(luad_rna)
cnv <- cnv/2
colnames(luad_rna) <- paste0("sample", 1:(ncol(luad_rna)))
colnames(cnv) <- colnames(luad_rna)

#Generate metadata#
metadata <- data.frame(patID = colnames(luad_rna),
                       condition = rep(c("A", "B"), each = 20))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 
metadata$condition <- as.factor(metadata$condition)


write.csv(cnv, file = "lung_cancer/LUAD/cnv_test_3.csv", row.names = T)
write.csv(luad_rna, file = "lung_cancer/LUAD/rna_test_3.csv", row.names = T)
write.csv(metadata, file = "lung_cancer/LUAD/metadata_3.csv", row.names = T)


# Scatter plot #
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/CN-aware-DGE/Python/")

res <- read.csv("results/res_CNnaive_test3.csv")
res_adj <- read.csv("results/res_CNaware_test3.csv")


res <- res %>% select(X,log2FoldChange)
res <- res %>% remove_rownames %>% column_to_rownames(var="X")
res_adj <- res_adj %>% select(X,log2FoldChange)
res_adj <- res_adj %>% remove_rownames %>% column_to_rownames(var="X")

common_genes <- intersect(rownames(cnv_mean), rownames(res))
cnv_mean <- cnv_mean[common_genes, ]
cnv_mean <- data.frame(cnv_mean)

plot_data <- merge(cnv_mean, res, by = "row.names")
plot_data_adj <- merge(cnv_mean, res_adj, by = "row.names")

colnames(plot_data) <- c("geneID", "cnv_mean", "logFC")
colnames(plot_data_adj) <- c("geneID", "cnv_mean", "logFC")

plot_data <- plot_data %>% dplyr::filter(logFC > -2.0,)
plot_data_adj <- plot_data_adj %>% dplyr::filter(logFC > -4.0,)

scatter1 = ggplot(plot_data, aes(x = log2(cnv_mean), y = logFC)) +
  geom_point(size=0.8) +
  #geom_smooth(method = "lm", color = "blue", se = T)+
  geom_smooth()+
  theme(legend.position = 'bottom') +
  labs(title = "Cluster 2 - no CN normalization", x ="log2(CN ratio)", y="log2FC") +
  theme_classic()
scatter1

scatter2 = ggplot(plot_data_adj, aes(x = log2(cnv_mean), y = logFC)) +
  geom_point(size=0.8) +
  #geom_smooth(method = "lm", color = "blue", se = T)+
  geom_smooth()+
  theme(legend.position = 'bottom') +
  labs(title = "Cluster 2 - with CN normalization", x ="log2(CN ratio)", y="log2FC") +
  theme_classic()
scatter2

gridExtra::grid.arrange(scatter1, scatter2, nrow = 1)



# Explortive analysis RNA z-score vs CNV #

# Normalization #
rna_normalized <- luad_rna %>%  as.matrix()
rna_log_normalized <- DESeq2::rlog(rna_normalized)

# Z-score trasformation
rna_zscore <- t(scale(t(rna_log_normalized)))
rna_zscore_normal <- rna_zscore %>% as.data.frame() %>% select(1:32)
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% select(33:64) 
rna_zscore_normal <- rna_zscore_normal %>% as.data.frame() %>%
  mutate(rna_mean = rowMeans(rna_zscore_normal)) %>% select(rna_mean)
rna_zscore_tumor <- rna_zscore_tumor %>% as.data.frame() %>% 
  mutate(rna_mean = rowMeans(rna_zscore_tumor)) %>% select(rna_mean)

#luad_cnv <- luad_cnv %>% as.data.frame() %>% 
#mutate(cnv_mean = rowMeans(luad_cnv))

#CNV factorization
cnv <- cnv_mean %>% 
  mutate(cnv = case_when(
    cnv_mean <= 0.5 ~ "0",
    cnv_mean > 0.5 & cnv_mean <= 1.8 ~ "1",
    cnv_mean > 1.8 & cnv_mean <= 2.4 ~ "2",
    cnv_mean > 2.4 & cnv_mean <= 3.5 ~ "3",
    cnv_mean > 3.5 & cnv_mean <= 4.4 ~ "4",
    cnv_mean > 4.4 ~ "5")) 
select(cnv)

cnv <- cnv[rownames(cnv) %in% rownames(rna_zscore_normal),]
cnv <- cnv %>% select(cnv)
