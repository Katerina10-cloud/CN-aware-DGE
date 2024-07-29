# Identify homogeneous groups of patients based on their CNV profiles #

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/")
pkgs <- c("dplyr", "ggplot2", "cluster", "factoextra", "heatmaply", "DESeq2")
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
kmeans_result <- kmeans(pca_data, centers = 5)

# Add cluster assignments to the original data
cnv_tumor <- as.data.frame(cnv_tumor)
cnv_tumor$Cluster <- kmeans_result$cluster

silhouette_result <- silhouette(kmeans_result$cluster, dist(pca_data))
avg_silhouette_width <- mean(silhouette_result[, "sil_width"])
print(paste("Average Silhouette Width:", avg_silhouette_width))

# Visualization
fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
  ggtitle("PCA Clustering of CNV Profiles")

cnv_filt <- cnv_tumor[cnv_tumor$Cluster %in% c(3,4,5),]

cnv_filt <- as.matrix(t(cnv_filt))


cnv_filt <- apply(cnv_filt, 2, function(x) ifelse(x > 10, 10, x)) 

hist(rowMeans(cnv_filt),
     main = "LUAD", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 6)

cnv_mean <- cnv_filt %>% 
  as.data.frame() %>% 
  mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  select(cnv_mean)

# RNA data processing #
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

luad_rna <- cbind(luad_norm, luad_tum)

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

save(luad_rna, file = "lung_cancer/LUAD/rna_test.Rdata")
save(cnv_filt, file = "lung_cancer/LUAD/cnv_filt.Rdata")
