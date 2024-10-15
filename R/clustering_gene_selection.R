### Identify homogeneous groups of patients based on their CNV profiles ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "DESeq2", "colorspace", 
          "ggpubr", "ggpointdensity", "ggeasy")
sapply(pkgs, require, character.only = TRUE)

source("CN-aware-DGE/R/utils.R")
#set.seed(12345)

# Input data
data_path <- "TCGA/lung/LUSC/cnv_tumor.RDS"
dataset_name <- "LUSC_cnv"
cnv_tumor <- readRDS(data_path)

# Clustering patients 
clustering_res <- clustering_patients_cnv(dataset_name, data_path)
clustering_res[[2]]
cnv_tumor <- clustering_res[[1]]

# Cluster selection and analysis
cnv_filt <- cnv_tumor[cnv_tumor$Cluster %in% c(1,2,3),]
cnv_filt <- subset(cnv_filt, select=-c(Cluster))
cnv_filt <- as.matrix(t(cnv_filt))
cnv_filt <- apply(cnv_filt, 2, function(x) ifelse(x > 10, 10, x))

cnv_tumor <- apply(cnv_tumor, 2, function(x) ifelse(x > 10, 10, x))

#saveRDS(cnv_filt, file = "TCGA/liver/test/cnv_filt_2.RDS")

hist(rowMeans(cnv_filt),
     main = "LUSC", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 15)

cnv_mean <- cnv_filt %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  dplyr::select(cnv_mean)

cnv_mean$geneID <- rownames(cnv_mean)


## RNA data processing ##

data_path <- "TCGA/colon/rna_counts.RDS"
dataset_name <- "COLON_rna"

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

# Filter DEGs
pval_cut <- 0.05
lfc_cut <- 1.0

res_deg <- res_naive_edge %>% 
  dplyr::filter(FDR < pval_cut, abs(logFC) > lfc_cut) %>% 
  dplyr::mutate(isDE = (abs(logFC) >= lfc_cut) & (FDR <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(logFC > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::select(logFC, FDR, isDE, DEtype) 

volcano <-  res_deg %>% 
  ggplot(mapping = aes(x=logFC, y=-log10(FDR), col=DEtype)) +
  geom_point(size=.8) +
  theme_bw() +
  scale_color_manual(values = de_gene_colors) +
  #ggh4x::facet_nested(tool~method, scale="free") +
  #facet_wrap(~method, nrow=1, scale = "free")+
  ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col="") +
  ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  ggplot2::geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'bottom')
volcano


# Calculate Pearson correlation between CNV and expression for a gene
rna_filt <- rna[ rownames(rna) %in% rownames(res_deg) , ]
rna_log_normalized <- rna_filt %>%  as.matrix() %>% DESeq2::vst()
rna_tumor <- rna_log_normalized[,21:40]

#common_genes <- intersect(rownames(cnv_filt), rownames(rna))
#cnv_filt <- cnv_filt[common_genes, ]
#rna <- rna[common_genes, ]

gene_correlations <- sapply(rownames(cnv_filt), function(gene) {
  calculate_correlation(cnv_filt[gene, ], rna_tumor[gene, ])
})

cor_results <- data.frame(
  Gene = rownames(cnv_mean),
  Correlation = gene_correlations
)

rownames(cor_results) <- cor_results$Gene

ggplot(data = cor_results, aes(x = Correlation)) +
  geom_histogram(color = "white", fill = "lightblue", bins = 30)+
  ggplot2::labs(x = "Correlation coefficient", y = "DEG count", col="")+
  ggplot2::geom_vline(xintercept = 0.4, linetype = 'dashed', colour="darkred")+
  theme_bw() 
  

#gene_p_values<- sapply(rownames(cnv_filt), function(gene) {
  #calculate_p_value(cnv_filt[gene, ], rna_tumor[gene, ])
#})

#cor_results$p_value <- gene_p_values
#cor_results$adjusted_p_value <- p.adjust(cor_results$p_value, method = "BH")

#pvalue_threshold <- 0.05
#correlation_threshold <- 0.3

#significant_genes <- cor_results%>%
  #filter(adjusted_p_value < pvalue_threshold & Correlation > correlation_threshold)


## Identify genes with linear dosage dependance (linear model fit) ##
rna_log_normalized <- rna %>% as.matrix %>% DESeq2::vst()
rna_tumor <- rna_log_normalized[,12:22]

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
slope_threshold <- 0.3

dosage_dependent_genes <- results_df %>%
  filter(adj.p.value < pvalue_threshold & abs(Slope) > slope_threshold)

dosage_independent_genes <- results_df %>%
  filter(adj.p.value >= pvalue_threshold | abs(Slope) <= slope_threshold)

# Subset both datasets to common genes #
common_genes <- intersect(rownames(cnv_tumor), rownames(rna_filt))
cnv_filt <- cnv_tumor[common_genes, ]
rna_filt <- rna_filt[common_genes, ]

cnv_normal <- matrix(2, nrow(rna_filt), 12)
rownames(cnv_normal) <- rownames(cnv_filt)
cnv <- cbind(cnv_normal, cnv_filt)
colnames(cnv) <- colnames(rna_filt)
cnv <- cnv/2
colnames(rna_filt) <- paste0("sample", 1:(ncol(rna_filt)))
colnames(cnv) <- colnames(rna_filt)

#Generate metadata#
metadata <- data.frame(patID = colnames(rna_filt),
                       condition = rep(c("A", "B"), each = 12))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 
metadata$condition <- as.factor(metadata$condition)

#saveRDS(cnv_filt, file = "TCGA/lung_cancer/LUAD/cnv_filt.RDS")
write.csv(cnv, file = "TCGA/colon/test/cnv_test.csv", row.names = T)
write.csv(rna_filt, file = "TCGA/colon/test/rna_test.csv", row.names = T)
write.csv(metadata, file = "TCGA/colon/test/metadata.csv", row.names = T)


