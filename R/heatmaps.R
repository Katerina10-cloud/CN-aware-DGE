setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "pheatmap", "grid")
sapply(pkgs, require, character.only = TRUE)

### Heatmap RNA ###

set.seed(123)  
n_genes <- 10  
n_tumor_samples <- 5  
n_normal_samples <- 5  
n_total_samples <- n_tumor_samples + n_normal_samples

normal_counts <- matrix(rpois(n_genes * n_normal_samples, lambda = 50), 
                        nrow = n_genes, ncol = n_normal_samples)

tumor_counts <- matrix(rpois(n_genes * n_tumor_samples, lambda = 50), 
                       nrow = n_genes, ncol = n_tumor_samples)

fold_change <- 2  
tumor_counts[1:(n_genes / 2), ] <- tumor_counts[1:(n_genes / 2), ] * fold_change
rna_counts <- cbind(normal_counts, tumor_counts)

rownames(rna_counts) <- paste("gene", 1:n_genes, sep = "_")
colnames(rna_counts) <- c(paste("normal", 1:n_normal_samples, sep = "_"), 
                          paste("tumor", 1:n_tumor_samples, sep = "_"))

#rna_counts_log <- log2(rna_counts + 1)
sample_annotation <- data.frame(Condition = rep(c("Normal", "Tumor"), 
                                             times = c(n_normal_samples, n_tumor_samples)))
rownames(sample_annotation) <- colnames(rna_counts)


heatmap_colors_rna <- colorRampPalette(c("white", "#889043"))(50)
status_colors <- c("Normal" = "#CCCCCC", "Tumor" = "#FF9896")
annotation_colors <- list(Condition = status_colors)

pheatmap(rna_counts, 
         scale = "row",  
         color = heatmap_colors_rna,  
         annotation_col = sample_annotation,  
         annotation_colors = annotation_colors,  
         cluster_rows = F,  
         cluster_cols = F, 
         show_rownames = F,
         show_colnames = F,
         main = "")


### Heatmap CN ###

rm(list = ls())
set.seed(456)  

n_genes <- 10  
n_samples <- 5

normal_data<- matrix(rnorm(n_genes * n_samples, mean = 2, sd = 0.1), 
                      nrow = n_genes, ncol = n_samples)

tumor_data <- normal_data + rnorm(n_genes * n_samples, mean = 0, sd = 1.5)

rownames(normal_data) <- rownames(tumor_data) <- paste("gene", 1:n_genes, sep = "_")
colnames(normal_data) <- paste("normal", 1:n_samples, sep = "_")
colnames(tumor_data) <- paste("tumor", 1:n_samples, sep = "_")

cn_data <- cbind(normal_data, tumor_data)

sample_annotation <- data.frame(Condition = rep(c("Normal", "Tumor"), 
                                                times = c(5, 5)))
rownames(sample_annotation) <- colnames(cn_data)

status_colors <- c("Normal" = "#CCCCCC", "Tumor" = "#FF9896")
annotation_colors <- list(Condition = status_colors)
cn_data <- abs(cn_data)

pheatmap(cn_data, 
         scale = "none",
         color = colorRampPalette(c("#3C5488B2", "white", "#E64B35B2"))(50),
         annotation_col = sample_annotation,  
         annotation_colors = annotation_colors,
         cluster_rows = F,  
         cluster_cols = F,  
         show_rownames = F,
         show_colnames = F,
         main = "")


# Heatmap simulation setting #

n_genes <- 10  
n_samples <- 10  
tumor_samples <- 5
healthy_samples <- 5

# Define gene groups
dosage_sensitive_genes <- 2   
dosage_insensitive_genes <- 4 
non_DE_genes <- n_genes - dosage_sensitive_genes - dosage_insensitive_genes # Remaining are non-DE

# Simulate baseline expression
baseline_expression <- rnorm(n_genes, mean = 5, sd = 1)  # Baseline expression for each gene

# Simulate CNV effects and differential expression factors
CNV_factors_matrix <- matrix(runif(dosage_sensitive_genes * tumor_samples, min = 0.5, max = 3), 
                             nrow = dosage_sensitive_genes, ncol = tumor_samples)

#diff_expr_factor <- 1.5  # Constant differential expression factor for dosage-insensitive genes

diff_expr_factors_matrix <- matrix(runif(dosage_insensitive_genes * tumor_samples, min = 1.5, max = 2), 
                                   nrow = dosage_insensitive_genes, ncol = tumor_samples)
# Create expression matrix
expression_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
colnames(expression_matrix) <- c(paste0("Healthy_", 1:healthy_samples), paste0("Tumor_", 1:tumor_samples))
rownames(expression_matrix) <- paste0("Gene_", 1:n_genes)

# Assign baseline expression to healthy samples
for (i in 1:healthy_samples) {
  expression_matrix[, i] <- baseline_expression
}

# Modify expression for tumor samples
# Dosage-sensitive genes: baseline * CNV_factors_matrix (unique per sample)
for (i in 1:tumor_samples) {
  expression_matrix[1:dosage_sensitive_genes, healthy_samples + i] <- 
    baseline_expression[1:dosage_sensitive_genes] * CNV_factors_matrix[, i]
}

# Dosage-insensitive genes: baseline * diff_expr_factors_matrix (unique per sample)
for (i in 1:tumor_samples) {
  expression_matrix[(dosage_sensitive_genes + 1):(dosage_sensitive_genes + dosage_insensitive_genes), healthy_samples + i] <- 
    baseline_expression[(dosage_sensitive_genes + 1):(dosage_sensitive_genes + dosage_insensitive_genes)] * diff_expr_factors_matrix[, i]
}

# Non-differentially expressed (non-DE) genes: baseline expression
for (i in 1:tumor_samples) {
  expression_matrix[(dosage_sensitive_genes + dosage_insensitive_genes + 1):n_genes, healthy_samples + i] <- 
    baseline_expression[(dosage_sensitive_genes + dosage_insensitive_genes + 1):n_genes]
}

# Scale data for better heatmap visualization
scaled_matrix <- t(scale(t(expression_matrix)))


# Replace NaN values resulting from zero variance with zero
scaled_matrix[is.na(scaled_matrix)] <- 0

# Create a heatmap
#annotation <- data.frame(Sample_Type = rep(c("Healthy", "Tumor"), each = 5))
#rownames(annotation) <- colnames(scaled_matrix)

sample_annotation <- data.frame(Condition = rep(c("Healthy", "Tumor"), 
                                                times = c(healthy_samples, tumor_samples)))
rownames(sample_annotation) <- colnames(scaled_matrix)


heatmap_colors_rna <- colorRampPalette(c("white", "#889043"))(50)
status_colors <- c("Healthy" = "#CCCCCC", "Tumor" = "#FF9896")
annotation_colors <- list(Condition = status_colors)

max_expression <- max(scaled_matrix, na.rm = TRUE)
breaks <- seq(0, max_expression, length.out = 51) 

pheatmap(scaled_matrix, 
         #scale = "row",
         cluster_rows = F, 
         cluster_cols = F,
         annotation_col = sample_annotation,
         annotation_colors = annotation_colors,
         show_rownames = F,
         show_colnames = F,
         color = heatmap_colors_rna,
         breaks = breaks,
         border_color = "black",
         main = "")

