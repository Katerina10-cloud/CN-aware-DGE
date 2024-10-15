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


heatmap_colors_rna <- colorRampPalette(c("white", "#889043"))
status_colors <- c("Normal" = "#CCCCCC", "Tumor" = "#FF9896")
annotation_colors <- list(Condition = status_colors)

pheatmap(rna_counts, 
         scale = "row",  
         color = heatmap_colors,  
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


