setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "DESeq2", "ggplot2", "gridExtra", "ggpubr", "ggrepel", "ggvenn", "ggpointdensity")
sapply(pkgs, require, character.only = TRUE)

# RNA counts simulation (DESeq2) #
n_genes = 1000
n_samples = 20
sizeFactors = rep(1, n_samples)

deseq_sim <- DESeq2::makeExampleDESeqDataSet(
  n = n_genes,
  m = n_samples,
  betaSD = 0.1,
  interceptMean = 8,
  interceptSD = 4,
  dispMeanRel = function(x) 8/x + 0.3,
  sizeFactors = sizeFactors
)

rna_counts <- data.frame(deseq_sim@assays@data@listData[["counts"]])

# Generate metadata #
metadata <- data.frame(patID = colnames(rna_counts),
                       condition = rep(c("A", "B"), each = 25))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 

# CNV simulation #

# Function to generate CNV data for different groups
generate_cnv_data <- function(n_cols, size, values, prob) {
  sapply(1:n_cols, function(x) sample(x = values, size = size, replace = TRUE, prob = prob))
}

# Generate CNV data for each group with varying parameters
cnv_0 <- generate_cnv_data(n_cols = 25, size = 1000, values = c(0.5, 1, 2, 3), prob = c(0.50, 0.30, 0.10, 0.10))
cnv_1 <- generate_cnv_data(n_cols = 25, size = 1000, values = c(0.5, 1, 2, 3), prob = c(0.10, 0.70, 0.10, 0.10))
cnv_2 <- generate_cnv_data(n_cols = 25, size = 1000, values = c(1, 2, 3, 4), prob = c(0.10, 0.70, 0.10, 0.10))
cnv_3 <- generate_cnv_data(n_cols = 25, size = 500, values = c(1, 2, 3, 4), prob = c(0.05, 0.50, 0.80, 0.10))
cnv_4 <- generate_cnv_data(n_cols = 25, size = 750, values = c(2, 3, 4, 5), prob = c(0.05, 0.50, 0.80, 0.10))
cnv_5 <- generate_cnv_data(n_cols = 25, size = 750, values = c(2, 3, 4, 5), prob = c(0.05, 0.50, 0.10, 0.80))

cnv_tumor <- rbind(cnv_0, cnv_1, cnv_2, cnv_3, cnv_4, cnv_5) %>% as.matrix()

cnv_normal <- matrix(2, nrow(rna_counts), 25)
cnv <- cbind(cnv_normal, cnv_tumor)

cnv <- cnv/2
#cnv <- apply(cnv, 2, function(x) x/2)

colnames(cnv) <- colnames(rna_counts)
rownames(cnv) <- paste0("G", 1:(nrow(cnv)))
rownames(rna_counts) <- paste0("G", 1:(nrow(rna_counts)))

# Simulate DGE induced by CN #
rna_cnv <- rna_counts * cnv
rna_cnv <- ceiling(rna_cnv)
rna_cnv <- rna_cnv + 10

write.csv(rna_cnv, file = "CN-aware-DGE/Python/datasets/rna_counts_cnv.csv", row.names = T)
write.csv(cnv, file = "CN-aware-DGE/Python/datasets/cnv.csv", row.names = T)
write.csv(metadata, file = "CN-aware-DGE/Python/datasets/metadata.csv", row.names = T)
