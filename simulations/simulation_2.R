setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "DESeq2", "ggplot2", "gridExtra", "ggpubr", "ggrepel", "ggvenn", "ggpointdensity")
sapply(pkgs, require, character.only = TRUE)

# RNA counts simulation (DESeq2) #
n_genes = 1000
n-samples = 20
sizeFactors = rep(1, m)

deseq_sim <- DESeq2::makeExampleDESeqDataSet(
  n = n_genes,
  m = n_samples,
  betaSD = 0.2,
  interceptMean = 8,
  interceptSD = 4,
  dispMeanRel = function(x) 8/x + 0.3,
  sizeFactors = sizeFactors
)


rna_counts <- data.frame(deseq_sim@assays@data@listData[["counts"]])

# Generate metadata #
metadata <- data.frame(patID = colnames(rna_counts),
                       condition = rep(c("A", "B"), each = 10))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 

# CNV simulation #
cnv_0 <- sapply(1:10, function(x) sample(x=c(0.5,1,2,3), size = 200, replace=TRUE, prob = c(.50, .30, .10, .10)))
cnv_1 <- sapply(1:10, function(x) sample(x=c(0.5,1,2,3), size = 200, replace=TRUE, prob = c(.10, .70, .10, .10)))
cnv_2 <- sapply(1:10, function(x) sample(x=c(1,2,3,4), size = 200, replace=TRUE, prob = c(.10, .70, .10, .10)))
cnv_3 <- sapply(1:10, function(x) sample(x=c(1,2,3,4), size = 100, replace=TRUE, prob = c(.05, .50, .80, .10)))
cnv_4 <- sapply(1:10, function(x) sample(x=c(2,3,4,5), size = 100, replace=TRUE, prob = c(.05, .50, .80, .10)))
cnv_5 <- sapply(1:10, function(x) sample(x=c(2,3,4,5), size = 200, replace=TRUE, prob = c(.05, .50, .10, .80)))
cnv_tumor <- rbind(cnv_0, cnv_1, cnv_2, cnv_3, cnv_4, cnv_5) %>% as.matrix()
cnv_normal <- matrix(2, nrow(rna_counts), 10)
cnv <- cbind(cnv_normal, cnv_tumor)

cnv <- cnv/2
#cnv <- apply(cnv, 2, function(x) x/2)

colnames(cnv) <- colnames(rna_counts)
rownames(cnv) <- paste0("G", 1:(nrow(cnv)))
rownames(rna_counts) <- paste0("G", 1:(nrow(rna_counts)))

# Simulate DGE induced by CN #
rna_cnv <- rna_counts * cnv
rna_cnv <- ceiling(rna_cnv)

write.csv(rna_cnv, file = "CN-aware-DGE/Python/datasets/rna_counts_cnv_v2.csv", row.names = T)

