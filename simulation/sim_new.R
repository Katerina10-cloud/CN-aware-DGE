# Simulate differential gene expression induced by copy number variation #

library(MASS)
library(ggplot2)

set.seed(123)

num_genes <- 1000
num_samples <- 20
mu <- rlnorm(num_genes, meanlog = 10, sdlog = 6)
dispersion <- 0.1

# Simulate gene expression count data
#rna_counts <- matrix(0, nrow = num_genes, ncol = num_samples)
#for (i in 1:num_genes) {
  #rna_counts[i, ] <- rnbinom(num_samples, size = size, mu = mu)
#}

rna_counts <- matrix(rnbinom(num_genes * num_samples, size = 1/dispersion, mu = mu), 
                          nrow = num_genes, ncol = num_samples)
colnames(gene_expression) <- paste0("Sample_", 1:n_samples)
rownames(gene_expression) <- paste0("Gene_", 1:n_genes)

#for (i in 1:num_genes) {
  #r <- mean_expression[i]   # Shape parameter for Negative Binomial
  #p <- 1 / (1 + mean_expression[i])  # Probability parameter for Negative Binomial
  #rna_counts[i, ] <- rnbinom(num_samples, size = r, prob = p)
#}

#rna_counts <- rna_counts + 1


# Simulate CN data #

#cnv_data <- matrix(sample(c(1, 2, 3, 4, 5), num_genes * num_samples, replace = TRUE), 
                   #nrow = num_genes, ncol = num_samples)

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

# Apply CNV effect to RNA-seq count data
# For simplicity, assume that CNV affects counts multiplicatively
for (i in 1:num_genes) {
  for (j in 1:num_samples) {
    rna_counts[i, j] <- rna_counts[i, j] * (cnv[i, j] /2)
  }
}

# Add row and column names
rna_counts <- as.data.frame(rna_counts)
cnv <- as.data.frame(cnv)
rownames(rna_counts) <- paste0("Gene_", 1:num_genes)
colnames(rna_counts) <- paste0("Sample_", 1:num_samples)
rownames(cnv) <- paste0("Gene_", 1:num_genes)
colnames(cnv) <- paste0("Sample_", 1:num_samples)

# Generate metadata #
metadata <- data.frame(patID = colnames(rna_counts),
                       condition = rep(c("0", "1"), each = 10))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 

