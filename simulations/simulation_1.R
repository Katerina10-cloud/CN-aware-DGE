### Simulate RNAseq counts ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "compcodeR", "DESeq2")
sapply(pkgs, require, character.only = TRUE)

#rna <- readRDS("TCGA/lung_cancer/LUAD/rna_counts.RDS")

rna <- read.csv("TCGA/lung_cancer/LUAD/rna_test_3.csv")
metadata <- read.csv("TCGA/lung_cancer/LUAD/metadata_3.csv")
rna <- rna %>% remove_rownames %>% column_to_rownames(var="X")
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="X")

dds <- DESeq2::DESeqDataSetFromMatrix(countData=rna, 
                              colData=metadata, 
                              design=~condition)
dds <- DESeq2::DESeq(dds)

dispersion <- environment(dds@dispersionFunction)[["fit"]][["model"]][["disps[good]"]]
mean <- environment(dds@dispersionFunction)[["fit"]][["data"]][["means"]]

dispersion <- dispersion[190:3189]
mean <- mean[190:3189]


## Using compcodeR simulator ##
rna_counts_sim <- compcodeR::generateSyntheticData(dataset = "rna_counts_sim", n.vars = 3000, 
                                                   samples.per.cond = 20, n.diffexp = 0, 
                                                   repl.id = 1, seqdepth = 1e7, 
                                                   fraction.upregulated = 0.0, 
                                                   between.group.diffdisp = F, 
                                                   filter.threshold.total = 1, 
                                                   filter.threshold.mediancpm = 0, 
                                                   fraction.non.overdispersed = 0, 
                                                   relmeans = mean,
                                                   dispersions = dispersion,
                                                   random.outlier.high.prob = 0,
                                                   random.outlier.low.prob = 0,
                                                   output.file = "rna_counts_sim.rds")

rna_counts_sim <- readRDS("TCGA/lung_cancer/rna_counts_sim.rds")

rna_counts_sim <- rna_counts_sim@count.matrix
rna_normal <- rna_counts_sim %>% as.data.frame() %>% select(1:20)
rna_tumor <- rna_counts_sim %>% as.data.frame() %>% select(21:40)
rna_counts <- cbind(rna_normal, rna_tumor)


# Generate metadata
metadata <- data.frame(patID = colnames(rna_counts),
                       condition = rep(c("A", "B"), each = 20))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 



### Simulate CN data ###

# Function to generate CNV data for different groups
generate_cnv_data <- function(n_cols, size, values, prob) {
  sapply(1:n_cols, function(x) sample(x = values, size = size, replace = TRUE, prob = prob))
}

# Generate CNV data for each group with varying parameters
cnv_0 <- generate_cnv_data(n_cols = 20, size = 400, values = c(0.5, 1, 2, 3), prob = c(0.50, 0.30, 0.10, 0.10))
cnv_1 <- generate_cnv_data(n_cols = 20, size = 400, values = c(0.5, 1, 2, 3), prob = c(0.10, 0.70, 0.10, 0.10))
cnv_2 <- generate_cnv_data(n_cols = 20, size = 800, values = c(1, 2, 3, 4), prob = c(0.10, 0.70, 0.10, 0.10))
cnv_3 <- generate_cnv_data(n_cols = 20, size = 400, values = c(1, 2, 3, 4), prob = c(0.05, 0.05, 0.80, 0.10))
cnv_4 <- generate_cnv_data(n_cols = 20, size = 500, values = c(2, 3, 4, 5), prob = c(0.05, 0.05, 0.80, 0.10))
cnv_5 <- generate_cnv_data(n_cols = 20, size = 500, values = c(2, 3, 4, 5), prob = c(0.05, 0.05, 0.10, 0.80))

cnv_tumor <- rbind(cnv_0, cnv_1, cnv_2, cnv_3, cnv_4, cnv_5) %>% as.matrix()

cnv_normal <- matrix(2, nrow(rna_counts), 20)
cnv <- cbind(cnv_normal, cnv_tumor)

cnv <- cnv/2

colnames(cnv) <- colnames(rna_counts)
rownames(cnv) <- rownames(rna_counts)

rna_cnv <- rna_counts * cnv
rna_cnv <- ceiling(rna_cnv)


write.csv(rna_cnv, file = "CN-aware-DGE/Python/datasets/rna_counts_cnv_v2.csv", row.names = T)
write.csv(cnv, file = "CN-aware-DGE/Python/datasets/cnv_v2.csv", row.names = T)
write.csv(metadata, file = "CN-aware-DGE/Python/datasets/metadata_v2.csv", row.names = T)

## Prepare CN data ##
#cnv <- apply(cnv, 1, function(x) x/2)
#cnv <- apply(cnv, 1, function(x) x+10e-9)


## Simulated CN data preprocessing OMICSSimLA ##

#cnv_sim_luad <- read.table(file = "data/GbmSim1.tsv", header = FALSE, sep = "")
#cnv_sim_luad <- cnv_sim_luad[-c(2:3)]
#cnv_sim_luad <- cnv_sim_luad[-which(duplicated(cnv_sim_luad$V1)), ]
#cnv_sim_luad <- cnv_sim_luad %>% remove_rownames %>% column_to_rownames(var="V1")
#colnames(cnv_sim_luad) <- paste0("G", 1:(ncol(cnv_sim_luad)-1))
#cnv_sim_luad <- t(cnv_sim_luad)

# Assigning CN states #
#cnv_sim_luad <- apply(cnv_sim_luad, 2, function(x) ifelse(x == "2,2", "6", x)) 

#cnv_sim_luad <- cnv_sim_luad %>% as.tibble() %>% mutate_if(is.character, as.numeric)
#hist(rowMeans(cnv_sim_luad))

### Frequency histogram ###
#lihc_cnv <- apply(lihc_cnv, 2, function(x) ifelse(x > 10, 10, x)) 
#hist(rowMeans(lihc_cnv),
     #main = "LIHC", 
     #xlab = "CN state",
     #ylab = "Proportion",
     #col = "#E1DEFC",
     #prob = TRUE,
     #breaks = 6)


### CN saturation curve ###
#cn <- seq(from = 0, to = 30, 1)
#cn <- cn/2
#sigma <- (2*exp(cn)) / (1+exp(cn))
#plot(x = log(cn), y = log(sigma), col = 4, type = "b", main = "CN saturation curve")

#cn <- seq(from = 0, to= 30, 1)
#rna <- seq(from = 20, to = 600, 50)
#rna_cn <- rna * sigma
#plot(x = log(cn), y = log(rna_cn), col = 4, type = "b", main = "RNA vs CN")
#cnv <- as.data.frame(cnv)
#rna_cnv <- rna_counts * cnv

