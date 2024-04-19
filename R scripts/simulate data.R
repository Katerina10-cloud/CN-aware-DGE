#Make a simulated DESeqDataSet 
#BiocManager::install("DESeq2")

setwd("C:/Users/rifug/Documents/CNV-informed-DGE-modelling")

library(DESeq2)
library(tidyverse)

### Generate RNA counts data ###
dds <- DESeq2::makeExampleDESeqDataSet(
  n = 20000,
  m = 100,
  betaSD = 0,
  interceptMean = 8,
  interceptSD = 2,
  dispMeanRel = function(x) 2/x + 1.6,
  #sizeFactors = rep(1, m)
  )

rna_counts <- dds@assays@data@listData[["counts"]]
rna_normal <- rna_counts %>% as.data.frame() %>% select(1:50)
rna_tumor <- rna_counts %>% as.data.frame() %>% select(51:100)
metadata <- data.frame(patID = colnames(rna_counts), 
                       condition = rep(c("A", "B"), each = 50)) 
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="patID") 

### Generate Copy Number homogeneous data ###
group0 <- rna_normal[1:1000,]
group1 <- rna_normal[1001:2500,]
group2 <- rna_normal[2501:12500,]
group3 <- rna_normal[12501:16500,]
group4 <- rna_normal[16501:18500,]
group5 <- rna_normal[18501:20000,]

cnv0 <- matrix(0.5, nrow(group0), 50)
cnv1 <- matrix(1, nrow(group1), 50)
cnv2 <- matrix(2, nrow(group2), 50)
cnv3 <- matrix(3, nrow(group3), 50)
cnv4 <- matrix(4, nrow(group4), 50)
cnv5 <- matrix(5, nrow(group5), 50)

cnv_tumor <- rbind(cnv0, cnv1, cnv2, cnv3, cnv4, cnv5)
rna_normal <- rbind(group0, group1, group2, group3, group4, group5)


### Heterogeneous CN data ### 
cnv_heterog_0 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3), size=1000, replace=TRUE, prob=c(.50,.30,.10,.10)))
cnv_heterog_1 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3), size=2000, replace=TRUE, prob=c(.20,.60,.10,.10)))
cnv_heterog_2 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3,4), size=10000, replace=TRUE, prob=c(.05,.05,.70,.15,.05)))
cnv_heterog_3 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3,4,5), size=3000, replace=TRUE, prob=c(.05,.05,.10,.60,.10,.10)))
cnv_heterog_4 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3,4,5), size=2500, replace=TRUE, prob=c(.05,.05,.10,.10,.60,.10)))
cnv_heterog_5 <- sapply(1:50, function(x) sample(x=c(1,2,3,4,5), size=1500, replace=TRUE, prob=c(.05,.05,.10,.10,.70)))

cnv_tumor <- rbind(cnv_heterog_0, cnv_heterog_1, cnv_heterog_2, cnv_heterog_3, cnv_heterog_4, cnv_heterog_5)

cnv_norm <- matrix(1, nrow(rna_normal), 50)

colnames(cnv_norm) <- colnames(rna_normal)
rownames(cnv_norm) <- rownames(rna_normal)

cnv_tumor <- cnv_tumor/2   
cnv <- cbind(cnv_tumor, cnv_norm)
cnv <- cnv + 10e-9
rna_mixed_cnv <- rna_mixed * cnv


## Reverse dosage simulation ##
rna_tum <- rna[,51:100]

gr0_rev <- rna_tum[1:50,]*4
gr1_rev <- rna_tum[501:550,]*4
gr3_rev <- rna_tum[1201:1300,]/4
gr4_rev <- rna_tum[1901:2000,]/4
gr5_rev <- rna_tum[2601:2650,]/4

## Other signal simulation ##
gr0_oth <- rna_tum[51:100,]*0.5
gr1_oth <- rna_tum[551:600,]*0.5
gr3_oth <- rna_tum[1301:1400,]*1.5
gr4_oth <- rna_tum[2001:2100,]*2
gr5_oth <- rna_tum[2651:2700,]*2

## Genes with CN signal ##
gr0_cn <- rna_tum[101:500,]
gr1_cn <- rna_tum[601:1200,]
gr3_cn <- rna_tum[1401:1900,]
gr4_cn <- rna_tum[2101:2600,]
gr5_cn <- rna_tum[2701:3000,]
gr2_cn <- rna_tum[3001:15000,]

deg_up <- subset(res1_nocnv, res1_nocnv$padj < 0.05 & res1_nocnv$log2FoldChange > 0.5)
deg_down <- subset(res1_nocnv, res1_nocnv$padj < 0.05 & res1_nocnv$log2FoldChange < -0.5)
deg <- rbind(deg_up, deg_down)
deg_cnv <- rownames(deg)

rna_tum <- rna %>% as.data.frame() %>% select(51:100)
genes_nocnv <- rna_tum[!rownames(rna_tum) %in% rownames(deg),]
genes_cnv <- rna_tum[rownames(rna_tum) %in% rownames(deg),]

## Diploid genes ##
genes_diploid <- genes_nocnv[1870:11747,]

#Genes with mixed signal (CN + other signal type)
rna_mixed <- rbind(gr0_rev, gr0_oth, gr0_cn,
                   gr1_rev, gr1_oth, gr1_cn,
                   gr3_rev, gr3_oth, gr3_cn,
                   gr4_rev, gr4_oth, gr4_cn,
                   gr5_rev, gr5_oth, gr5_cn, 
                   gr2_cn)

genes <- rownames(rna_mixed)
rna_normal <- rna_normal[genes,]
cnv_tumor <- cnv_tumor[genes,]
rna_mixed <- cbind(rna_normal, rna_mixed)


## Simulated data preprocessing OMICSSimLA ##

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

cnv_sim_lusc <- read.table(file = "LuscSim_cnv.txt", header = FALSE, sep = "")
cnv_sim_lusc <- cnv_sim_lusc[-c(2:3)]
cnv_sim_lusc <- cnv_sim_lusc %>% remove_rownames %>% column_to_rownames(var="V1")
colnames(cnv_sim_lusc) <- paste0("G", 1:(ncol(cnv_sim_lusc)-1))
cnv_sim_lusc <- t(cnv_sim_lusc)

# Assigning CN states #
cnv_sim_lusc <- apply(cnv_sim_lusc, 2, function(x) ifelse(x == "2,2", "6", x)) 

cnv_sim_lusc <- cnv_sim_lusc %>% as.tibble() %>% mutate_if(is.character, as.numeric)
hist(rowMeans(cnv_sim_lusc))

# Sampling data generation #

cnv_0 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3), size = 500, replace=TRUE, prob = c(.50, .30, .10, .10)))
cnv_1 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3), size = 700, replace=TRUE, prob = c(.20, .60, .10, .10)))
cnv_3 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3,4), size = 700, replace=TRUE, prob = c(.05, .05, .15, .70, .05)))
cnv_4 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3,4,5), size = 700, replace=TRUE, prob = c(.05, .05, 0.10, 0.10, .60, .10)))
cnv_5 <- sapply(1:50, function(x) sample(x=c(1,2,3,4,5), size = 400, replace=TRUE, prob = c(.05, .05, 0.10, .10, .70)))
cnv_2 <- sapply(1:50, function(x) sample(x=c(1,2,3,4,5), size = 12000, replace=TRUE, prob = c(.05, .80, .05, .05, .05)))

cnv_tumor <- rbind(cnv_0, cnv_1, cnv_3, cnv_4, cnv_5, cnv_2)


hist(rowMeans(cnv_2),
     main = "CNV simulation frequency (2884 genes)", 
     xlab = "CN state",
     breaks = 6)

# Generate RNAseq counts data #
library(compcodeR)

rna_counts_sim <- compcodeR::generateSyntheticData(dataset = "rna_counts_sim", n.vars = 15000, 
                                                   samples.per.cond = 50, n.diffexp = 0, 
                                                   repl.id = 1, seqdepth = 1e7, 
                                                   fraction.upregulated = 0.0, 
                                                   between.group.diffdisp = FALSE, 
                                                   filter.threshold.total = 1, 
                                                   filter.threshold.mediancpm = 0, 
                                                   fraction.non.overdispersed = 0, 
                                                   relmeans = "auto",
                                                   dispersions = "auto",
                                                   output.file = "rna_counts_sim.rds")

rna_counts_sim <- rna_counts_sim@count.matrix
rna_normal <- rna_counts_sim %>% as.data.frame() %>% select(1:50)
rna_tumor <- rna_counts_sim %>% as.data.frame() %>% select(51:100)
metadata <- data.frame(patID = colnames(rna_counts_sim),
                       condition = rep(c("A", "B"), each = 50))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID")    

write.csv(metadata, file = "/Users/katsiarynadavydzenka/Documents/PhD_AI/model_fit_Python/data_simulation/metadata.csv")
write.csv(rna_mixed_cnv, file = "/Users/katsiarynadavydzenka/Documents/PhD_AI/model_fit_Python/data_simulation/sim_1/rna_mixed_cnv.csv")
save(rna_mixed, file = "rna_mixed_sim.Rdata")



