###------------------------------------------------------------###
### Simulate RNAseq counts ###
###------------------------------------------------------------###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/de_fit_Python")

library(tidyverse)
library(compcodeR)

rna_counts_sim <- compcodeR::generateSyntheticData(dataset = "rna_counts_sim", n.vars = 15000, 
                                                   samples.per.cond = 180, n.diffexp = 0, 
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
rna_normal <- rna_counts_sim %>% as.data.frame() %>% select(1:90)
rna_tumor <- rna_counts_sim %>% as.data.frame() %>% select(91:180)
metadata <- data.frame(patID = colnames(rna_nocnv),
                       condition = rep(c("A", "B"), each = 90))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 


###-----------------------------------------------------------###
### Simulate CN data ###
###-----------------------------------------------------------###

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


### Generate heterogeneous CN data ### 
cnv_0 <- sapply(1:36, function(x) sample(x=c(0.5,1,2,3), size = 500, replace=TRUE, prob = c(.50, .30, .10, .10)))
cnv_1 <- sapply(1:36, function(x) sample(x=c(0.5,1,2,3), size = 500, replace=TRUE, prob = c(.20, .60, .10, .10)))
cnv_3 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3,4), size = 700, replace=TRUE, prob = c(.05, .05, .15, .70, .05)))
cnv_4 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3,4,5), size = 700, replace=TRUE, prob = c(.05, .05, 0.10, 0.10, .60, .10)))
cnv_5 <- sapply(1:50, function(x) sample(x=c(1,2,3,4,5), size = 400, replace=TRUE, prob = c(.05, .05, 0.10, .10, .70)))
cnv_2 <- sapply(1:36, function(x) sample(x=c(1,2,3,4,5), size = 11000, replace=TRUE, prob = c(.05, .80, .05, .05, .05)))

cnv_tumor <- rbind(laml_cnv_sd1, cnv_heterog_2)
cnv_norm <- matrix(1, nrow(rna_normal), 90)

colnames(cnv_heterog_2) <- colnames(rna_tumor)
rownames(cnv_heterog_2) <- rownames(rna_tumor)
 
cnv_tumor <- cnv_tumor/2   
cnv <- cbind(cnv_norm, cnv_tumor)
cnv <- cnv + 10e-9
rna_nocnv <- rna_counts_sim * cnv


## Reverse dosage simulation ##
rna_tum <- rna_nocnv[,37:72]

gr0_rev <- rna_tum[3001:3050,]*4
gr1_rev <- rna_tum[3501:3550,]*4
gr3_rev <- rna_tum[1:350,]/4
#gr4_rev <- rna_tum[1901:2000,]/4
#gr5_rev <- rna_tum[2601:2650,]/4

## Other signal simulation ##
gr0_oth <- rna_tum[3051:3100,]*0.5
gr1_oth <- rna_tum[3551:3600,]*0.5
gr3_oth <- rna_tum[351:450,]*1.5
gr4_oth <- rna_tum[451:550,]*2
#gr5_oth <- rna_tum[2651:2700,]*2

## Genes with CN signal ##
gr_cn <- rna_tum[551:3000,]
gr0_cn <- rna_tum[3101:3500,]
gr1_cn <- rna_tum[3601:4000,]
#gr3_cn <- rna_tum[1401:1900,]
#gr4_cn <- rna_tum[2101:2600,]
#gr5_cn <- rna_tum[2701:3000,]
gr2_cn <- rna_tum[4001:15000,]


#Genes with mixed signal (CN + other signal type)
rna_mixed <- rbind(gr3_rev, gr3_oth, gr4_oth, gr_cn, gr0_rev, gr0_oth, gr0_cn, gr1_rev, gr1_oth, gr1_cn, gr2_cn)
genes <- rownames(rna_mixed)
rna_normal <- rna_normal[genes,]
cnv_tumor <- cnv_tumor[genes,]
rna_mix <- cbind(rna_normal, rna_mixed)


## Simulated CN data preprocessing OMICSSimLA ##

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/CNV-informed-DGE-modelling")

cnv_sim_luad <- read.table(file = "data/GbmSim1.tsv", header = FALSE, sep = "")
cnv_sim_luad <- cnv_sim_luad[-c(2:3)]
cnv_sim_luad <- cnv_sim_luad[-which(duplicated(cnv_sim_luad$V1)), ]
cnv_sim_luad <- cnv_sim_luad %>% remove_rownames %>% column_to_rownames(var="V1")
colnames(cnv_sim_luad) <- paste0("G", 1:(ncol(cnv_sim_luad)-1))
cnv_sim_luad <- t(cnv_sim_luad)

# Assigning CN states #
cnv_sim_luad <- apply(cnv_sim_luad, 2, function(x) ifelse(x == "2,2", "6", x)) 

cnv_sim_luad <- cnv_sim_luad %>% as.tibble() %>% mutate_if(is.character, as.numeric)
hist(rowMeans(cnv_sim_luad))

### Frequency histogram ###
hist(rowMeans(cnv_tumor),
     main = "CNV frequency TCGA-LAML (3500 genes)", 
     xlab = "CN state",
     breaks = 20)


write.csv(metadata, file = "de_fit_Python/data_simulation/sim2_real_cnv/aml/metadata.csv")
write.csv(rna_nocnv, file = "de_fit_Python/data_simulation/sim2_real_cnv/aml/rna_nocnv.csv")

#laml_cnv_sd1 <- laml_cnv_sd1[1:3000,]
#colnames(cnv_heterog_2) <- colnames(rna_tumor)
#cnv_heterog_2 <- as.data.frame(cnv_heterog_2)
#rownames(cnv_heterog_2) <- paste0(1:(nrow(cnv_heterog_2)))


### CN saturation curve ###
cn <- seq(from = 0, to = 30, 1)
cn <- cn/2
sigma <- (2*exp(cn)) / (1+exp(cn))
plot(x = log(cn), y = log(sigma), col = 4, type = "b", main = "CN saturation curve")

cn <- seq(from = 0, to= 30, 1)
rna <- seq(from = 0, to = 600, 20)
rna_cn <- rna * cn
plot(x = log(cn), y = log(rna), col = 4, type = "b", main = "RNA vs CN")






