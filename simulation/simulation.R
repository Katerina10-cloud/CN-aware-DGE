###------------------------------------------------------------###
### Simulate RNAseq counts ###
###------------------------------------------------------------###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA")

library(tidyverse)
library(compcodeR)

# Using compcodeR simulator
rna_counts_sim <- compcodeR::generateSyntheticData(dataset = "rna_counts_sim", n.vars = 15000, 
                                                   samples.per.cond = 36, n.diffexp = 0, 
                                                   repl.id = 1, seqdepth = 1e7, 
                                                   fraction.upregulated = 0.0, 
                                                   between.group.diffdisp = TRUE, 
                                                   filter.threshold.total = 1, 
                                                   filter.threshold.mediancpm = 0, 
                                                   fraction.non.overdispersed = 0, 
                                                   relmeans = "auto",
                                                   dispersions = "auto",
                                                   random.outlier.high.prob = 20,
                                                   random.outlier.low.prob = 20,
                                                   output.file = "rna_counts_sim.rds")

rna_counts_sim <- rna_counts_sim@count.matrix
rna_normal <- luad_rna_sd %>% as.data.frame() %>% select(1:45)
rna_tumor <- luad_rna_sd %>% as.data.frame() %>% select(46:90)
rna_counts <- cbind(brca_rna_normal, brca_rna_tum)

# Using OmicsSIMLA
rna_counts_sim1 <- read.table("rna_luad1.tsv", header = TRUE, sep = "")
rna_counts_sim2 <- read.table("rna_luad2.tsv", header = TRUE, sep = "")

rna_counts_sim2 <- rna_counts_sim2[,c(1:5000)]
rna_counts_sim2 <- rna_counts_sim2 %>% remove_rownames %>% column_to_rownames(var="FAM")
colnames(rna_counts_sim) <- paste0("S", 1:(ncol(rna_counts_sim)))
rownames(rna_counts_sim) <- paste("G", 1:(nrow(rna_counts_sim)))

# Generate metadata
metadata <- data.frame(patID = colnames(rna_counts_sim),
                       condition = rep(c("A", "B"), each = 36))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 


###-----------------------------------------------------------###
### Simulate CN data ###
###-----------------------------------------------------------###

### Generate heterogeneous CN data ### 
cnv_0 <- sapply(1:36, function(x) sample(x=c(0.5,1,2,3), size = 500, replace=TRUE, prob = c(.50, .30, .10, .10)))
cnv_1 <- sapply(1:36, function(x) sample(x=c(0.5,1,2,3), size = 500, replace=TRUE, prob = c(.20, .60, .10, .10)))
#cnv_3 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3,4), size = 700, replace=TRUE, prob = c(.05, .05, .15, .70, .05)))
#cnv_4 <- sapply(1:50, function(x) sample(x=c(0.5,1,2,3,4,5), size = 700, replace=TRUE, prob = c(.05, .05, 0.10, 0.10, .60, .10)))
#cnv_5 <- sapply(1:50, function(x) sample(x=c(1,2,3,4,5), size = 400, replace=TRUE, prob = c(.05, .05, 0.10, .10, .70)))
cnv_2 <- sapply(1:36, function(x) sample(x=c(1,2,3,4), size = 11000, replace=TRUE, prob = c(.05, .70, .10, .10)))

rna_normal <- rna_counts_sim[1:36]
rna_tumor <- rna_counts_sim[37:72]

cnv_luad <- cnv_luad[1:3000,]
cnv_tumor <- rbind(cnv_luad, cnv_0, cnv_1, cnv_2) %>% as.matrix()
cnv_norm <- matrix(2, nrow(rna_normal), 36)

colnames(cnv_norm) <- colnames(rna_normal)
rownames(cnv_norm) <- rownames(rna_normal)

#cnv <- list(cnv_norm, cnv_tumor) 
cnv <- cbind(cnv_tumor, cnv_norm)


## Prepare CN data ##
cnv <- apply(cnv, 1, function(x) x/2)
cnv <- apply(cnv, 1, function(x) x+10e-9)

rna_cnv <- rna_nocnv * cnv


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
lihc_cnv <- apply(lihc_cnv, 2, function(x) ifelse(x > 10, 10, x)) 
hist(rowMeans(lihc_cnv),
     main = "LIHC", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 6)

write.csv(metadata, file = "metadata.csv")
write.csv(rna_cnv, file = "rna_cnv.csv")

save(cnv, file = "sim_OmicsSIMLA/cnv_sim.Rdata")


### CN saturation curve ###
cn <- seq(from = 0, to = 30, 1)
cn <- cn/2
sigma <- (2*exp(cn)) / (1+exp(cn))
plot(x = log(cn), y = log(sigma), col = 4, type = "b", main = "CN saturation curve")

cn <- seq(from = 0, to= 30, 1)
rna <- seq(from = 20, to = 600, 50)
rna_cn <- rna * sigma
plot(x = log(cn), y = log(rna_cn), col = 4, type = "b", main = "RNA vs CN")


cnv <- as.data.frame(cnv)
rna_cnv <- rna_counts * cnv


