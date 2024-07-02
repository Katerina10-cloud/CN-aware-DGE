### CN-aware ABCD-DNA test ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "gridExtra", "ggpubr", "ggrepel", "ggvenn", "ggpointdensity",
          "DESeq2", "edgeR")
sapply(pkgs, require, character.only = TRUE)

### Test on simulated data ###

# RNA counts simulation #
m = 20
sizeFactors = rep(1, m)
deseq_sim <- DESeq2::makeExampleDESeqDataSet(
  n = 1000,
  m = 20,
  betaSD = 0,
  interceptMean = 6,
  interceptSD = 2,
  dispMeanRel = function(x) 6/x + 0.6,
  sizeFactors = sizeFactors
)
rna_counts <- data.frame(deseq_sim@assays@data@listData[["counts"]])

# Generate metadata #
metadata <- data.frame(patID = colnames(rna_counts),
                       condition = rep(c("0", "1"), each = 10))
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

#lib.size <- colSums(rna_cnv, na.rm=T)

# Fit the NB GLMs #
design <- model.matrix(~1+condition, data=metadata)
edger.obj <- edgeR::DGEList(rna_cnv)
edger.obj <- edgeR::calcNormFactors(edger.obj, method="TMM")
edger.obj <- edgeR::estimateDisp(edger.obj, design)
#fit <- edgeR::glmFit(edger.obj, design, dispersion=0.05)
fit <- edgeR::glmFit(edger.obj, design)
lrt <- edgeR::glmLRT(fit, coef=2)
res <- topTags(lrt, n=1000)$table
#glmz <- -sign(lrt$table$logFC)*abs(qnorm(lrt$table$PValue/2))

#res <- lrt[["table"]]
#saveRDS(res, file = "CN-aware-DGE/results/rna_edge.RDS")

#lib.size <- edger.obj[["samples"]][["lib.size"]]
#norm.factors <- edger.obj[["samples"]][["norm.factors"]]
#lib.size <- lib.size*norm.factors
#lib.size <- log(lib.size)


# Fit the CN-aware NB GLMs #
data_obj <- edgeR::DGEList(rna_cnv, group = metadata$condition)
data_obj <- edgeR::calcNormFactors(data_obj)
offset <- outer( rep(1,nrow(data_obj)), getOffset(data_obj)) + log(cnv)
#data_obj <- edgeR::estimateDisp(data_obj, design)
#fit_adj <- edgeR::glmFit(y=rna_cnv, design=design, offset=offset,dispersion = 0.05)
fit_adj <- edgeR::glmFit(y=rna_cnv, design=design, offset=offset, dispersion = edger.obj[["tagwise.dispersion"]])
#test <- edgeR::glmTreat(fit, coef=2)
lrt_adj <- edgeR::glmLRT(fit_adj, coef=2)
res_adj <- topTags(lrt_adj, n=1000)$table
#res_adj <- lrt_adj[["table"]]

rownames_idx <- match(rownames(res), rownames(res_adj))
res_adj <- res_adj[rownames_idx,]

#saveRDS(res_adj, file = "CN-aware-DGE/results/rna_edge_adj.RDS")

### Test real data ###
load("TCGA/lung_cancer/LUAD/data/clinical_luad.Rdata")
load("TCGA/lung_cancer/LUAD/data/luad_cnv_tumor.Rdata")
load("TCGA/lung_cancer/LUAD/data/luad_rna_tumor.Rdata")
load("TCGA/lung_cancer/LUAD/data/luad_rna_normal.Rdata")

clinical_luad <- clinical_luad %>% 
  dplyr::rename(patID=bcr_patient_barcode, stage=stage_event_pathologic_stage) %>% 
  dplyr::select(patID,stage)

clinical_luad <- clinical_luad[ clinical_luad$stage %in% c("Stage IIIA"),]
luad_cnv_tumor <- luad_cnv_tumor[ ,colnames(luad_cnv_tumor) %in% clinical_luad$patID]

colnames(luad_rna_norm) <- substr(colnames(luad_rna_tum), 1, 12)
colnames(luad_rna_tum) <- substr(colnames(luad_rna_tum), 1, 12)

luad_rna_norm <- luad_rna_norm[ ,colnames(luad_rna_norm) %in% colnames(luad_cnv_tumor)]
luad_rna_tum <- luad_rna_tum[ ,colnames(luad_rna_tum) %in% colnames(luad_cnv_tumor)]

x <- colnames(luad_rna_norm)
names(luad_rna_norm) <- paste(x,"-11A")

x <- colnames(luad_rna_tum)
names(luad_rna_tum) <- paste(x,"-01A")

luad_rna <- cbind(luad_rna_norm, luad_rna_tum)

luad_rna <- luad_rna[which(rowSums(luad_rna)>10000),]

luad_cnv_tumor <- luad_cnv_tumor[ rownames(luad_cnv_tumor) %in% rownames(luad_rna),]

rownames_idx <- match(rownames(luad_rna), rownames(luad_cnv_tumor))
luad_cnv_tumor <- luad_cnv_tumor[rownames_idx,]

luad_cnv_normal <- matrix(2, nrow(luad_cnv_tumor), 10) %>% as.data.frame()
colnames(luad_cnv_normal) <- colnames(luad_cnv_tumor)
rownames(luad_cnv_normal) <- rownames(luad_cnv_tumor)

x <- colnames(luad_cnv_normal)
names(luad_cnv_normal) <- paste(x,"-11A")

x <- colnames(luad_cnv_tumor)
names(luad_cnv_tumor) <- paste(x,"-01A")

luad_cnv <- cbind(luad_cnv_normal, luad_cnv_tumor)

luad_cnv <- apply(luad_cnv, 2, function(x) ifelse(x > 10, 10, x)) 
luad_cnv <- luad_cnv/2

#luad_rna <- apply(luad_rna, 2, function(x) ifelse(x == 0, 1, x)) 

hist(rowMeans(luad_cnv_tumor),
     main = "LUAD", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 10)

colnames(luad_rna) <- paste0("sample", 1:(ncol(luad_rna)))
colnames(luad_cnv) <- colnames(luad_rna)

metadata <- data.frame(patID = colnames(luad_rna),
                       condition = rep(c("0", "1"), each = 10))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 

# DE test #

design <- model.matrix(~1+condition, data=metadata)
edger.obj <- edgeR::DGEList(luad_rna)
edger.obj <- edgeR::calcNormFactors(edger.obj, method="TMM")
edger.obj <- edgeR::estimateDisp(edger.obj, design)
fit <- edgeR::glmFit(edger.obj, design)
lrt <- edgeR::glmLRT(fit, coef=2)
res <- topTags(lrt, n=10057)$table

# Fit the CN-aware NB GLMs #
design <- model.matrix(~1+condition, data=metadata)
data_obj <- edgeR::DGEList(luad_rna)
data_obj <- edgeR::calcNormFactors(data_obj, method="TMM")
offset <- outer( rep(1,nrow(luad_rna)), getOffset(data_obj)) + log(luad_cnv) %>% as.data.frame()

offset <- offset %>% filter_all(all_vars(!is.infinite(.))) %>% as.matrix()
luad_rna <- luad_rna[ rownames(luad_rna) %in% rownames(offset),]
luad_cnv <- luad_cnv[ rownames(luad_cnv) %in% rownames(luad_rna),]

#Perform some checks #
#any(sapply(offset, is.infinite))
#sum(offset < 0)
#sf <- data_obj[["samples"]][["norm.factors"]]
#any(sapply(sf, is.infinite))

fit_adj <- edgeR::glmFit(luad_rna, design, offset, dispersion = edger.obj[["tagwise.dispersion"]])
lrt_adj <- edgeR::glmLRT(fit_adj, coef=2)
res_adj <- edgeR::topTags(lrt_adj, n=10057)$table

rownames_idx <- match(rownames(res), rownames(res_adj))
res_adj <- res_adj[rownames_idx,]


# Plot results #

plotMD(lrt, main="no CN normalization")
plotMD(lrt_adj, main = "wit CN normalization")
abline(h=c(-0.5, 0.5), col="blue")

# Trended biases of CN on DGE #
cnv_tumor <- luad_cnv[,11:20]
cnv_tumor <- cnv_tumor %>% 
  as.data.frame() %>% 
  mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  select(cnv_mean)

res_lfc <- res %>% select(logFC)
res_lfc_adj <- res_adj %>% select(logFC)
plot_data <- merge(cnv_tumor, res_lfc, by = "row.names")
plot_data_adj <- merge(cnv_tumor, res_lfc_adj, by = "row.names")

scatter1 = ggplot(plot_data, aes(x = log2(cnv_mean), y = logFC)) +
  geom_point(size=0.8) +
  #geom_smooth(method = "lm", color = "blue", se = T)+
  geom_smooth()+
  theme(legend.position = 'bottom') +
  labs(title = "no CN normalization", x ="log2(CN ratio)", y="log2FC") +
  theme_classic()
scatter1

scatter2 = ggplot(plot_data_adj, aes(x = log2(cnv_mean), y = logFC)) +
  geom_point(size=0.8) +
  #geom_smooth(method = "lm", color = "blue", se = T)+
  geom_smooth()+
  theme(legend.position = 'bottom') +
  labs(title = "with CN normalization", x ="log2(CN ratio)", y="log2FC") +
  theme_classic()

gridExtra::grid.arrange(scatter1, scatter2, nrow = 1)

# Log2FC comparison #
res_lfc <- res %>% 
  dplyr::select(logFC) %>% 
  dplyr::rename(logFC=logFC)
res_lfc_adj <- res_adj %>% 
  dplyr::select(logFC) %>% 
  dplyr::rename(logFC_adj = logFC)

plot_lfc <- merge(res_lfc, res_lfc_adj, by = "row.names")

comparison_lfc <- ggplot(plot_lfc, aes(x=logFC_adj, y=logFC)) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "Log2FC comparison")+
  xlab("log2FC with CN normalization") +
  ylab ("log2FC no CN normalization") +
  theme_classic()+
  theme(legend.position="none")
comparison_lfc

# p-value comparison #
res_pval <- res %>% 
  dplyr::select(FDR) %>% 
  dplyr::rename(FDR=FDR)
res_pval_adj <- res_adj %>% 
  dplyr::select(FDR) %>% 
  dplyr::rename(FDR_adj = FDR)

plot_pval <- merge(res_pval, res_pval_adj, by = "row.names")

comparison_pval <- ggplot(plot_pval, aes(x=-log10(FDR_adj), y=-log10(FDR))) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "p-value comparison")+
  xlab("FDR with CN normalization") +
  ylab ("FDR no CN normalization") +
  theme_classic()+
  theme(legend.position="none")
comparison_pval

gridExtra::grid.arrange(comparison_lfc, comparison_pval, nrow = 1)


