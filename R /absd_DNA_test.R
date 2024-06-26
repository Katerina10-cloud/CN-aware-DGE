### CN-aware ABCD-DNA test ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "gridExtra", "ggpubr", "ggrepel", "ggvenn", "ggpointdensity",
          "DESeq2", "edgeR")
sapply(pkgs, require, character.only = TRUE)


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

# Fit the NB GLMs #
design <- model.matrix(~1+condition, data=metadata)
edger.obj <- edgeR::DGEList(rna_cnv)
edger.obj <- edgeR::calcNormFactors(edger.obj, method="TMM")
#edger.obj <- edgeR::estimateDisp(edger.obj, design)
fit <- edgeR::glmFit(edger.obj, design, dispersion=0.05)
lrt <- edgeR::glmLRT(fit, coef=2)
res <- topTags(lrt, n=1000)$table
#glmz <- -sign(lrt$table$logFC)*abs(qnorm(lrt$table$PValue/2))

#res <- lrt[["table"]]
#saveRDS(res, file = "CN-aware-DGE/results/rna_edge.RDS")

# Fit the CN-aware NB GLMs #
data_obj <- edgeR::DGEList(rna_cnv, group = metadata$condition)
data_obj <- edgeR::calcNormFactors(data_obj)
offset <- outer( rep(1,nrow(data_obj)), getOffset(data_obj)) + log(cnv)
#data_obj <- edgeR::estimateDisp(data_obj, design)
fit_adj <- edgeR::glmFit(y=rna_cnv, design=design, offset=offset,dispersion = 0.05)
#test <- edgeR::glmTreat(fit, coef=2)
lrt_adj <- edgeR::glmLRT(fit_adj, coef=2)
res_adj <- topTags(lrt_adj, n=1000)$table
#res_adj <- lrt_adj[["table"]]

rownames_idx <- match(rownames(res), rownames(res_adj))
res_adj <- res_adj[rownames_idx,]

#saveRDS(res_adj, file = "CN-aware-DGE/results/rna_edge_adj.RDS")



# Plot results #

plotMD(lrt, main="no CN normalization")
plotMD(lrt_adj, main = "wit CN normalization")
abline(h=c(-0.5, 0.5), col="blue")

# Trended biases of CN on DGE #
cnv_tumor <- cnv[,11:20]
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
  geom_smooth(method = "lm", color = "blue", se = T)+
  theme(legend.position = 'bottom') +
  labs(title = "no CN normalization", x ="log2(CN ratio)", y="log2FC") +
  theme_classic()
scatter1

scatter2 = ggplot(plot_data_adj, aes(x = log2(cnv_mean), y = logFC)) +
  geom_point(size=0.8) +
  geom_smooth(method = "lm", color = "blue", se = T)+
  theme(legend.position = 'bottom') +
  labs(title = "with CN normalization", x ="log2(CN ratio)", y="log2FC") +
  theme_classic()
scatter1+scatter2


# Log2FC comparison #
res_lfc <- res %>% 
  dplyr::select(logFC) %>% 
  dplyr::rename(logFC=logFC)
res_lfc_adj <- res_adj %>% 
  dplyr::select(logFC) %>% 
  dplyr::rename(logFC_adj = logFC)

plot_data <- merge(res_lfc, res_lfc_adj, by = "row.names")

comparison <- ggplot(plot_data, aes(x=logFC_adj, y=logFC)) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  labs(title = "Log2FC comparison")+
  xlab("log2FC with CN normalization") +
  ylab ("log2FC no CN normalization") +
  theme_classic()+
  theme(legend.position="none")
comparison

