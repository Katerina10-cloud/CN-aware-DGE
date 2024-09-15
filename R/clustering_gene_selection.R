### Identify homogeneous groups of patients based on their CNV profiles ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "DESeq2", "colorspace", 
          "ggpubr", "ggpointdensity", "ggeasy")
sapply(pkgs, require, character.only = TRUE)

source("CN-aware-DGE/R/utils.R")
#set.seed(12345)

# Input data
data_path <- "TCGA/lung/LUAD/cnv_tumor.RDS"
dataset_name <- "LUAD_cnv"

# Clustering patients 
clustering_res <- clustering_patients_cnv(dataset_name, data_path)
clustering_res[[2]]
cnv_tumor <- clustering_res[[1]]

#Select patients for boxplot
#cnv_tumor <- as.matrix(t(cnv_tumor))
#cnv_tumor <- cnv_tumor %>% as.data.frame() %>% dplyr::select(7,5,13,14,15,18,22,24,25,28)
#cnv_filt <- apply(cnv_tumor, 2, function(x) ifelse(x > 15, 15, x)) 

# Cluster selection and analysis
cnv_filt <- cnv_tumor[cnv_tumor$Cluster %in% c(3),]
cnv_filt <- subset(cnv_filt, select=-c(Cluster))
cnv_filt <- as.matrix(t(cnv_filt))
cnv_filt <- apply(cnv_filt, 2, function(x) ifelse(x > 10, 10, x))

saveRDS(cnv_filt, file = "TCGA/lung/cnv_filt.RDS")

hist(rowMeans(cnv_filt),
     main = "LUAD", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 15)

cnv_mean <- cnv_filt %>% 
  as.data.frame() %>% 
  mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  select(cnv_mean)

cnv_mean$geneID <- rownames(cnv_mean)


## RNA data processing ##

data_path <- "TCGA/lung/LUAD/rna_counts.RDS"
dataset_name <- "LUAD_rna"

rna <- rna_processing(dataset_name, data_path, cnv_filt)
rna_norm <- rna[[1]]
rna_tum <- rna[[2]]

# Gene selection #

# Exclude genes with low expression in normal tissue #
low_expression_threshold <- 15
expression_summary <- data.frame(
  Gene = rownames(rna_norm),
  MeanExpression = rowMeans(rna_norm)
)

summary(expression_summary$MeanExpression)

filtered_genes <- expression_summary %>%
  filter(MeanExpression > low_expression_threshold)

rna_norm <- rna_norm[filtered_genes$Gene, ]
rna_tum <- rna_tum[filtered_genes$Gene, ]
rna <- cbind(rna_norm, rna_tum)

# Expression variability check #
gene_iqr <- apply(rna, 1, IQR)
gene_sd <- apply(rna, 1, sd)

variability_summary <- data.frame(
  Gene = rownames(rna),
  IQR = gene_iqr,
  SD = gene_sd
)

iqr_threshold <- quantile(variability_summary$IQR, 0.75)
#sd_threshold <- quantile(variability_summary$SD, 0.75)

filtered_genes_iqr <- variability_summary %>%
  filter(IQR > iqr_threshold)

#filtered_genes_sd <- variability_summary %>%
#filter(SD > sd_threshold)

rna_filt <- rna[filtered_genes_iqr$Gene, ]


## Correlation analysis ##

# Filter DEGs
pval_cut <- 0.05
lfc_cut <- 1.0

res_deg <- res_naive_edge %>% 
  dplyr::filter(FDR < pval_cut, abs(logFC) > lfc_cut) %>% 
  dplyr::mutate(isDE = (abs(logFC) >= lfc_cut) & (FDR <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(logFC > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::select(logFC, FDR, isDE, DEtype) 

volcano <-  res_deg %>% 
  ggplot(mapping = aes(x=logFC, y=-log10(FDR), col=DEtype)) +
  geom_point(size=.8) +
  theme_bw() +
  scale_color_manual(values = de_gene_colors) +
  #ggh4x::facet_nested(tool~method, scale="free") +
  #facet_wrap(~method, nrow=1, scale = "free")+
  ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col="") +
  ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  ggplot2::geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'bottom')
volcano


# Calculate Pearson correlation between CNV and expression for a gene
rna_filt <- rna[ rownames(rna) %in% rownames(res_deg) , ]
rna_log_normalized <- rna_filt %>%  as.matrix() %>% DESeq2::vst()
rna_tumor <- rna_log_normalized[,21:40]

#common_genes <- intersect(rownames(cnv_filt), rownames(rna))
#cnv_filt <- cnv_filt[common_genes, ]
#rna <- rna[common_genes, ]

gene_correlations <- sapply(rownames(cnv_filt), function(gene) {
  calculate_correlation(cnv_filt[gene, ], rna_tumor[gene, ])
})

cor_results <- data.frame(
  Gene = rownames(cnv_mean),
  Correlation = gene_correlations
)

rownames(cor_results) <- cor_results$Gene

ggplot(data = cor_results, aes(x = Correlation)) +
  geom_histogram(color = "white", fill = "lightblue", bins = 30)+
  ggplot2::labs(x = "Correlation coefficient", y = "DEG count", col="")+
  ggplot2::geom_vline(xintercept = 0.4, linetype = 'dashed', colour="darkred")+
  theme_bw() 
  

#gene_p_values<- sapply(rownames(cnv_filt), function(gene) {
  #calculate_p_value(cnv_filt[gene, ], rna_tumor[gene, ])
#})

#cor_results$p_value <- gene_p_values
#cor_results$adjusted_p_value <- p.adjust(cor_results$p_value, method = "BH")

#pvalue_threshold <- 0.05
#correlation_threshold <- 0.3

#significant_genes <- cor_results%>%
  #filter(adjusted_p_value < pvalue_threshold & Correlation > correlation_threshold)


## Identify genes with linear dosage dependance (linear model fit) ##
rna_log_normalized <- rna %>% as.matrix %>% DESeq2::vst()
rna_tumor <- rna_log_normalized[,12:22]

common_genes <- intersect(rownames(cnv_filt), rownames(rna_tumor))
cnv_filt <- cnv_filt[common_genes, ]
rna_tumor <- rna_tumor[common_genes, ]

model_results <- t(sapply(rownames(cnv_filt), function(gene) {
  fit_linear_model(cnv_filt[gene, ], rna_tumor[gene, ])
}))

results_df <- data.frame(Gene = rownames(model_results), Slope = model_results[, 1], p.value = model_results[, 2])
results_df$adj.p.value <- p.adjust(results_df$p.value, method = "BH")

# Filter for Linear Dosage Dependence #
pvalue_threshold <- 0.05
slope_threshold <- 0.3

dosage_dependent_genes <- results_df %>%
  filter(adj.p.value < pvalue_threshold & abs(Slope) > slope_threshold)

dosage_independent_genes <- results_df %>%
  filter(adj.p.value >= pvalue_threshold | abs(Slope) <= slope_threshold)

# Subset both datasets to common genes #
common_genes <- intersect(rownames(cnv_filt), rownames(rna_filt))
cnv_filt <- cnv_filt[common_genes, ]
rna_filt <- rna_filt[common_genes, ]

cnv_normal <- matrix(2, nrow(rna), 20)
rownames(cnv_normal) <- rownames(cnv_filt)
cnv <- cbind(cnv_normal, cnv_filt)
colnames(cnv) <- colnames(rna)
cnv <- cnv/2
colnames(rna) <- paste0("sample", 1:(ncol(rna)))
colnames(cnv) <- colnames(rna)

#Generate metadata#
metadata <- data.frame(patID = colnames(rna),
                       condition = rep(c("A", "B"), each = 20))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 
metadata$condition <- as.factor(metadata$condition)

#saveRDS(cnv_filt, file = "TCGA/lung_cancer/LUAD/cnv_filt.RDS")
write.csv(cnv, file = "TCGA/lung_cancer/LUAD/cnv_test_4.csv", row.names = T)
write.csv(rna, file = "TCGA/lung_cancer/LUAD/rna_test_4.csv", row.names = T)
write.csv(metadata, file = "TCGA/lung_cancer/LUAD/metadata_4.csv", row.names = T)


### Explortive analysis RNA (z-score) vs CNV relationship ###

# Normalization 
rna_log_normalized <- rna %>% as.matrix() %>% DESeq2::rlog()

# Z-score trasformation
rna_zscore <- t(scale(t(rna_log_normalized)))
rna_zscore_normal <- rna_zscore %>% as.data.frame() %>% dplyr::select(1:10)
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% dplyr::select(11:20) %>% as.matrix()

rna_zscore_normal <- rna_zscore_normal %>% 
  as.data.frame() %>%
  dplyr::mutate(zscore_mean = rowMeans(rna_zscore_normal)) %>% 
  dplyr::select(zscore_mean)

rna_zscore_tumor <- rna_zscore_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(zscore_mean = rowMeans(rna_zscore_tumor)) %>% 
  dplyr::select(zscore_mean)

#CNV factorization
#cnv_mean <- cnv_filt %>% 
  #as.data.frame() %>% 
  #mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  #select(cnv_mean)

cnv <- cnv_mean %>% 
  dplyr::mutate(cnv = case_when(
    cnv_mean > 0.5 & cnv_mean <= 1.7 ~ "1",
    cnv_mean > 1.7 & cnv_mean <= 2.5 ~ "2",
    cnv_mean > 2.5 & cnv_mean <= 3.5 ~ "3",
    cnv_mean > 3.3 & cnv_mean <= 4.2 ~ "4",
    cnv_mean > 4.2 ~ "5")) %>% 
  dplyr::select(cnv)

cnv <- cnv[rownames(cnv) %in% rownames(rna_zscore_normal),]

### Preparing data for boxplot ###
rna_zscore_tumor <- rna_zscore_tumor %>% dplyr::mutate(sample_type = "Tumor")
rna_zscore_normal <- rna_zscore_normal %>% dplyr::mutate(sample_type = "Normal")

p_luad_n <- cbind(rna_zscore_normal, cnv)  
p_luad_n$cnv <- as.factor(p_luad_n$cnv)
p_luad_t <- cbind(rna_zscore_tumor, cnv) 
p_luad_t$cnv <- as.factor(p_luad_t$cnv)
plot_all_luad <- rbind(p_luad_n, p_luad_t) %>% dplyr::mutate(cancer_type = "LUAD")

#Apply filtering
plot_all_luad_filt <- plot_all_luad %>% 
  dplyr::filter(cnv=="2" & abs(zscore_mean) <= 0.5 | cnv=="1" | cnv=="3" | cnv=="4" | cnv=="5")
p_luad_t <- p_luad_t %>% 
  dplyr::filter(cnv=="2" & abs(zscore_mean) <= 0.5 | cnv=="1" | cnv=="3" | cnv=="4" | cnv=="5")

#Median zscore
median_z_scores <- p_luad_t %>%
  group_by(cnv) %>%
  summarise(zscore_median = median(zscore_mean, na.rm = TRUE)) %>% 
  as.data.frame()

m_zscore <- ggplot(median_z_scores, aes(x = cnv, y = zscore_median)) +
  geom_point() +  
  geom_smooth(method='lm') + 
  labs(title = "TCGA", x = "Copy Number", y = "Median Z score") +
  theme_bw()+
  scale_y_continuous(breaks = seq(-0.8, 0.8, 0.1))
m_zscore

#Boxplot #

summary.stats <- p_luad_t %>%
  group_by(cnv) %>%
  get_summary_stats() %>%
  select(cnv, n)
summary.stats

summary.plot <- ggsummarytable(
  summary.stats, x = "cnv", y = c("n"),
  ggtheme = theme_bw()
)
summary.plot

col <- qualitative_hcl(5, palette = "Warm")

bxp_t <- ggplot(p_luad_t, aes(x = cnv, y = zscore_mean, fill = cnv)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch = F, show.legend = F)+
  geom_smooth(method = "loess", formula = y ~ x, se=FALSE, color="darkred", aes(group=1), linetype = 'dashed')+
  labs(x="CN group", y = "mRNA Z-score")+
  #facet_wrap(~cancer_type)+
  #theme(strip.text.x = element_text(size=12, color="black", face="bold.italic"))+
  ggplot2::theme(legend.position = 'none')+
  theme_bw()+
  scale_fill_manual(values=col)+
  font("xy.text", size = 12, color = "black", face = "bold")+
  font("title", size = 12, color = "black")+
  font("xlab", size = 12)+
  font("ylab", size = 12)+
  ggtitle("TCGA - LUAD")+
  theme(plot.title = element_text(hjust = 0.5))
bxp_t

#Comparison boxplot Tumor vs Normal
bxp_all <- ggplot(plot_all_luad_filt, aes(x = cnv, y = zscore_mean, fill = sample_type)) + 
  geom_boxplot(position = position_dodge())+
  labs(x="CN group", y = "mRNA Z-score", title = "LUAD")+
  theme_bw()+
  ggplot2::theme(legend.position = 'bottom')+
  #facet_wrap(~cancer_type)
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  font("xy.text", size = 12, color = "black", face = "bold")+
  font("title", size = 12, color = "black")+
  font("xlab", size = 12)+
  font("ylab", size = 12)+
  ggtitle("TCGA - LUAD")+
  theme(plot.title = element_text(hjust = 0.5))
bxp_all

#ggarrange(bxp, summary.plot, ncol = 1, align = "v", heights = c(0.80, 0.20))
gridExtra::grid.arrange(bxp_t, bxp_all, nrow = 1)

#cnv <- cnv %>% 
  #mutate(cn_group = case_when(
    #cnv == "2" ~ "diploid",
    #cnv == "0" ~ "cn_loss",
    #cnv == "1" ~ "cn_loss",
    #cnv == "3" ~ "cn_gain",
    #cnv == "4" ~ "cn_gain",
    #cnv == "5" ~ "cn_amplification"))

#gene group factorization
#res_noCNV <- res_noCNV %>% 
  #mutate(gene_group = case_when(
    #B1_1 < -0.6 & padj < 0.05 ~ "DEG",
    #B1_1 > 0.6 & padj < 0.05 ~ "DEG",
    #B1_1 <= 0.6 & B1_1 >= -0.6 ~ "noDEG",
    #B1_1 >= 0.6 & padj > 0.05 ~ "not significant",
    #B1_1 <= -0.6 & padj > 0.05 ~ "not significant"))

#Gene group facrorization based on Effect size difference
#deg <- deg %>% 
  #mutate(gene_group = case_when(
    #Difference <= -1.0 ~ "super-dosage",
    #Difference >= 1.0 ~ "super-dosage",
    #Difference > -1.0 & Difference < -0.4 ~ "dosage-sensitive",
    #Difference < 1.0 & Difference > 0.4 ~ "dosage-sensitive",
    #Difference >= -0.4 & Difference <= 0.4 ~ "dosage_insensitive"
  #))

#Calculate SD across rows and columns 
#library(matrixStats)

#rna <- transform(rna, sd=apply(rna, 1, sd, na.rm=TRUE)) %>% 
  #subset(sd > 25.0) %>% 
  #select(1:192) %>% 
  #as.matrix()

#cnv_tumor <- as.matrix(cnv_tum)
#cnv_tumor_sd <- matrixStats::colSds(cnv_tumor) %>% 
  #as.data.frame() %>% 
  #setNames("sd")  %>% 
  #subset(sd > 1)