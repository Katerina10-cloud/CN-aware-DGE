### Explortive analysis RNA (z-score) vs CNV relationship ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "cluster", "factoextra", "heatmaply", "DESeq2", "tidyverse", "DESeq2", "colorspace", 
          "ggpubr", "ggpointdensity", "ggeasy")
sapply(pkgs, require, character.only = TRUE)

source("CN-aware-DGE/R/utils.R")

# Input data
data_path <- "TCGA/liver/cnv_tumor.RDS"
dataset_name <- "LIHC_cnv"
cnv_tumor <- readRDS(data_path)

# Clustering patients 
clustering_res <- clustering_patients_cnv(dataset_name, data_path)
clustering_res[[2]]
cnv_tumor <- clustering_res[[1]]

#Select patients for boxplot
cnv_tumor <- as.matrix(t(cnv_tumor))

#LUAD
cnv_tumor <- cnv_tumor %>% as.data.frame() %>% dplyr::select(7,5,13,14,15,18,22,24,25,28) 

#BRCA
#cnv_filt <- cnv_tumor[cnv_tumor$Cluster %in% c(5,4,1),]

#LUSC
#cnv_filt <- cnv_tumor[cnv_tumor$Cluster %in% c(1,2),]

cnv_filt <- cnv_tumor[cnv_tumor$Cluster %in% c(1,2),]
cnv_filt <- subset(cnv_filt, select=-c(Cluster))
cnv_filt <- as.matrix(t(cnv_filt))
cnv_filt <- apply(cnv_filt, 2, function(x) ifelse(x > 15, 15, x)) 

hist(rowMeans(cnv_filt),
     main = "LUSC", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 15)

cnv_mean <- cnv_filt %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  dplyr::select(cnv_mean)

cnv_mean$geneID <- rownames(cnv_mean)


## RNA data ##

data_path <- "TCGA/liver/rna_counts.RDS"
dataset_name <- "LIHC_rna"

rna <- rna_processing(dataset_name, data_path, cnv_filt)
rna_norm <- rna[[1]]
rna_tum <- rna[[2]]

# Exclude genes with low expression in normal tissue #
low_expression_threshold <- 15
expression_summary <- data.frame(
  Gene = rownames(rna_norm),
  MeanExpression = rowMeans(rna_norm)
)

filtered_genes <- expression_summary %>%
  filter(MeanExpression > low_expression_threshold)

rna_norm <- rna_norm[filtered_genes$Gene, ]
rna_tum <- rna_tum[filtered_genes$Gene, ]
rna <- cbind(rna_norm, rna_tum)

# Normalization 
rna_log_normalized <- rna %>% as.matrix() %>% DESeq2::vst()

# Z-score trasformation
rna_zscore <- t(scale(t(rna_log_normalized)))
rna_zscore_normal <- rna_zscore %>% as.data.frame() %>% dplyr::select(1:50)
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% dplyr::select(51:100) %>% as.matrix()

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
    cnv_mean > 0.5 & cnv_mean <= 1.8 ~ "1",
    cnv_mean > 1.8 & cnv_mean <= 2.5 ~ "2",
    cnv_mean > 2.5 & cnv_mean <= 3.5 ~ "3",
    cnv_mean > 3.3 & cnv_mean <= 4.2 ~ "4",
    cnv_mean > 4.2 ~ "5")) %>% 
  dplyr::select(cnv)

cnv <- cnv[rownames(cnv) %in% rownames(rna_zscore_normal),]

### Preparing data for boxplot ###
rna_zscore_tumor <- rna_zscore_tumor %>% dplyr::mutate(sample_type = "Tumor")
rna_zscore_normal <- rna_zscore_normal %>% dplyr::mutate(sample_type = "Normal")

p_lihc_n <- cbind(rna_zscore_normal, cnv)  
p_lihc_n$cnv <- as.factor(p_lihc_n$cnv)
p_lihc_t <- cbind(rna_zscore_tumor, cnv) 
p_lihc_t$cnv <- as.factor(p_lihc_t$cnv)
plot_all_lihc <- rbind(p_lihc_n, p_lihc_t) %>% dplyr::mutate(cancer_type = "LIHC")

#Apply filtering
plot_all_lihc_filt <- plot_all_lihc %>% 
  dplyr::filter(cnv=="2" & abs(zscore_mean) <= 0.5 | cnv=="1" | cnv=="3" | cnv=="4" | cnv=="5")
p_lihc_t <- p_lihc_t %>% 
  dplyr::filter(cnv=="2" & abs(zscore_mean) <= 0.5 | cnv=="1" | cnv=="3" | cnv=="4" | cnv=="5")

p_brca_t <- p_brca_t %>% dplyr::mutate(cancer_type = "BRCA")
p_lusc_t <- p_lusc_t %>% dplyr::mutate(cancer_type = "LUSC")
p_lihc_t <- p_lihc_t %>% dplyr::mutate(cancer_type = "LIHC")

p_tumor <- rbind(p_brca_t, p_lusc_t, p_lihc_t)

#Median zscore
#median_z_scores <- p_luad_t %>%
  #group_by(cnv) %>%
  #summarise(zscore_median = median(zscore_mean, na.rm = TRUE)) %>% 
  #as.data.frame()

#m_zscore <- ggplot(median_z_scores, aes(x = cnv, y = zscore_median)) +
  #geom_point() +  
  #geom_smooth(method='lm') + 
  #labs(title = "TCGA", x = "Copy Number", y = "Median Z score") +
  #theme_bw()+
  #scale_y_continuous(breaks = seq(-0.8, 0.8, 0.1))
#m_zscore

#Boxplot #

summary.stats <- p_lihc_t %>%
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

p_tumor$cnv <- factor(p_tumor$cnv, levels = c(1, 2, 3, 4, 5))

bxp_t <- ggplot(p_tumor, aes(x = cnv, y = zscore_mean, fill = cnv)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch = F, show.legend = F)+
  geom_smooth(method = "loess", formula = y ~ x, se=FALSE, color="darkblue", aes(group=1), linetype = 'dashed')+
  labs(x="CN group", y = "mRNA Z-score")+
  facet_wrap(~cancer_type)+
  theme(strip.text.x = element_text(size=12, color="black", face="bold.italic"))+
  ggplot2::theme(legend.position = 'none')+
  theme_bw()+
  scale_fill_manual(values=col)+
  font("xy.text", size = 12, color = "black", face = "plain")+
  font("title", size = 12, color = "black")+
  font("xlab", size = 12)+
  font("ylab", size = 12)+
  #ggtitle("LIHC")+
  theme(legend.position = "")+
  theme(plot.title = element_text(hjust = 0.5))
bxp_t

plot_all_filt <- rbind(plot_all_brca_filt, plot_all_lihc_filt, plot_all_lusc_filt)
plot_all_filt$cnv <- factor(plot_all_filt$cnv, levels = c(1, 2, 3, 4, 5))

#Comparison boxplot Tumor vs Normal
bxp_all <- ggplot(plot_all_filt, aes(x = cnv, y = zscore_mean, fill = sample_type)) + 
  geom_boxplot(position = position_dodge())+
  labs(x="CN group", y = "mRNA Z-score", title = "")+
  labs(fill = "sample type")+
  theme_bw()+
  ggplot2::theme(legend.position = 'bottom')+
  facet_wrap(~cancer_type)+
  scale_fill_manual(values=c("lightgray", "#D7C78B"))+
  font("xy.text", size = 12, color = "black", face = "plain")+
  font("title", size = 12, color = "black")+
  font("xlab", size = 12)+
  font("ylab", size = 12)+
  #ggtitle("LUSC")+
  theme(plot.title = element_text(hjust = 0.5))
bxp_all

#ggarrange(bxp, summary.plot, ncol = 1, align = "v", heights = c(0.80, 0.20))
gridExtra::grid.arrange(bxp_t, bxp_all, nrow = 2)





