# Relationship CNV vs RNA

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(DESeq2)

#Exploration of CNV and RNAseq data 
statRes_map_CNV = read.csv('~/model_fit_Python/model_results/results_luad/results_3/statRes_map_CNV.csv',header=TRUE)
statRes_map_noCNV = read.csv('~/model_fit_Python/model_results/results_brca/statRes_map_noCNV.csv',header=TRUE)
#deg_merged$group <- as.factor(deg_merged$group)

save(plot_data_2, file = "~/model_fit_Python/model_results/results_luad/results_3/plot_data_2.Rdata")
save(metadata, file = "~/model_data/TCGA/lung_cancer/LUAD/metadata.Rdata")
write.csv(rna, file = "~/model_data/TCGA/lung_cancer/LUAD/rna_cnv_10.csv")

statRes_map_noCNV <- statRes_map_noCNV %>% remove_rownames %>% column_to_rownames(var="GeneID")
rna_counts <- rna_counts[(rownames(rna_counts) %in% rownames(res_allGenes)),]
statRes_map_noCNV <- cbind(statRes_map_noCNV, cnv)


#res_allGenes$GeneID <- rownames(res_allGenes)
#res_allGenes <- res_allGenes %>% mutate(difference = B1_2 - B1_1)
#cnv <- cnv[(rownames(cnv) %in% rownames(res_allGenes)),]
#cnv <- cnv %>% mutate(cnv_mean = rowMeans(cnv))

#rna_tumor = read.delim("model_data/TCGA/lung_cancer/LUSC/s1/s1_rna_tumor.tsv", header=TRUE, sep="\t")
#rna <- read.delim("~/model_data/TCGA/lung_cancer/LUSC/s10/s10_rna_tumor.tsv", header=TRUE, sep="\t")

#cnv_brca <- brca_cnv_tumor %>% select(7:10, 13, 14, 18, 19, 22, 23, 25:27, 29, 34, 35, 37, 41, 45, 47, 49, 53, 54, 56:58, 61, 62, 64, 69, 72, 73, 76, 83, 90, 92:96, 102, 107, 109)
#cnv_3 <- cnv_brca %>% select(23,24,36)
#rna_norm_3 <- brca_rna_norm %>% select(55,57,93)
#rna_tum_3 <- brca_rna_tum %>% select(55,57,92)

statRes_map_cnv <- statRes_map_cnv %>% select(1,3,7)
statRes_map_NOcnv <- na.omit(statRes_map_NOcnv)
colnames(statRes_map_NOcnv)[2] <- "padj_1"

cnv_3$GeneID <- row.names(cnv_3)
cnv_3 <- merge(statRes_map_NOcnv, cnv_3, by="GeneID")
cnv_3 <- cnv_3 %>% select(1, 4:6)
rna_cnv <- rna_cnv %>% remove_rownames %>% column_to_rownames(var="Row.names")

rna_cnv <- merge(rna_cnv, cnv, by = "row.names")


rna_tumor <- rna_tumor %>% select(2,4) %>% na.omit() %>% 
  colnames(rna_tumor)[2] <- "s1_tumor" %>% 
  rna_tumor[!duplicated(rna_tumor$GeneID), ] %>% #remove dublicates 
  remove_rownames %>% column_to_rownames(var="GeneID") 

rna_normal <- rna_normal %>% select(2,4) %>% na.omit() %>% 
  colnames(rna_normal)[2] <- "s1_tumor" %>% 
  rna_normal[!duplicated(rna_normal$GeneID), ] %>% #remove dublicates
  remove_rownames %>% column_to_rownames(var="GeneID")

rna_normal_tumor <- cbind(rna_normal, rna_tumor)

rna_cnv = read.csv('~/model_fit_Python/model_data/last_test/rna_cnv.csv',header=TRUE)
#colnames(cnv)[9] <- "cnv_mean"

#cnv <- replace(cnv, cnv>5,6)  
cnv_3_luad <- cnv_3_luad %>% mutate(cnv_mean = rowMeans(cnv_3_luad)) 


cnv_tumor <- luad_cnv_tumor %>% replace(cnv_tumor, cnv_tumor>5, 6) %>% 
  mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  subset(cnv_tumor, cnv_mean <= 0.9 | cnv_mean >= 1.5) %>% 
  select(1:10) %>%
  remove_rownames %>% column_to_rownames(var="Row.names")


cnv_normal <- data.frame("TCGA-50-5932-11A" = rep(1, nrow(luad_cnv)), "TCGA-44-6147-11A" = rep(1, nrow(luad_cnv)),
                         "TCGA-55-6979-11A" = rep(1, nrow(luad_cnv)), "TCGA-50-5931-11A" = rep(1, nrow(luad_cnv)),
                         "TCGA-91-6835-11A" = rep(1, nrow(luad_cnv)), "TCGA-44-6776-11A" = rep(1, nrow(luad_cnv)),
                         "TCGA-44-6778-11A" = rep(1, nrow(luad_cnv)), "TCGA-44-6145-11A" = rep(1, nrow(luad_cnv)),
                         "TCGA-50-5939-11A" = rep(1, nrow(luad_cnv)), "TCGA-50-5935-11A" = rep(1, nrow(luad_cnv)))
cnv <- cbind(luad_cnv, cnv_normal)

#transform log2 ratio to integer
#round( (2^1.5) * 2)
log2fc_integer <- function(x){ 
  round((2^x)*2)
}
cnv_tumor[1:10] <- lapply(cnv_tumor[1:10], FUN = log2fc_integer)

#Metadata generation
metadata <- data.frame(patID = colnames(luad_rna), 
                       condition = rep(c("A", "B"), each = 10)) 
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="patID")  
metadata$condition <- as.factor(metadata$condition)

#Correction of normal RNAseq counts for CNV
cnv <- cbind(cnv, cnv_normal_3)
cnv <- cnv[1:3]/2
cnv_3 <- cnv_3/2

#cnv_normal <- cnv[11:20]
#cnv <- cnv %>% cnv[(rownames(cnv %in% rownames(res_allGenes)),] #delete rows by name

cnv <- cnv + 10e-9
rna <- luad_rna * cnv


#Making barplot
#row.names(cnv) <- 1:nrow(cnv) cnv[1:3]
#cvn$GeneID <- rownames(cvn)
df <- data.frame(dge_groups=rep(c("DEG", "DEG_CNV"), each=3),
                 gene_groups=rep(c("up_gains", "down_loss", "neutral_cnv"),2),
                 frequency=c(4.6, 13.0, 60.6, 13.2, 3.7, 8.4))

plot_1 <- ggplot(data=df, aes(x=gene_groups, y=frequency, fill=dge_groups)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=frequency, label=frequency), vjust=1.6, 
            color="black", size=3.5)+
  labs(title="DGE and CNV relationship")+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

#Plots

#Counts normalization
load("~/model_data/TCGA/lung_cancer/LUSC/cnv_lusc.Rdata")
load("~/model_data/TCGA/lung_cancer/LUSC/rna_lusc.Rdata")

cnv <- cnv/2
cnv <- cnv + 10e-9
rna_normal <- rna_normal * cnv
#rna_normal <- round(rna_normal, 0)

#Counts normalization
rna_normalized <- luad_rna %>%  as.matrix()
#rna_normalized <- DESeq2::varianceStabilizingTransformation(rna_normalized)
rna_log_normalized <- DESeq2::rlog(rna_normalized)

#rna.vst <- rna.vst[(rownames(rna.vst) %in% rownames(cnv)),]
#rna <- rna.vst %>% as.data.frame() %>% select(11:20) 
#rna_tumor <- rna.vst %>% as.data.frame() %>% select(4:6) 
#rna_normal <- rna.vst %>% as.data.frame() %>% select(1:3)

# Manually scaling
#(x - mean(x)) / sd(x)
#z-score calculation Gene Expression
#dim(rna.vst)

rna_tum_norm <- cbind(rna_tum, rna_norm)

rna_zscore <- t(scale(t(rna_log_normalized)))

rna_zscore_normal <- rna_zscore %>% as.data.frame() %>% select(1:10)
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% select(11:20) 
rna_zscore_normal <- rna_zscore_normal %>% as.data.frame %>% mutate(rna_mean = rowMeans(rna_zscore_normal)) %>% select(11)
rna_zscore_tumor <- rna_zscore_tumor %>% as.data.frame %>% mutate(rna_mean = rowMeans(rna_zscore_tumor)) %>% select(11)


#Selecting most variable genes
#nTop = 10000
#sds <- genefilter::rowSds(rna.vst)
#rna_filt <- rna.vst[order(sds, decreasing = T)[1:nTop],]

#cnv factorization
#cnv_tumor_3 <- replace(cnv_tumor_3, cnv_tumor_3>5, 6)
#cnv <- cnv %>% select(5,6,7)
#cnv <- cnv[(rownames(cnv) %in% rownames(rna_filt)),]

luad_cnv <- luad_cnv %>%
  mutate(cnv_mean = rowMeans(luad_cnv))

cnv <- luad_cnv %>% select(11)

#CNV factorization
cnv <- cnv %>% 
  mutate(cnv = case_when(
  cnv_mean <= 0.5 ~ "0",
  cnv_mean > 0.5 & cnv_mean <= 1.6 ~ "1",
  cnv_mean > 1.6 & cnv_mean <= 2.4 ~ "2",
  cnv_mean > 2.4 & cnv_mean <= 3.5 ~ "3",
  cnv_mean > 3.5 & cnv_mean <= 4.5 ~ "4",
  cnv_mean > 4.5 ~ "5")) 

cnv <- cnv %>% 
  mutate(cn_group = case_when(
  cnv == "2" ~ "diploid",
  cnv == "0" ~ "cn_loss",
  cnv == "1" ~ "cn_loss",
  cnv == "3" ~ "cn_gain",
  cnv == "4" ~ "cn_gain",
  cnv == "5" ~ "cn_amplification"))

#gene group factorization
res_noCNV <- res_noCNV %>% 
  mutate(gene_group = case_when(
    B1_1 < -0.6 & padj < 0.05 ~ "DEG",
    B1_1 > 0.6 & padj < 0.05 ~ "DEG",
    B1_1 <= 0.6 & B1_1 >= -0.6 ~ "noDEG",
    B1_1 >= 0.6 & padj > 0.05 ~ "not significant",
    B1_1 <= -0.6 & padj > 0.05 ~ "not significant"))

#Gene group facrorization based on Effect size difference
deg <- deg %>% 
  mutate(gene_group = case_when(
    Difference <= -1.0 ~ "super-dosage",
    Difference >= 1.0 ~ "super-dosage",
    Difference > -1.0 & Difference < -0.4 ~ "dosage-sensitive",
    Difference < 1.0 & Difference > 0.4 ~ "dosage-sensitive",
    Difference >= -0.4 & Difference <= 0.4 ~ "dosage_insensitive"
  ))


cnv <- res_allGenes %>% select(5)
deg_b1 <- deg_merged %>% select(1,3)
deg_b2 <- deg_merged %>% select(2,3)
plot_data_1 <- deg_b1 %>% mutate(effect_size = "B1_1")
plot_data_2 <- deg_b2 %>% mutate(effect_size = "B1_2")
plot_data <- rbind(plot_data_1, plot_data_2)

colnames(plot_data_2)[1] <- "B1"


rna_zscore_tumor <- rna_zscore_tumor %>% mutate(sample_type = "Tumor")
rna_zscore_normal <- rna_zscore_normal %>% mutate(sample_type = "Normal")


plot_data_2 <- merge(rna_zscore_tumor, cnv, by = "row.names")
plot_data_2 <- plot_data_2 %>% remove_rownames %>% column_to_rownames(var="Row.names")
plot_data_2$cnv <- as.factor(plot_data_2$cnv)
plot_data <- rbind(plot_data_2, plot_data_1)


#Boxplot
# Compute summary statistics
summary.stats <- plot_data_2 %>%
  group_by(cnv) %>%
  get_summary_stats() %>%
  select(cnv, n)
summary.stats <- summary.stats[c(1,3,5,7,9),]

summary.plot <- ggsummarytable(
  summary.stats, x = "cnv", y = c("n"),
  ggtheme = theme_bw()
)
summary.plot



#Create boxplot
bxp <- ggplot(plot_data_2, aes(x = cnv, y = rna_mean, fill = cnv)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch = TRUE)+
  labs(title="CNV patterns and mRNA expression (LUAD)",x="CNV group", y = "mRNA Z-score")+
  theme_classic()
bxp

bxp <- ggplot(plot_data, aes(x = cnv, y = rna_mean, fill = sample_type)) + 
  geom_boxplot(position = position_dodge())+
  labs(title="CNV patterns and mRNA expression (BRCA, tumor samples (3))",x="CNV group", y = "mRNA Z-score")+
  theme_classic()+
  facet_wrap(~sample_type, ncol=6)
bxp

ggarrange(
  bxp, summary.plot, ncol = 1, align = "v",
  heights = c(0.80, 0.20)
)

#Comparison boxplot
ggplot(plot_data, aes(x = effect_size, y = B1, fill = effect_size)) + 
  geom_boxplot(position = position_dodge()) +
  labs(title="CNV patterns and effect size (DEG (3482), Tumor vs Normal)",x="CNV group", y = "B1")+
  facet_wrap(~group, ncol=6) +
  theme_classic()

#Comparison boxplot
ggplot(plot_data, aes(x = sample_type, y = rna_mean, fill = sample_type)) + 
  geom_boxplot(position = position_dodge()) +
  labs(title="CNV patterns and mRNA expression (LUAD,Tumor vs Normal)",x="CNV group", y = "mRNA Z-score")+
  facet_wrap(~cnv, ncol=6) +
  theme_classic()

#load("~/model_fit_Python/model_results/lusc_fit/")

#Violin plot
violin_plot <- ggplot(res_allGenes, aes(x = group, y = difference, fill = group))+
  geom_violin(trim=FALSE)+
  #geom_jitter(shape=10, position=position_jitter(0.1))+
  labs(title="CNV patterns and Effect size difference ((|B1_2| - |B1_1|), (n = 27 489 genes))",x="CNV group", y = "Effect size difference (log2FC)")+
  geom_hline(yintercept = 0, linetype='dashed', color='blue')+
  geom_boxplot(width=0.1)+
  theme_classic()
violin_plot  

#Scatter plot
scatterplot <- ggplot(res_allGenes, aes(x=cnv_mean, y=difference)) + 
  geom_point()+
  geom_smooth()+
  labs(title="CNV patterns and Effect size difference relationship ((|B1_2| - |B1_1|), (n = 27489 genes))",x="CNV mean", y = "Effect size difference")+
  theme_classic()+
  geom_hline(yintercept = 0, linetype='dashed', color='red')+
  geom_vline(xintercept = 2, linetype='dashed', color='blue')
scatterplot  

#Stacked Barplot
#sum(res_allGenes$gene_group == "other" & res_allGenes$cn_group == "Diploid")
data_barplot <- data.frame(
  gene_group = rep(c("DEG", "no_DEG", "not significant"), each = 4),
  cn_group = rep(c("cn_loss", "diploid", "cn_gain", "cn_amplificatin"), 3),
  number_of_genes = c(314, 3977, 3173, 58, 528, 7862, 4981, 46, 85, 1189, 860, 7)
)

data_barplot_geneDosage <- data.frame(
  gene_group = rep(c("super-dosage", "d-sensitive", "d-insensitive"), each = 6),
  cn_group = rep(c("0", "1", "2", "3", "4", "5"), 3),
  number_of_genes = c(8, 26, 30, 16, 33, 139, 0, 34, 544, 203, 204, 13, 1, 0, 1225, 483, 6, 2)
) 

barplot <- ggplot(data_barplot, aes(fill = cn_group, y = number_of_genes, x = gene_group))+
  geom_bar(stat = "identity")+
  labs(x='Gene group', y='Frequency', title='CNV informed Gene Expression, LUAD (n_genes=23 080)')+
  scale_fill_manual('Position', values=c('red', 'pink', 'steelblue', 'jellow'))+
  #geom_text(aes(gene_group, label = number_of_genes), size = 3, position=position_dodge2(width=0.5))+
  theme_minimal()+
  #facet_wrap("cn_group")+
  guides(fill=guide_legend("CN group"))
barplot + scale_fill_brewer(palette = "Set2")
                 