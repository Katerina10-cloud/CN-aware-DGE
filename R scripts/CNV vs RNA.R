# Relationship CNV vs RNA

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(DESeq2)

#Exploration of CNV and RNAseq data 
#rna = read.csv('model_data/TCGA/lung_cancer/last_test/rna_counts_3.csv',header=TRUE)
#deg_merged$group <- as.factor(deg_merged$group)

cnv <- cnv %>% select(4)
deg_merged <- cbind(deg_merged, cnv)

rna_tumor = read.delim("model_data/TCGA/lung_cancer/LUSC/s1/s1_rna_tumor.tsv", header=TRUE, sep="\t")
#rna <- read.delim("~/model_data/TCGA/lung_cancer/LUSC/s10/s10_rna_tumor.tsv", header=TRUE, sep="\t")

rna_tumor <- rna_tumor %>% select(2,4) %>% na.omit() %>% 
  colnames(rna_tumor)[2] <- "s1_tumor" %>% 
  rna_tumor[!duplicated(rna_tumor$GeneID), ] %>% #remove dublicates 
  remove_rownames %>% column_to_rownames(var="GeneID") 

rna_normal <- rna_normal %>% select(2,4) %>% na.omit() %>% 
  colnames(rna_normal)[2] <- "s1_tumor" %>% 
  rna_normal[!duplicated(rna_normal$GeneID), ] %>% #remove dublicates
  remove_rownames %>% column_to_rownames(var="GeneID")

rna_normal_tumor <- cbind(rna_normal, rna_tumor)

#save(rna_normal_tumor, file = "model_data/TCGA/lung_cancer/LUSC/rna_normal_tumor.Rdata")

cnv_tumor = read.csv('model_fit_Python/model_data/test3/cnv.csv',header=TRUE)
#colnames(cnv)[9] <- "cnv_mean"
cnv_tumor <- luad_cnv_tumor 

cnv <- replace(cnv, cnv>5,6)  
cnv <- cnv %>% mutate(cnv_mean = rowMeans(cnv)) 

cnv_tumor <- luad_cnv_tumor %>% replace(cnv_tumor, cnv_tumor>5, 6) %>% 
  mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  subset(cnv_tumor, cnv_mean <= 0.9 | cnv_mean >= 1.5) %>% 
  select(1:10) %>%
  remove_rownames %>% column_to_rownames(var="Row.names")

cnv_tumor <- cnv_tumor %>% mutate(cnv_type = ifelse(cnv_mean == 2,"neutral", ifelse(cnv_mean>2, "gain", "loss")))

cnv_normal <- data.frame(s8_normal = rep(1, nrow(cnv_tumor)), s16_normal = rep(1, nrow(cnv_tumor)), 
                         s27_normal = rep(1, nrow(cnv_tumor)))
cnv <- cbind(cnv_tumor, cnv_normal)

#transform log2 ratio to integer
#round( (2^1.5) * 2)
log2fc_integer <- function(x){ 
  round((2^x)*2)
}
cnv_tumor[1:10] <- lapply(cnv_tumor[1:10], FUN = log2fc_integer)

#Metadata generation
metadata <- data.frame(patID = c("s8_normal", "s16_normal", "s27_normal", 
                                 "s8_tumor", "s16_tumor", "s27_tumor"), 
                       condition = rep(c("A", "B"), each = 3)) 
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="patID") %>% 
  as.factor(metadata$condition)


#Correction of normal RNAseq counts for CNV
cnv <- cbind(cnv, cnv_normal_3)
cnv <- cnv[1:3]/2

#cnv_normal <- cnv[11:20]
#cnv <- cnv %>% cnv[(rownames(cnv %in% rownames(rna)) ,] #delete rows by name

cnv <- cnv + 10e-9
rna_normal_tumor <- rna_normal_tumor * cnv


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
rna_normal <- round(rna_normal, 0)

#Counts normalization
rna_lusc <- rna_lusc %>% select(5,6,7,15,16,17) 
rna_normalized <- rna_normal_tumor_3 %>%  as.matrix()
rna_normalized <- DESeq2::varianceStabilizingTransformation(rna_normalized)
#rna.log <- DESeq2::rlog(rna_tumor)

rna.vst <- rna.vst[(rownames(rna.vst) %in% rownames(cnv)),]
rna <- rna.vst %>% as.data.frame() %>% select(11:20) 
rna_tumor <- rna.vst %>% as.data.frame() %>% select(4:6) 
rna_normal <- rna.vst %>% as.data.frame() %>% select(1:3)

# Manually scaling
(x - mean(x)) / sd(x)
#z-score calculation Gene Expression
#dim(rna.vst)
rna_zscore <- t(scale(t(rna_normalized)))

rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% select(4:6) 
rna_zscore_tumor <- rna_zscore_tumor %>% mutate(rna_mean = rowMeans(rna_zscore_tumor)) %>% na.omit()

rna_zscore_tumor <- rna_zscore_tumor %>% as.data.frame %>% mutate(rna_mean = rowMeans(rna_zscore_tumor)) %>% select(4)
rna_normal <- rna_normal %>% as.data.frame %>% mutate(rna_mean = rowMeans(rna_normal)) %>% select(4)

#Selecting most variable genes
nTop = 10000
sds <- genefilter::rowSds(rna.vst)
rna_filt <- rna.vst[order(sds, decreasing = T)[1:nTop],]

#cnv factorization
cnv_tumor_3 <- replace(cnv_tumor_3, cnv_tumor_3>5, 6)
cnv <- cnv %>% select(5,6,7)
cnv <- cnv[(rownames(cnv) %in% rownames(rna_filt)),]

cnv_tumor_3 <- cnv_tumor_3 %>%
  #select(1:10) %>% 
  mutate(cnv_mean = rowMeans(cnv_tumor_3)) %>% 
  select(4)

#cnv factorization
cnv <- cnv %>% 
  mutate(group = case_when(
  cnv_mean < 0.5 ~ "0",
  cnv_mean >= 0.5 & cnv_mean <= 1.5 ~ "1",
  cnv_mean > 1.5 & cnv_mean < 2.5 ~ "2",
  cnv_mean >= 2.5 & cnv_mean <= 3.5 ~ "3",
  cnv_mean > 3.5 & cnv_mean <= 4.5 ~ "4",
  cnv_mean > 4.5 ~ "5")) %>% 
  select(5)

deg_b1 <- deg_merged %>% select(1,3)
deg_b2 <- deg_merged %>% select(2,3)
plot_data_1 <- deg_b1 %>% mutate(effect_size = "B1_1")
plot_data_2 <- deg_b2 %>% mutate(effect_size = "B1_2")
plot_data <- rbind(plot_data_1, plot_data_2)

colnames(plot_data_2)[1] <- "B1"

#Boxplot
# Compute summary statistics
summary.stats <- plot_data_1 %>%
  group_by(group) %>%
  get_summary_stats() %>%
  select(group, n)

summary.plot <- ggsummarytable(
  summary.stats, x = "group", y = c("n"),
  ggtheme = theme_bw()
)
summary.plot

#Create boxplot
bxp <- ggplot(plot_data_1, aes(x = group, y = rna_mean, fill = group)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch = TRUE)+
  labs(title="CNV patterns and mRNA expression (LUAD, tumor samples (3))",x="CNV group", y = "mRNA Z-score")+
  theme_classic()

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


#load("~/model_fit_Python/model_results/lusc_fit/")

#Violin plot
violin_plot <- ggplot(deg_merged, aes(x = group, y = difference, fill = group))+
  geom_violin(trim=FALSE)+
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")+
  geom_jitter(shape=10, position=position_jitter(0.1))+
  labs(title="CNV patterns and Effect size difference ((|B1_2| - |B1_1|), (n = 3482 DE genes))",x="CNV group", y = "Effect size difference (log2FC)")+
  geom_hline(yintercept = 0, linetype='dashed', color='blue')+
  theme_classic()
violin_plot  

#Scatter plot
scatterplot <- ggplot(deg_merged, aes(x=cnv_mean, y=difference)) + 
  geom_point()+
  geom_smooth()+
  labs(title="CNV patterns and Effect size difference relationship ((|B1_2| - |B1_1|), (n = 3482 DE genes))",x="CNV mean", y = "Effect size difference (log2FC)")+
  theme_classic()
scatterplot  

