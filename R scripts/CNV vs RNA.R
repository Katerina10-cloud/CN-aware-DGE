# Relationship CNV vs RNA

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(DESeq2)

#Exploration of CNV and RNAseq data

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
cnv_tumor <- cnv_tumor %>% replace(cnv_tumor, cnv_tumor>4, 5) %>% 
  mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  subset(cnv_tumor, cnv_mean <= 0.9 | cnv_mean >= 1.5) %>% 
  select(1:10) %>%
  remove_rownames %>% column_to_rownames(var="Row.names")

cnv_tumor <- cnv_tumor %>% mutate(cnv_type = ifelse(cnv_mean == 2,"neutral", ifelse(cnv_mean>2, "gain", "loss")))

cnv_normal <- data.frame(s1_normal = rep(1, nrow(cnv)), s2_normal = rep(1, nrow(cnv)), 
                         s3_normal = rep(1, nrow(cnv)), s4_normal = rep(1, nrow(cnv)), 
                         s5_normal = rep(1, nrow(cnv)), s6_normal = rep(1, nrow(cnv)), 
                         s7_normal = rep(1, nrow(cnv)), s8_normal = rep(1, nrow(cnv)), 
                         s9_normal = rep(1, nrow(cnv)), s10_normal = rep(1, nrow(cnv))) 

#transform log2 ratio to integer
#round( (2^1.5) * 2)
log2fc_integer <- function(x){ 
  round((2^x)*2)
}
cnv_tumor[1:10] <- lapply(cnv_tumor[1:10], FUN = log2fc_integer)

#Metadata generation
metadata <- data.frame(patID = c("s1_normal", "s2_normal", "s3_normal", "s4_normal", "s5_normal",
                                 "s6_normal", "s7_normal", "s8_normal", "s9_normal", "s10_normal",
                                 "s1_tumor", "s2_tumor", "s3_tumor", "s4_tumor", "s5_tumor",
                                 "s6_tumor", "s7_tumor", "s8_tumor", "s9_tumor", "s10_tumor"), 
                       condition = rep(c("A", "B"), each = 10)) 
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="patID") %>% 
  as.factor(metadata$condition)


#Correction of normal RNAseq counts for CNV
cnv <- cbind(cnv_tumor, cnv_normal)
cnv_tumor <- cnv[1:10] %>% 
  cnv_tumor/2

#cnv_normal <- cnv[11:20]
rna_normal_tumor <- rna_normal_tumor[(rownames(rna_normal_tumor) %in% rownames(cnv)),] #delete rows by name

cnv <- cnv + 10e-9
rna_normal_tumor <- rna_normal_tumor * cnv


#Making barplot
#row.names(cnv) <- 1:nrow(cnv)
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
rna_lusc <- as.matrix(rna_lusc)
rna.vst <- DESeq2::varianceStabilizingTransformation(rna_lusc)
#rna.log <- DESeq2::rlog(rna_tumor)

rna_tumor <- rna_tumor[(rownames(rna_tumor) %in% rownames(cnv)),]
rna <- rna.vst %>% as.data.frame() %>% select(11:20) 
rna_tumor <- rna_filt %>% as.data.frame() %>% select(11:20) 
rna_normal <- rna_filt %>% as.data.frame() %>% select(1:10)

# Manually scaling
(x - mean(x)) / sd(x)
#z-score calculation Gene Expression
#dim(rna.vst)
#rna_zscore <- t(scale(t(rna.vst)))

rna_tumor <- rna_tumor %>% as.data.frame %>% mutate(rna_mean = rowMeans(rna_tumor)) %>% select(11)

#Selecting most variable genes
nTop = 20000
sds <- genefilter::rowSds(rna.vst)
rna_filt <- rna.vst[order(sds, decreasing = T)[1:nTop],]

#cnv factorization
cnv <- replace(cnv, cnv>10, 10)
cnv <- cnv %>% select(1:10)
cnv <- cnv[(rownames(cnv) %in% rownames(rna_tumor)),]

cnv <- cnv %>%
  #select(1:10) %>% 
  mutate(cnv_mean = rowMeans(cnv)) %>% 
  select(11)

#cnv factorization
cnv <- cnv %>% 
  mutate(group = case_when(
  cnv_mean <= 1.6 ~ "1",
  cnv_mean > 1.6 & cnv_mean <= 2.5 ~ "2",
  cnv_mean > 2.5 & cnv_mean <= 4 ~ "3",
  cnv_mean > 4 ~ "4")) %>% 
  #as.factor(cnv$group) %>% 
  select(2)

plot_data <- cbind(rna_tumor, cnv)

#Boxplot
# Compute summary statistics
summary.stats <- plot_data %>%
  group_by(group) %>%
  get_summary_stats() %>%
  select(group, n)

summary.plot <- ggsummarytable(
  summary.stats, x = "group", y = c("n"),
  ggtheme = theme_bw()
)
summary.plot

#Create boxplot
bxp <- ggplot(plot_data, aes(x = group, y = rna_mean, fill = group)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch = TRUE)+
  labs(title="CNV patterns and mRNA expression (LUSC, tumor samples)",x="CNV group", y = "Normalized gene expression")+
  theme_classic()

ggarrange(
  bxp, summary.plot, ncol = 1, align = "v",
  heights = c(0.80, 0.20)
)


