###--------------------------------------------------------###
### Plots ###
###--------------------------------------------------------###

library(ggplot2)
library(ggpubr)
library(tidyverse)

### Making barplot ###
#row.names(cnv) <- 1:nrow(cnv) cnv[1:3]
#cvn$GeneID <- rownames(cvn)
df <- data.frame(dge_groups=rep(c("DEG", "DEG_CNV"), each=3),
                 gene_groups=rep(c("up_gains", "down_loss", "neutral_cnv"),2),
                 frequency=c(4.6, 13.0, 60.6, 13.2, 3.7, 8.4))

plot <- ggplot(data=df, aes(x=gene_groups, y=frequency, fill=dge_groups)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=frequency, label=frequency), vjust=1.6, 
            color="black", size=3.5)+
  labs(title="DGE and CNV relationship")+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()


#preparing data for boxplot
rna_zscore_tumor <- rna_zscore_tumor %>% mutate(sample_type = "Tumor")
rna_zscore_normal <- rna_zscore_normal %>% mutate(sample_type = "Normal")

plot_data_1 <- cbind(rna_zscore_normal, cnv)  
plot_data_1$cnv <- as.factor(plot_data_1$cnv)

plot_data_2 <- cbind(rna_zscore_tumor, cnv) 
plot_data_2$cnv <- as.factor(plot_data_2$cnv)
plot_data <- rbind(plot_data_1, plot_data_2)


### Boxplot ###

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
# The palette with grey:
cbPalette <- c("#0072B2", "#D55E00", "#999999", "#E69F00", "#CC79A7", "#D55E00")

bxp1 <- ggplot(plot_data, aes(x = cnv, y = rna_mean, fill = cnv)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch = TRUE)+
  labs(title="CNV patterns and mRNA expression (LUAD)",x="CNV group", y = "mRNA Z-score")+
  theme_classic()
bxp1 <- bxp + scale_fill_manual(values=cbPalette)+
  font("xy.text", size = 18, color = "black", face = "bold")+
  font("title", size = 18, color = "black", face = "bold.italic")+
  font("xlab", size = 16)+
  font("ylab", size = 16)
bxp1

#Comparison boxplot Tumor vs Normal
bxp2 <- ggplot(plot_data, aes(x = cnv, y = rna_mean, fill = sample_type)) + 
  geom_boxplot(position = position_dodge())+
  labs(title="CNV patterns and mRNA expression (LUAD)",x="CNV group", y = "mRNA Z-score")+
  #theme_classic()+
  facet_wrap(~sample_type, ncol=5)
bxp2

ggarrange(
  bxp, summary.plot, ncol = 1, align = "v",
  heights = c(0.80, 0.20)
)

###-------------------------------------------------------------------###
### Violin plot ###
###-------------------------------------------------------------------###
  
plot_1 <- ggplot(deg, aes(x = cnv, y = difference, fill = cnv))+
  geom_violin(trim=FALSE)+
  #geom_jitter(shape=10, position=position_jitter(0.1))+
  labs(title="Copy-number-informed DEG (n=7522)",x="CN group", y = "Effect size difference (log2)")+
  geom_hline(yintercept = 0, linetype='dashed', color='blue')+
  geom_boxplot(width=0.1)+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 12, color = "black")+
  font("xlab", size = 12)+
  font("ylab", size = 12)

plot_2 <- ggplot(genes_cnvReg, aes(x = cnv, y = difference, fill = cnv))+
  geom_violin(trim=FALSE)+
  #geom_jitter(shape=10, position=position_jitter(0.1))+
  labs(title="Copy-number regulated genes (n=1309)",x="CN group", y = "Effect size difference (log2)")+
  geom_hline(yintercept = 0, linetype='dashed', color='blue')+
  geom_boxplot(width=0.1)+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 12, color = "black")+
  font("xlab", size = 12)+
  font("ylab", size = 12)

grid.arrange(barplot, plot_1, plot_2, nrow = 1)


###-------------------------------------------------------------------###
### Scatter plot ###
###-------------------------------------------------------------------###

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

barplot <- ggplot(data_barplot, aes(y = number_of_genes, x = gene_group, fill = cn_group))+
  geom_bar(stat = "identity")+
  labs(x='Gene group', y='Frequency', title='CN-informed Gene Expression, LUAD (n=23 080)')+
  scale_fill_brewer(palette="")+
  #geom_text(aes(gene_group, label = number_of_genes), size = 3, position=position_dodge2(width=0.5))+
  theme_bw()+
  #facet_wrap("cn_group")+
  guides(fill=guide_legend("CN group"))+
  font("xy.text", size = 12, color = "black")+
  font("xlab", size = 12)+
  font("ylab", size = 12)
barplot


###-------------------------------------------------------------------###
### Volcano Plot ###
###-------------------------------------------------------------------###

library(gridExtra)
library(tidyverse)

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
#colnames(res3_nocnv)[3] <- "B1_1"
#colnames(res4_cnv)[3] <- "B1_2"
res1$diffexpressed <- "NO"
res1$diffexpressed[res1$log2FoldChange > 0.5 & res1$padj < 0.05] <- "UP"
res1$diffexpressed[res4$log2FoldChange < -0.5 & res1$padj < 0.05] <- "DOWN"

#Make simple graphics
p1 <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "black", "red"))+
  scale_x_continuous(breaks = seq(-8, 8, 1))+
  labs(title="DESeq2: CN signal",x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 10, color = "black")+
  font("xlab", size = 10)+
  font("ylab", size = 10)+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p1

p2 <- ggplot(data = res2, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("black", "blue", "red"))+
  scale_x_continuous(breaks = seq(-8, 8, 1))+
  labs(title="DESeqCN: CN signal",x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 10, color = "black")+
  font("xlab", size = 10)+
  font("ylab", size = 10)+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p2

p3 <- ggplot(data = res3, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "black", "red"))+
  scale_x_continuous(breaks = seq(-10, 10, 2))+
  labs(title="DESeq2: mixed signals", x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 10, color = "black")+
  font("xlab", size = 10)+
  font("ylab", size = 10)+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p3

p4 <- ggplot(data = res4, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "black", "red"))+
  scale_x_continuous(breaks = seq(-8, 8, 1))+
  labs(title="DESeqCN: mixed signals", x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 10, color = "black")+
  font("xlab", size = 10)+
  font("ylab", size = 10)+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p4

res_nocnv <- stat_res_luad %>% select(B1_1, padj_1) %>% rename("padj" = "padj_1")
res_cnv <- stat_res_luad %>% select(B1_2, padj_2) %>% rename("padj" = "padj_2")
  
p5 <- ggplot(data = res_map_tum_norm, aes(x = B1_1, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "black", "red"))+
  scale_x_continuous(breaks = seq(-14, 10, 4))+
  labs(x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 8, color = "black")+
  font("xlab", size = 8)+
  font("ylab", size = 8)
p5

#delete rows by name
res_nocnv <- res_nocnv[!(row.names(res_nocnv) %in% c("PYCR1")),]

p6 <- ggplot(data = res_cnv, aes(x = B1_2, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "black", "red"))+
  scale_x_continuous(breaks = seq(-14, 10, 4))+
  labs(x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 8, color = "black")+
  font("xlab", size = 8)+
  font("ylab", size = 8)
p6

#Plots
gridExtra::grid.arrange(p1, p2, p3,p4, nrow = 2)
#grid.arrange(g2, arrangeGrob(g3, g4, ncol=2), nrow = 2)


###------------------------------------------------------------###
### Histogram ###
###------------------------------------------------------------###

plot2 <- ggplot(res2, aes(x=padj))+
  geom_histogram(color="darkblue", fill="lightblue", bins = 100)+
  geom_vline(xintercept = c(0.05), col = "red", linetype = 'dashed')+
  labs(title="DESeqCN : CN signal", x="FDR")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
plot1

    
###----------------------------------------------------------------###
### Density plot ###
###----------------------------------------------------------------###  

ggplot(resFit_merged, aes(Difference)) +
  geom_histogram(bins = 1000) +
  xlim(-4,4) +
  labs(title="B1 difference") +
  theme_classic()

ggplot(resFit_plot,aes(x=B1, fill = Group)) + 
  geom_histogram(data=subset(resFit_plot,Group == 'B1_1'),alpha =0.6, bins = 1000) +
  geom_histogram(data=subset(resFit_plot,Group == 'B1_2'), alpha =0.6, bins = 1000) +
  xlim(-7,12) +
  labs(title="B1 parameters density (MAP) test 2") +
  theme(legend.position="bottom", legend.box = "horizontal") + 
  theme_classic()

#Scatter plot
# scatter NV vs DP
scatter_plot = ggplot(statRes_map_CNV, aes(x = cnv_mean, y = Log2FC)) +
  geom_point(aes(color = cnv_type), size = .8) +
  theme(legend.position = 'bottom') +
  labs(title = "CNV correcrtion") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_hline(yintercept=1, linetype="dashed", color = "blue") +
  geom_hline(yintercept=-1, linetype="dashed", color = "blue") +
  geom_vline(xintercept = 1.8, color = "blue", size = 0.4) +
  geom_vline(xintercept = 2.5, color = "blue", size = 0.4) +
  geom_vline(xintercept = 4, color = "blue", size = 0.4) +
  theme_classic()
#facet_wrap(~FILTER) +
scatter_plot


###-------------------------------------------------------------###
### ROC curve ###
###-------------------------------------------------------------###

#BiocManager::install("metaseqR2")
library(metaseqR2)

p1 <- matrix(res3$padj)
colnames(p1) <- "DESeq2"
p2 <- matrix(res4$padj)
colnames(p2) <- "DESeqCN"
p <- cbind(p1,p2)

res3 <- res3 %>% 
  mutate(truth = case_when(
    log2FoldChange >= 1.0 & padj < 0.05 ~ "1",
    log2FoldChange <= -1.0 & padj < 0.05 ~ "1",
    log2FoldChange < 1 | log2FoldChange > -1 & padj >= 0.05 ~ "0")) 

truth <- as.vector(as.numeric(res4$truth))
names(truth) <- res4$X

metaseqR2::diagplotRoc(truth = truth, p = p2, sig = 0.05, x = "fpr",
                       y = "tpr", path = NULL, draw = TRUE)




