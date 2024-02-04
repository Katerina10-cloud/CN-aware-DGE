#Plots
library(ggplot2)
library(ggpubr)

#Making barplot
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

#Violin plot
violin_plot <- ggplot(deg, aes(x = cnv, y = difference, fill = cnv))+
  geom_violin(trim=FALSE)+
  #geom_jitter(shape=10, position=position_jitter(0.1))+
  labs(title="CNV patterns and Effect size difference (|B1_2| - |B1_1|)",x="CNV group", y = "Effect size difference (log2)")+
  geom_hline(yintercept = 0, linetype='dashed', color='blue')+
  geom_boxplot(width=0.1)+
  #theme_classic()
  theme(legend.position="none")+
  font("xy.text", size = 16, color = "black", face = "bold")+
  font("xlab", size = 16)+
  font("ylab", size = 16)
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


# Volcano Plot
library(gridExtra)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
colnames(res_nocnv2)[3] <- "B1_1"
colnames(res_cnv2)[3] <- "B1_2"
res_nocnv2$diffexpressed <- "NO"
res_nocnv2$diffexpressed[res_nocnv2$B1_1 >= 0.6 & res_nocnv2$padj < 0.05] <- "UP"
res_nocnv2$diffexpressed[res_nocnv2$B1_1 <= -0.6 & res_nocnv2$padj < 0.05] <- "DOWN"

#Make simple graphics
p1 <- ggplot(data = res_nocnv1, aes(x = B1_1, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "#bb0c00"))+
  scale_x_continuous(breaks = seq(-14, 10, 4))+
  labs(title="DGE Tum vs Norm (M1)",x="Effect size (B1_1)")+
  theme(legend.position="none")+
  #theme_minimal()+
  font("xy.text", size = 14, color = "black", face = "bold")+
  font("xlab", size = 14)+
  font("ylab", size = 14)
p1

p2 <- ggplot(data = res_cnv1, aes(x = B1_2, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "#bb0c00"))+
  scale_x_continuous(breaks = seq(-14, 10, 4))+
  labs(title="DGE Tum vs Norm (M2)",x="Effect size (B1_2)")+
  theme(legend.position="none")+
  #theme_minimal()+
  font("xy.text", size = 14, color = "black", face = "bold")+
  font("xlab", size = 14)+
  font("ylab", size = 14)
p2

p3 <- ggplot(data = res_nocnv2, aes(x = B1_1, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "#bb0c00"))+
  scale_x_continuous(breaks = seq(-14, 10, 4))+
  labs(title="DGE Tum vs Norm (M1)",x="Effect size (B1_1)")+
  theme(legend.position="none")+
  #theme_minimal()+
  font("xy.text", size = 14, color = "black", face = "bold")+
  font("xlab", size = 14)+
  font("ylab", size = 14)
p3

p4 <- ggplot(data = res_cnv2, aes(x = B1_2, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("black", "#00AFBB", "#bb0c00"))+
  scale_x_continuous(breaks = seq(-10, 10, 4))+
  labs(title="DGE Tum vs Norm (M2)",x="Effect size (B1_2)")+
  theme(legend.position="none")+
  #theme_minimal()+
  font("xy.text", size = 14, color = "black", face = "bold")+
  font("xlab", size = 14)+
  font("ylab", size = 14)
p4

p5 <- ggplot(data = stat_res_luad, aes(x = B1_1, y = -log10(padj_1), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "#bb0c00"))+
  scale_x_continuous(breaks = seq(-10, 10, 4))+
  labs(title="DGE Tum vs Norm (M1)",x="Effect size (B1_1)")+
  theme(legend.position="none")+
  #theme_minimal()+
  font("xy.text", size = 14, color = "black", face = "bold")+
  font("xlab", size = 14)+
  font("ylab", size = 14)
p5

#delete rows by name
stat_res_luad <- stat_res_luad[!(row.names(stat_res_luad) %in% c("NXPH3")),]

p6 <- ggplot(data = stat_res_luad, aes(x = B1_2, y = -log10(padj_2), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "#bb0c00"))+
  scale_x_continuous(breaks = seq(-10, 10, 4))+
  labs(title="DGE Tum vs Norm (M2)",x="Effect size (B1_2)")+
  theme(legend.position="none")+
  #theme_minimal()+
  font("xy.text", size = 14, color = "black", face = "bold")+
  font("xlab", size = 14)+
  font("ylab", size = 14)
p6

#Plots
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)
#grid.arrange(g2, arrangeGrob(g3, g4, ncol=2), nrow = 2)

#Density plot

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




