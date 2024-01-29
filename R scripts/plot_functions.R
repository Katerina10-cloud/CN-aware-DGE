#Plots
library(ggplot2)
library(ggpubr)

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


plot_data_1 <- deg_b1 %>% mutate(effect_size = "B1_1")
plot_data_2 <- deg_b2 %>% mutate(effect_size = "B1_2")
plot_data <- rbind(plot_data_1, plot_data_2)

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
violin_plot <- ggplot(deg, aes(x = cnv, y = difference, fill = cnv))+
  geom_violin(trim=FALSE)+
  #geom_jitter(shape=10, position=position_jitter(0.1))+
  labs(title="CNV patterns and Effect size difference ((|B1_2| - |B1_1|), (n = 7522 genes))",x="CNV group", y = "Effect size difference (log2FC)")+
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


# Volcano Plot
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
colnames(res_sint_cnv_heterog)[3] <- "B1_2"
res_sint_cnv_heterog$diffexpressed <- "NO"
res_sint_cnv_heterog$diffexpressed[res_sint_cnv_heterog$B1_1 >= 0.6 & res_sint_cnv_heterog$padj < 0.05] <- "UP"
res_sint_cnv_heterog$diffexpressed[res_sint_cnv_heterog$B1_1 <= -0.6 & res_sint_cnv_heterog$padj < 0.05] <- "DOWN"

ggplot(data = res_sint_nocnv, aes(x = B1_1, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "#bb0c00"), labels = c("Down-regulated", "Not significant", "Up-regulated"))+
  scale_x_continuous(breaks = seq(-14, 10, 2))+
  labs(title="Tumor vs Normal (BRCA, sample size = 100)",x="Effect size (B1_1)")+
  theme_classic()

ggplot(data = res_sint_cnv, aes(x = B1_2, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "#bb0c00"), labels = c("Down-regulated", "Not significant", "Up-regulated"))+
  scale_x_continuous(breaks = seq(-14, 10, 2))+
  labs(title="Tumor vs Normal (BRCA, sample size = 100)",x="Effect size (B1_2)")+
  theme_classic()

ggplot(data = res_sint_nocnv_heterog, aes(x = B1_1, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "#bb0c00"), labels = c("Down-regulated", "Not significant", "Up-regulated"))+
  scale_x_continuous(breaks = seq(-14, 10, 2))+
  labs(title="Tumor vs Normal (BRCA, sample size = 100)",x="Effect size (B1_1)")+
  theme_classic()

ggplot(data = res_sint_control, aes(x = B1_1, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "#bb0c00"), labels = c("Down-regulated", "Not significant", "Up-regulated"))+
  scale_x_continuous(breaks = seq(-14, 10, 2))+
  labs(title="Normal_1 vs Normal_2 (BRCA, sample size = 100)",x="Effect size (B1_1)")+
  theme_classic()


ggarrange(
  plot_1, plot_2, plot_3, plot_4, ncol = 4)

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




