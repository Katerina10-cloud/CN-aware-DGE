setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "colorspace", "ggpubr", "gridExtra", "ggpointdensity", "metaseqR2", "ggalluvial",
          "ggridges", "ggforce", "ggparallel", "alluvial")
sapply(pkgs, require, character.only = TRUE)

### Scatter plot (RNA LFC and CN relationship) ###

res_naive <- read.csv("CN-aware-DGE/Python/results/COAD/res_CNnaive_test.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/COAD/res_CNaware_test.csv")
cnv_filt <- read.csv("TCGA/colon/test/cnv_test.csv")

res_naive <- res_naive %>% dplyr::select(X,log2FoldChange, padj) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_adj <- res_adj %>% dplyr::select(X,log2FoldChange, padj) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

cnv_filt <- cnv_filt %>% remove_rownames %>% column_to_rownames(var="X")
cnv_filt <- cnv_filt[21:40] * 2 # LUAD
cnv_filt <- cnv_filt[58:114] * 2 # BRCA
cnv_filt <- cnv_filt[17:34] * 2 # LIHC
cnv_filt <- cnv_filt[22:42] * 2 # HNSC
cnv_filt <- cnv_filt[13:24] * 2 # COAD

cnv_mean <- cnv_filt %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  dplyr::select(cnv_mean) 
  #dplyr::mutate(tumor_type = "COAD")

cnv_mean <- rbind(cnv_mean_luad, cnv_mean_brca, cnv_mean_lihc, cnv_mean_hnsc, cnv_mean_coad)

hist <- ggplot(cnv_mean, aes(x = cnv_mean)) +
  geom_histogram(binwidth = 0.4, fill = "#DF8F4499", color = "black") +
  labs(
    title = "",
    x = "CN state",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 17, color = "black", hjust = 1),
    axis.text.y = element_text(size = 17, color = "black"),
    axis.title.x = element_text(size = 16, face = "plain"),
    axis.title.y = element_text(size = 16, face = "plain"),
    strip.text = element_text(size = 16, face = "plain")
  ) +
  facet_wrap(~factor(tumor_type, levels = c("LUAD", "BRCA", "LIHC", "HNSC", "COAD")), nrow = 1)
hist

ggsave("CN-aware-DGE/plots/supplementary/hist_tumor_types.png", dpi = 400, width = 12.0, height = 3.5, plot = hist)

common_genes <- intersect(rownames(cnv_mean), rownames(res_naive))
cnv_mean <- cnv_mean[common_genes, ] %>% data.frame()

colnames(cnv_mean) <- c("cnv_mean")
cnv_mean <- cnv_mean %>% 
  dplyr::mutate(cnv_group = case_when(
    cnv_mean > 1.5 & cnv_mean < 2.5 ~ "diploid",
    cnv_mean >= 2.5 & cnv_mean < 4.5 ~ "gain",
    cnv_mean >= 4.5  ~ "amplification")) 

res_naive <- res_naive %>% dplyr::select(log2FoldChange, padj)
res_adj <- res_adj %>% dplyr::select(log2FoldChange, padj)

res_naive_edge <- res_naive_edge %>% dplyr::select(logFC, padj)
res_aware_edge <- res_aware_edge %>% dplyr::select(logFC, padj)

p_naive <- merge(cnv_mean, res_naive, by = "row.names")
p_adj <- merge(cnv_mean, res_adj, by = "row.names")

colnames(p_naive) <- c("geneID", "cnv_mean", "cnv_group", "logFC", "padj")
colnames(p_adj) <- colnames(p_naive)


d_scatter_pydeseq <- rbind(p_naive %>% dplyr::mutate(method = "CN naive"), 
                   p_adj %>% dplyr::mutate(method = "CN aware"))

d_scatter_pydeseq <- d_scatter_pydeseq %>% dplyr::mutate(tool = "PyDESeq2")
d_scatter_pydeseq <- d_scatter_pydeseq %>% dplyr::filter(abs(logFC) < 5.0 ,)

d_scatter_edge <- rbind(p_naive %>% dplyr::mutate(method = "CN naive"), 
                           p_adj %>% dplyr::mutate(method = "CN aware"))
d_scatter_edge <- d_scatter_edge %>% dplyr::mutate(tool = "edgeR")
d_scatter_edge <- d_scatter_edge %>% dplyr::filter(logFC > -6.0 ,)

d_scatter <- rbind(d_scatter_edge, d_scatter_pydeseq)

# Tumor types datasets
d_scatter_luad <- d_scatter_pydeseq %>% dplyr::mutate(tumor_type = "LUAD")
d_scatter_brca <- d_scatter_pydeseq %>% dplyr::mutate(tumor_type = "BRCA")
d_scatter_lihc <- d_scatter_pydeseq %>% dplyr::mutate(tumor_type = "LIHC")
d_scatter_hnsc <- d_scatter_pydeseq %>% dplyr::mutate(tumor_type = "HNSC")
d_scatter_coad <- d_scatter_pydeseq %>% dplyr::mutate(tumor_type = "COAD")

d_scatter <- rbind(d_scatter_luad, d_scatter_brca, d_scatter_lihc, d_scatter_hnsc, d_scatter_coad)

cn_group_colors <- c("diploid" = "#ADB6B6", "gain" = "#d1ab75", "amplification" = "#B48EAD")

scatter = ggplot(d_scatter, aes(x = log2(cnv_mean), y = logFC)) +
  geom_point(size=0.6, aes(colour = factor(cnv_group))) +
  #geom_smooth(method = "lm", color = "blue", se = T)+
  geom_smooth()+
  labs(x ="log2(CN ratio)", y="log2FC", colour = "CN group") +
  theme(plot.title=element_text(hjust=0.7, vjust=0.7))+
  #ggh4x::facet_nested(factor(tool, levels = c("PyDESeq2", "edgeR"))~factor(method, levels = c("CN naive", "CN aware")), scales ="free", independent = "y")+
  #facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow=1, scale = "free_y")+
  ggh4x::facet_nested(factor(method, levels = c("CN naive", "CN aware"))~factor(tumor_type, levels = c("LUAD", "BRCA", "LIHC", "HNSC", "COAD")), scales ="free", independent = "y")+
  guides(color = guide_legend(override.aes = list(size=2)))+
  theme_bw()+
  scale_color_manual(values = c(cn_group_colors))+
  ggplot2::theme(legend.position = 'bottom',
                 legend.text = element_text(size = 12),
                 legend.title = element_text(size = 14),
                 strip.text = element_text(size = 14, face = "plain"))
scatter


ggsave("CN-aware-DGE/plots/supplementary/scatter_cancer_types.png", dpi = 400, width = 10.0, height = 4.5, plot = scatter)


# Log2FC comparison  - scatter plot #

lfc_naive <- res_naive_edge %>% dplyr::select(logFC) %>% dplyr::rename(logFC=logFC)
lfc_adj <- res_adj_edge %>% dplyr::select(logFC) %>% dplyr::rename(logFC_adj=logFC)

res_naive_edge$geneID <- rownames(res_naive_edge)
res_adj_edge$geneID <- rownames(res_adj_edge)

common_genes <- intersect(rownames(res_naive_edge), rownames(res_adj_edge))
res_naive_edge <- res_naive_edge[common_genes, ]
res_adj_edge <- res_adj_edge[common_genes, ]

lfc_naive <- res_naive_edge %>% dplyr::select(logFC) %>% dplyr::rename(logFC=logFC)
lfc_adj <- res_adj_edge %>% dplyr::select(logFC) %>% dplyr::rename(logFC_adj=logFC)

plot_lfc_edge <- merge(lfc_naive, lfc_adj, by = "row.names")
plot_lfc_pydeseq <- merge(lfc_naive, lfc_adj, by = "row.names")

#plot_lfc_edge <- subset(plot_lfc_edge, plot_lfc_edge$logFC < 2.0 & plot_lfc_edge$logFC > -2.0)
colnames(plot_lfc_edge) <- colnames(plot_lfc_pydeseq)

d_lfc <- rbind(plot_lfc_edge %>% dplyr::mutate(method = "edgeR"), 
                   plot_lfc_pydeseq %>% dplyr::mutate(method = "PyDESeq2"))


comparison_lfc <- ggplot(d_lfc, aes(x=logFC_adj, y=logFC)) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=TRUE, linetype="dashed",color="darkred")+
  xlab("log2FC CN aware") +
  ylab ("log2FC CN naive") +
  facet_wrap(~method)+
  theme_classic()+
  theme(legend.position="none")
comparison_lfc


# p-value comparison #

pval_naive <- res_naive_edge %>% dplyr::select(padj) %>% dplyr::rename(padj=padj)
pval_adj <- res_adj_edge %>% dplyr::select(padj) %>% dplyr::rename(padj = padj)
plot_pval_edge <- merge(pval_naive, pval_adj, by = "row.names")

pval_naive <- res_naive_pydeseq %>% dplyr::select(padj) %>% dplyr::rename(padj=padj)
pval_adj <- res_adj_pydeseq %>% dplyr::select(padj) %>% dplyr::rename(padj_adj = padj)
plot_pval_pydeseq <- merge(pval_naive, pval_adj, by = "row.names")

colnames(plot_pval_pydeseq) <- c("geneID", "padj", "padj_adj")
colnames(plot_pval_edge) <- colnames(plot_pval_pydeseq)

d_pval <- rbind(plot_pval_edge %>% dplyr::mutate(method = "edgeR"), 
               plot_pval_pydeseq %>% dplyr::mutate(method = "PyDESeq2"))

comparison_pval <- ggplot(d_pval, aes(x=-log10(padj), y=-log10(padj_adj))) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",color="darkred")+
  xlab("-log10(padj) CN aware") +
  ylab ("-log10(padj) CN naive") +
  facet_wrap(~method)+
  theme_classic()+
  theme(legend.position="none")
comparison_pval

gridExtra::grid.arrange(comparison_lfc, comparison_pval, nrow = 2)


### Volcano plot ###

res_naive_pydeseq <- read.csv("CN-aware-DGE/Python/results/LIHC/res_CNnaive_test2.csv")
res_adj_pydeseq <- read.csv("CN-aware-DGE/Python/results/LIHC/res_CNaware_test2.csv")

lfc_cut <- 1.0
pval_cut <- .05
de_gene_colors <- c("Not significant" = "gray", "Down-reg" = "#8491B4B2", "Up-reg"="#E7969C")

res_adj_pydeseq <- res_adj_pydeseq %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  #dplyr::mutate(tool = "PyDESeq2") %>% 
  dplyr::mutate(tumor_type = "LIHC") %>% 
  dplyr::mutate(method = "CN aware") %>%
  #dplyr::mutate(n_genes = "3000 genes") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, tumor_type, method) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_naive_pydeseq <- res_naive_pydeseq %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN naive") %>% 
  dplyr::mutate(tumor_type = "LIHC") %>% 
  #dplyr::mutate(tool = "PyDESeq2") %>%
  #dplyr::mutate(n_genes = "3000 genes") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, tumor_type, method) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")


res_naive_pydeseq <- res_naive_pydeseq[!(row.names(res_naive_pydeseq) %in% c("PYCR1")),] #LUAD
res_adj_pydeseq <- res_adj_pydeseq[!(row.names(res_adj_pydeseq) %in% c("CAV1")),] #LUAD
#res_naive_edge <- res_naive_edge[!(row.names(res_naive_edge) %in% c("PYCR1")),]

colnames(res_naive_pydeseq) <- c("logFC", "padj", "isDE", "DEtype", "tumor_type", "method")
colnames(res_adj_pydeseq) <- colnames(res_naive_pydeseq) 

d_volcano_luad <- rbind(res_naive_pydeseq, res_adj_pydeseq)
d_volcano_brca <- rbind(res_naive_pydeseq, res_adj_pydeseq)
d_volcano_lihc <- rbind(res_naive_pydeseq, res_adj_pydeseq)
d_volcano_hnsc <- rbind(res_naive_pydeseq, res_adj_pydeseq)
d_volcano_coad <- rbind(res_naive_pydeseq, res_adj_pydeseq)

d_volcano_brca <- d_volcano_brca %>% dplyr::filter(logFC < 5.0 ,)
d_volcano_lihc <- d_volcano_lihc %>% dplyr::filter(abs(logFC) < 5.0 ,) 
d_volcano_hnsc <- d_volcano_hnsc %>% dplyr::filter(abs(logFC) < 6.0 ,) 
d_volcano_coad <- d_volcano_coad %>% dplyr::filter(abs(logFC) < 5.0 ,) 

d_volcano <- rbind(d_volcano_luad, d_volcano_brca, d_volcano_lihc, d_volcano_hnsc, d_volcano_coad)

#colnames(d_volcano_luad) <- c("logFC", "padj", "isDE", "DEtype", "method", "tumor_type")

#d_volcano_3000 <- d_volcano_3000[!(row.names(d_volcano_3000) %in% c("g259")),]
#d_volcano_3000 <- d_volcano_3000 %>% dplyr::filter(log2FoldChange < 1.9 ,)
#d_volcano_1000 <- rbind(res_naive_pydeseq, res_adj_pydeseq)
#d_volcano_3000 <- rbind(res_naive_pydeseq, res_adj_pydeseq)
#d_volcano_sim <- rbind(d_volcano_1000, d_volcano_3000)
#colnames(d_volcano_sim) <- c("logFC", "padj", "isDE", "DEtype", "method", "n_genes")

res_adj_edge <- res_adj_edge %>%
  dplyr::mutate(isDE = (abs(logFC) >= lfc_cut) & (FDR <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(logFC > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(tool = "edgeR") %>% 
  dplyr::mutate(method = "CN aware") %>% 
  dplyr::select(logFC, FDR, isDE, DEtype, method, tool) 

res_naive_edge <- res_naive_edge %>%
  dplyr::mutate(isDE = (abs(logFC) >= lfc_cut) & (FDR <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(logFC > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN naive") %>% 
  dplyr::mutate(tool = "edgeR") %>% 
  dplyr::select(logFC, FDR, isDE, DEtype, method, tool) 

colnames(res_naive_edge) <- c("logFC", "padj", "isDE", "DEtype", "method", "tool")

d_volcano_edge <- rbind(res_naive_edge, res_aware_edge)
d_volcano_pydeseq <- rbind(res_naive_pydeseq, res_adj_pydeseq)
colnames(d_volcano_edge) <- c("logFC", "padj", "isDE", "DEtype", "method", "tool")
colnames(d_volcano_pydeseq) <- colnames(d_volcano_edge)

common_genes <- intersect(rownames(d_volcano_edge), rownames(d_volcano_pydeseq))
d_volcano_pydeseq <- d_volcano_pydeseq[common_genes, ]
d_volcano_edge <- d_volcano_edge[common_genes, ]
d_volcano <- rbind(d_volcano_pydeseq, d_volcano_edge)

#d_volcano <- rbind(res_naive_pydeseq, res_adj_pydeseq, res_naive_edge, res_adj_edge)

p_volcanos <-  d_volcano %>% 
  ggplot(mapping = aes(x=logFC, y=-log10(padj), col=DEtype)) +
  geom_point(size=.9) +
  theme_bw() +
  scale_color_manual(values = de_gene_colors) +
  #ggh4x::facet_nested(factor(tool, levels = c("PyDESeq2", "edgeR"))~factor(method, levels = c("CN naive", "CN aware")), scales ="free", independent = "y")+
  #ggh4x::facet_nested(factor(n_genes, levels = c("1000 genes", "3000 genes"))~factor(method, levels = c("CN naive", "CN aware")), scales ="free", independent = "y")+
  ggh4x::facet_nested(factor(method, levels = c("CN naive", "CN aware"))~factor(tumor_type, levels = c("LUAD", "BRCA", "LIHC", "HNSC", "COAD")), scales ="free", independent = "y")+
  #facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow=1, scale = "free")+
  ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col="") +
  ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  ggplot2::geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'bottom',
                 legend.text = element_text(size = 14),
                 strip.text = element_text(size = 14, face = "plain"))
p_volcanos

ggsave("CN-aware-DGE/plots/supplementary/volcano_cancer_types.png", dpi = 400, width = 10.0, height = 4.5, plot = p_volcanos)


### Impact of CN normalization on DGE analysis ###

res_naive <- read.csv("CN-aware-DGE/Python/results/COAD/res_CNnaive_test.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/COAD/res_CNaware_test.csv")

lfc_cut <- 1.0
pval_cut <- .05

res_adj <- res_adj %>%
  #dplyr::filter(padj < pval_cut, abs(log2FoldChange) > lfc_cut) %>% 
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DE_type = if_else(!isDE, "n.s", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN aware") %>% 
  dplyr::mutate(cancer_type = "COAD") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DE_type, method, cancer_type) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_naive <- res_naive %>%
  #dplyr::filter(padj < pval_cut, abs(log2FoldChange) > lfc_cut) %>% 
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DE_type = if_else(!isDE, "n.s", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN naive") %>% 
  dplyr::mutate(cancer_type = "COAD") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DE_type, method, cancer_type) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

common_genes <- intersect(rownames(res_naive), rownames(res_adj))
res_naive <- res_naive[common_genes, ] 
res_adj <- res_adj[common_genes, ] 

p_data_luad <- rbind(res_naive, res_adj)
p_data_brca <- rbind(res_naive, res_adj)
p_data_lihc <- rbind(res_naive, res_adj)
p_data_hnsc <- rbind(res_naive, res_adj)
p_data_coad <- rbind(res_naive, res_adj)

p_data_coad <- p_data_coad %>%
  group_by(method, DE_type) %>%
  summarise(count = n()) %>%  # Count occurrences of each DEtype per method
  mutate(percentage = count / sum(count) * 100) %>% 
  dplyr::mutate(cancer_type = "COAD")

p_data_coad$method <- factor(p_data_coad$method, levels = c("CN naive", "CN aware"))
p_data <- rbind(p_data_luad, p_data_brca, p_data_lihc, p_data_hnsc, p_data_coad)

#de_gene_colors <- c("Not significant" = "darkgray", "Down-reg" = "lightskyblue2", "Up-reg"="plum")

#p_bar <- ggplot2::ggplot(p_data, aes(x = method, y = percentage, fill = DE_type)) +
  #ggplot2::geom_bar(stat = "identity", width=0.5) +
  #ggplot2::geom_text(aes(label = count), 
            #position = position_stack(vjust = 0.5), size = 3)+
  #scale_fill_manual(values = de_gene_colors) +
  #labs(x = "", y = "") +
  #facet_wrap(~factor(cancer_type, levels = c("LUAD", "BRCA", "LIHC", "HNSC", "COAD")), nrow=1)+
  #theme_bw()+
  #ggplot2::theme(legend.position = 'bottom')
#p_bar


# Sankey diagram #

CN_naive = c("Down-reg", "n.d.", "Up-reg")
CN_aware = c("Down-reg", "n.d.", "Up-reg")

luad_naive <- c(638, 3179, 598)
luad_aware <- c(1298, 2897, 220)
freq_luad <- c(638,1298,3179,2897,698,220)

brca_naive <- c(558,4042,514)
brca_aware <- c(1215,3752,147)
freq_brca <- c(558,1215,4042,3752,514,147)

lihc_naive <- c(456,3135,354)
lihc_aware <- c(1241, 2640,64)
freq_lihc <- c(456,1241,3135,2640,354,64)

hnsc_naive <- c(449,3828,358)
hnsc_aware <- c(869,3582,184)
freq_hnsc <- c(449,869,3828,3582,358,184)

coad_naive <- c(449,3454,588)
coad_aware <- c(689,3561,279)
freq_coad <- c(449,689,3454,3561,588,279)
  
  
data <- data.frame(
  CN_naive = CN_naive,
  CN_aware = CN_aware,
  freq_CN_naive = coad_naive,  
  freq_CN_aware = coad_aware  
)

data_coad <- data.frame(
  CN_naive = rep(data$CN_naive, each = 2),
  CN_aware = rep(data$CN_aware, each = 2),
  freq = freq_coad,
  method = rep(c("CN_naive", "CN_aware"), times = 3)
) %>% dplyr::mutate(DE_group = CN_aware) %>% dplyr::mutate(cancer_type = "COAD")

data_long <- rbind(data_luad, data_brca, data_lihc, data_hnsc, data_coad)

`DE group` <- c("Down-reg" = "#8491B4B2", "n.s." = "lightgray", "Up-reg" = "#F5A2A2")

ggplot(data_luad,
       aes(x = method, stratum = CN_naive, alluvium = CN_aware, y = freq)) +
  geom_alluvium(aes(fill = DE_group)) +
  geom_stratum(aes(fill = DE_group)) +
  geom_text(stat = "stratum", aes(label = freq), size = 3, vjust = 0.5) +
  scale_x_discrete(limits = c("CN_naive", "CN_aware"), expand = c(0.2, 0.2)) +
  #facet_wrap(~cancer_type, nrow=1)+
  facet_wrap(~factor(cancer_type, levels = c("LUAD", "BRCA", "LIHC", "HNSC", "COAD")), nrow=1)+
  scale_fill_manual(values = `DE group`)+
  labs(y = "", x = "", fill = "DE group") +
  theme_classic() +
  theme(legend.position = "bottom") 


# Sankey dynamic gene groups transitions 

res_naive <- res_naive %>% dplyr::select(DE_type) %>% dplyr::rename(CN_naive = DE_type)
res_adj <- res_adj %>% dplyr::select(DE_type) %>% dplyr::rename(CN_aware = DE_type)

res_join <- cbind(res_naive, res_adj)

data_flow <- res_join %>%
  group_by(CN_naive,CN_aware) %>%
  summarise(freq = n()) %>%
  ungroup()

data_ggforce <- data_flow  %>%
  gather_set_data(1:2) %>%        
  arrange(x,CN_naive,desc(CN_aware))

data_ggforce$CN_naive <- factor(data_ggforce$CN_naive)
data_ggforce$CN_aware <- factor(data_ggforce$CN_aware)

data_ggforce <- data_ggforce %>%
  group_by(CN_naive, CN_aware) %>%
  mutate(y_mid = freq / 2)

g_group_colors <- c("Down-reg" = "#8491B4B2", "n.s" = "lightgray", "Up-reg" = "#F5A2A2")

sankey_coad <- ggplot(data_ggforce, aes(x = x, id = id, split = y, value = freq)) +
  geom_parallel_sets(aes(fill = CN_naive), alpha = 0.9, axis.width = 0.2,
                     n = 4415, strength = 0.5, color = "black", linewidth = 0.3) +
  geom_parallel_sets_axes(axis.width = 0.25, fill = "gray93",
                          color = "gray80", linewidth = 0.5) +  
  #geom_parallel_sets_labels(colour = 'gray35', size = 2.0, angle = 0, fontface = "plain") +
  scale_fill_manual(values = g_group_colors, name = "Gene group") +
  scale_color_manual(values = g_group_colors) +
  scale_x_continuous(breaks = 1:2, labels = c("CN-naive", "CN-aware"))+
  theme_minimal() +
  theme(
    legend.position = "",
    legend.title = element_text(size = 15, face = "plain"),
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x  = element_blank()
  )
sankey_coad



ggsave("CN-aware-DGE/plots/supplementary/sankey_brca.png", dpi = 400, width = 6.0, height = 7.0, plot = sankey_brca)

# Define gene groups selection #

res_naive <- read.csv("CN-aware-DGE/Python/results/LUAD/res_CNnaive_test1.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/LUAD/res_CNaware_test1.csv")

lfc_cut <- 1.0
pval_cut <- .05

res_naive <- res_naive %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DE_type = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::select(X,log2FoldChange, padj, isDE, DE_type)

colnames(res_naive) <- c("geneID", "log2FC_naive", "padj_naive", "isDE_naive", "DE_type_naive")

res_adj <- res_adj %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DE_type = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::select(X,log2FoldChange, padj, isDE, DE_type)

colnames(res_adj) <- c("geneID", "log2FC_adj", "padj_adj", "isDE_adj", "DE_type_adj")

res_joint <- dplyr::left_join(res_naive, res_adj, by = "geneID")

dosage_insensitive_up <- res_joint %>%
  dplyr::filter(padj_adj < pval_cut & log2FC_adj > lfc_cut) %>% 
  dplyr::mutate(lfc_diff = abs(log2FC_adj - abs(log2FC_naive)))
  
dosage_insensitive_down <- res_joint %>% 
  dplyr::filter(padj_adj < pval_cut & log2FC_adj < -lfc_cut) %>% 
  dplyr::mutate(lfc_diff = log2FC_adj + abs(log2FC_naive))

not_significant_adj <- res_joint %>% 
  dplyr::filter(padj_adj > pval_cut & log2FC_adj < lfc_cut | padj_adj > pval_cut & log2FC_adj > -lfc_cut) %>% 
  dplyr::mutate(lfc_diff = log2FC_adj + abs(log2FC_naive)) %>% 
  dplyr::select(geneID, log2FC_adj, padj_adj) %>% 
  dplyr::mutate(method = "CN aware")
colnames(not_significant_adj) <- c("geneID", "log2FC", "padj", "method")

not_significant_naive <- res_joint %>% 
  dplyr::filter(padj_naive > pval_cut & log2FC_naive < lfc_cut | padj_naive > pval_cut & log2FC_naive > -lfc_cut) %>% 
  dplyr::select(geneID, log2FC_naive, padj_naive) %>% 
  dplyr::mutate(method = "CN naive")
colnames(not_significant_naive) <- c("geneID", "log2FC", "padj", "method")

not_significant <- rbind(not_significant_naive, not_significant_adj)

deg_adj <- rbind(dosage_insensitive_up, dosage_insensitive_down)
deg_adj <- deg_adj %>% 
  dplyr::select(geneID, log2FC_adj, padj_adj) %>% 
  dplyr::mutate(method = "CN aware")
colnames(deg_adj) <- c("geneID", "log2FC", "padj", "method")
  
deg_naive <- res_joint %>% 
  dplyr::filter(padj_naive < pval_cut & log2FC_naive > lfc_cut | padj_naive < pval_cut & log2FC_naive < -lfc_cut) %>% 
  dplyr::select(geneID, log2FC_naive, padj_naive) %>% 
  dplyr::mutate(method = "CN naive")
colnames(deg_naive) <- c("geneID", "log2FC", "padj", "method")

deg <- rbind(deg_naive, deg_adj)

p1 <- ggplot(not_significant, aes(x = padj, y = method, fill = method)) +
  geom_density_ridges(alpha = 0.8, scale = 0.8) +  # Adjust alpha for transparency and scale for height
  scale_x_log10() +  
  theme_classic() +  
  theme(axis.title.y = element_blank(),  
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = "not DE genes", x = "p-value") +
  scale_fill_manual(values = c("#EEBD91", "#C9d175"))+
  theme(legend.position = "bottom")
p1


p2 <- ggplot(deg, aes(x = -log10(padj), y = method, fill = method)) +
  geom_density_ridges(alpha = 0.8, scale = 1.0) +  
  #scale_x_log10() +  
  theme_classic() +  
  theme(axis.title.y = element_blank(),  
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = "DE genes", x = "-log10(p-value)") +
  scale_fill_manual(values = c("#EEBD91", "#C9d175"))+
  theme(legend.position = "bottom")
p2

p1+p2
