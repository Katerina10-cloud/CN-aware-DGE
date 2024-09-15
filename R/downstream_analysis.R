setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "colorspace", "ggpubr", "gridExtra", "ggpointdensity", "metaseqR2", "ggalluvial")
sapply(pkgs, require, character.only = TRUE)

### Scatter plot (RNA LFC and CN relationship) ###

res_naive <- read.csv("CN-aware-DGE/Python/results/res_CNnaive_test1.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/res_CNaware_test1.csv")
cnv_filt <- read.csv("TCGA/lung/LUAD/cnv_test_1.csv")

res_naive <- res_naive %>% dplyr::select(X,log2FoldChange, padj) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_adj <- res_adj %>% dplyr::select(X,log2FoldChange, padj) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

cnv_filt <- cnv_filt %>% 
  as.data.frame() %>% 
  remove_rownames %>% 
  column_to_rownames(var="X") 

cnv_mean <- cnv_filt %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  dplyr::select(cnv_mean)

common_genes <- intersect(rownames(cnv_mean), rownames(res_naive))
cnv_mean <- cnv_mean[common_genes, ] %>% data.frame()

res_naive <- res_naive_edge %>% dplyr::select(logFC, FDR)
res_adj <- res_adj_edge %>% dplyr::select(logFC, FDR)

p_naive <- merge(cnv_mean, res_naive, by = "row.names")
p_adj <- merge(cnv_mean, res_adj, by = "row.names")

colnames(p_naive) <- c("geneID", "cnv_mean", "logFC", "padj")
colnames(p_adj) <- colnames(p_naive)

#p_naive <- plot_data %>% dplyr::filter(logFC > -2.5,)
#p_adj <- p_adj %>% dplyr::filter(logFC > -4.0,)

d_scatter_pydeseq <- rbind(p_naive %>% dplyr::mutate(method = "CN naive"), 
                   p_adj %>% dplyr::mutate(method = "CN aware"))
d_scatter_pydeseq <- d_scatter_pydeseq %>% 
  dplyr::mutate(tool = "PyDESeq2")
d_scatter_pydeseq <- d_scatter_pydeseq %>% dplyr::filter(padj > 0.000000e+00,)



d_scatter_edge <- rbind(p_naive %>% dplyr::mutate(method = "CN naive"), 
                           p_adj %>% dplyr::mutate(method = "CN aware"))
d_scatter_edge <- d_scatter_edge %>% 
  dplyr::mutate(tool = "edgeR")

d_scatter <- rbind(d_scatter_edge, d_scatter_pydeseq)

d_scatter <- d_scatter %>% dplyr::filter(logFC < 6.0,)


scatter = ggplot(d_scatter, aes(x = log2(cnv_mean), y = logFC)) +
  geom_point(size=0.8) +
  #geom_smooth(method = "lm", color = "blue", se = T)+
  geom_smooth()+
  #theme(legend.position = 'bottom') +
  labs(x ="log2(CN ratio)", y="log2FC") +
  theme(plot.title=element_text(hjust=0.7, vjust=0.7))+
  #facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow=1, scale = "free_y")+
  ggh4x::facet_nested(factor(tool, levels = c("PyDESeq2", "edgeR"))~factor(method, levels = c("CN naive", "CN aware")), scales ="free", independent = "y")+
  guides(color = guide_legend(override.aes = list(size=2)))+
  theme_bw()
scatter

#gridExtra::grid.arrange(scatter1, scatter2, nrow = 2)


# Log2FC comparison  - scatter plot #

lfc_naive <- res_adj_edge %>% dplyr::select(logFC) %>% dplyr::rename(logFC=logFC)
lfc_adj <- res_naive_edge %>% dplyr::select(logFC) %>% dplyr::rename(logFC_adj=logFC)

res_naive_pydeseq$geneID <- rownames(res_naive_pydeseq)
res_adj_pydeseq$geneID <- rownames(res_adj_pydeseq)

common_genes <- intersect(rownames(res_naive_pydeseq), rownames(res_adj_pydeseq))
res_naive_pydeseq <- res_naive_pydeseq[common_genes, ]
res_adj_pydeseq <- res_adj_pydeseq[common_genes, ]

lfc_naive <- res_adj_pydeseq %>% dplyr::select(log2FoldChange) %>% dplyr::rename(logFC=log2FoldChange)
lfc_adj <- res_naive_pydeseq %>% dplyr::select(log2FoldChange) %>% dplyr::rename(logFC_adj=log2FoldChange)

plot_lfc_edge <- merge(lfc_naive, lfc_adj, by = "row.names")
plot_lfc_pydeseq <- merge(lfc_naive, lfc_adj, by = "row.names")

#plot_lfc_edge <- subset(plot_lfc_edge, plot_lfc_edge$logFC < 2.0 & plot_lfc_edge$logFC > -2.0)
colnames(plot_lfc_edge) <- colnames(plot_lfc_pydeseq)

d_lfc <- rbind(plot_lfc_edge %>% dplyr::mutate(method = "edgeR"), 
                   plot_lfc_pydeseq %>% dplyr::mutate(method = "PyDESeq2"))


comparison_lfc <- ggplot(d_lfc, aes(x=logFC_adj, y=logFC)) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=TRUE, linetype="dashed",color="darkred")+
  xlab("log2FC with CN normalization") +
  ylab ("log2FC no CN normalization") +
  facet_wrap(~method, scales = "free")+
  theme_bw()+
  theme(legend.position="none")
comparison_lfc


# p-value comparison #

pval_naive <- res_naive_edge %>% dplyr::select(FDR) %>% dplyr::rename(padj=FDR)
pval_adj <- res_adj_edge %>% dplyr::select(FDR) %>% dplyr::rename(padj = FDR)
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
  xlab("padj with CN normalization") +
  ylab ("padj no CN normalization") +
  facet_wrap(~method)+
  theme_bw()+
  theme(legend.position="none")
comparison_pval

gridExtra::grid.arrange(comparison_lfc, comparison_pval, nrow = 1)


### Volcano plot ###

res_naive_pydeseq <- read.csv("CN-aware-DGE/Python/results/res_CNnaive_test2.csv")
res_adj_pydeseq <- read.csv("CN-aware-DGE/Python/results/res_CNaware_test2.csv")

lfc_cut <- 1.0
pval_cut <- .05
de_gene_colors <- c("Not significant" = "gray", "Down-reg" = "darkblue", "Up-reg"="darkred")

res_adj_pydeseq <- res_adj_pydeseq %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(tool = "PyDESeq2") %>% 
  dplyr::mutate(method = "CN aware") %>%
  #dplyr::mutate(n_genes = "1000 genes") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, method, tool) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_naive_pydeseq <- res_naive_pydeseq %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN naive") %>% 
  dplyr::mutate(tool = "PyDESeq2") %>%
  #dplyr::mutate(n_genes = "1000 genes") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, method, tool) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_adj_pydeseq <- res_adj_pydeseq %>% dplyr::filter(padj > 0.000000e+00,)
res_naive_pydeseq <- res_naive_pydeseq %>% dplyr::filter(padj > 2.041247e-79,)


#res1_naive_pydeseq <- res1_naive_pydeseq[-c(115),]
#res2_naive_pydeseq <- res2_naive_pydeseq %>% dplyr::filter(log2FoldChange < 1.5,)
#res2_adj_pydeseq <- res2_adj_pydeseq %>% dplyr::filter(log2FoldChange < 1.0 ,)

d_volcano_pydeseq <- rbind(res_naive_pydeseq, res_adj_pydeseq)

colnames(d_volcano_pydeseq) <- c("logFC", "padj", "isDE", "DEtype", "method", "tool")

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

res_naive_edge <- res_naive_edge[-c(1),]

d_volcano_edge <- rbind(res_naive_edge, res_adj_edge)
colnames(d_volcano_edge) <- c("logFC", "padj", "isDE", "DEtype", "method", "tool")
#colnames(d_volcano_pydeseq) <- colnames(d_volcano_edge)

common_genes <- intersect(rownames(d_volcano_edge), rownames(d_volcano_pydeseq))
d_volcano_pydeseq <- d_volcano_pydeseq[common_genes, ]
d_volcano_edge <- d_volcano_edge[common_genes, ]
d_volcano <- rbind(d_volcano_pydeseq, d_volcano_edge)

p_volcanos <-  d_volcano_pydeseq %>% 
  ggplot(mapping = aes(x=logFC, y=-log10(padj), col=DEtype)) +
  geom_point(size=.8) +
  theme_bw() +
  scale_color_manual(values = de_gene_colors) +
  #ggh4x::facet_nested(factor(tool, levels = c("PyDESeq2", "edgeR"))~factor(method, levels = c("CN naive", "CN aware")), scales ="free_y", independent = "y")+
  facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow=1, scale = "free")+
  ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col="") +
  ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  ggplot2::geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'bottom')
p_volcanos


# ROC curve #

res_naive <- read.csv("CN-aware-DGE/Python/results/res_CNnaive_test3.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/res_CNaware_test3.csv")

p1 <- matrix(res_naive$padj)
colnames(p1) <- "Pydeseq2_CN_naive"
p2 <- matrix(res_adj$padj)
colnames(p2) <- "Pydeseq2_CN_aware"
p <- cbind(p1,p2)

res <- res_adj%>% 
  mutate(truth = case_when(
    log2FoldChange >= 1.0 & padj <= 0.05 ~ "1",
    log2FoldChange <= -1.0 & padj <= 0.05 ~ "1",
    log2FoldChange < 1.0 | log2FoldChange > -1.0 & padj > 0.05 ~ "0")) 

truth <- as.vector(as.numeric(res$truth))
names(truth) <- res$geneID

rna_roc <- metaseqR2::diagplotRoc(truth = truth, p = p2, sig = 0.05, x = "fpr",
                                   y = "tpr", path = NULL, draw = TRUE)



### Impact of CN normalization on DGE analysis ###

# Stacked barplot #

res_naive <- read.csv("CN-aware-DGE/Python/results/res_CNnaive_test3.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/res_CNaware_test3.csv")

lfc_cut <- 1.0
pval_cut <- .05
#de_gene_colors <- c("Not significant" = "gray", "Down-reg" = "darkblue", "Up-reg"="darkred")

res_adj <- res_adj %>%
  #dplyr::filter(padj < pval_cut, abs(log2FoldChange) > lfc_cut) %>% 
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN aware") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, method) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_naive <- res_naive %>%
  #dplyr::filter(padj < pval_cut, abs(log2FoldChange) > lfc_cut) %>% 
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::mutate(method = "CN naive") %>% 
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype, method) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X")

res_adj<- res_adj %>% dplyr::filter(padj > 0.000000e+00,)

common_genes <- intersect(rownames(res_naive), rownames(res_adj))
res_naive <- res_naive[common_genes, ] 

p_data <- rbind(res_adj, res_naive)

p_data <- p_data %>%
  group_by(method, DEtype) %>%
  summarise(count = n()) %>%  # Count occurrences of each DEtype per method
  mutate(percentage = count / sum(count) * 100)

de_gene_colors <- c("Not significant" = "gray", "Down-reg" = "lightblue", "Up-reg"="pink")

p_bar <- ggplot2::ggplot(p_data, aes(x = method, y = percentage, fill = DEtype)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), size = 3)+
  scale_fill_manual(values = de_gene_colors) +
  labs(x = "", y = "") +
  theme_bw()+
  ggplot2::theme(legend.position = 'bottom')
p_bar


# Sankey diagram #

alluvial_data <- data.frame(
  CN_naive = c("Down-reg", "n.d.", "Up-reg"),
  CN_aware = c("Down-reg", "n.d.", "Up-reg"),
  freq_CN_naive = c(91, 2843, 246),  
  freq_CN_aware = c(670, 2456, 54)  
)

# Prepare the data in long format for ggalluvial
wide_alluvial_data <- data.frame(
  CN_naive = c("Down-reg", "n.d", "Up-reg", "n.d", "n.d", "n.d"),       
  CN_aware = c("n.d", "n.d", "n.d", "Down-reg", "n.d", "Up-reg"),
  frequency = c(alluvial_data$freq_CN_naive, alluvial_data$freq_CN_aware), 
  DE_group = c("Down-reg", "n.d", "Up-reg", "Down-reg", "n.d", "Up-reg") 
)

ggplot(data = wide_alluvial_data,
       aes(axis1 = CN_naive, axis2 = CN_aware, y = frequency)) +
  scale_x_discrete(limits = c("CN_naive", "CN_aware"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = DE_group)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "gray", "pink")) +
  theme(legend.position = "bottom") 


#library(networkD3)

#nodes <- data.frame(name = c("Down-reg", "n.d.", "Up-reg",
                             #"Down-reg", "n.d.", "Up-reg"))

#links <- data.frame(source = c(0, 0, 1, 1, 2, 2),  
                    #target = c(3, 3, 3, 4, 4, 5),  
                    #value  = c(91, 2843, 246, 670, 2456, 54))

#sankey <- sankeyNetwork(Links = links, Nodes = nodes,
                        #Source = "source", Target = "target", Value = "value",
                        #NodeID = "name", units = "Count",
                        #fontSize = 15, nodeWidth = 30, 
                        #sinksRight = TRUE)
#sankey
  
  
  



