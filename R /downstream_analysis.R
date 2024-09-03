setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "colorspace", "ggpubr", "gridExtra", "ggpointdensity")
sapply(pkgs, require, character.only = TRUE)

### Scatter plot (RNA LFC and CN relationship) ###

res <- read.csv("CN-aware-DGE/Python/results/res_CNnaive_test4.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/res_CNaware_test4.csv")
cnv_filt <- readRDS("TCGA/lung_cancer/LUAD/cnv_filt.RDS")

res <- res %>% select(logFC,FDR)
#res <- res %>% remove_rownames %>% column_to_rownames(var="X")
res_adj <- res_adj %>% select(logFC,FDR)
#res_adj <- res_adj %>% remove_rownames %>% column_to_rownames(var="X")

cnv_mean <- cnv_filt %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_filt)) %>% 
  dplyr::select(cnv_mean)

common_genes <- intersect(rownames(cnv_mean), rownames(res))
cnv_mean <- cnv_mean[common_genes, ] %>% data.frame()

p_naive <- merge(cnv_mean, res, by = "row.names")
p_adj <- merge(cnv_mean, res_adj, by = "row.names")

colnames(p_naive) <- c("geneID", "cnv_mean", "logFC", "FDR")
colnames(p_adj) <- colnames(p_naive)

p_naive <- plot_data %>% dplyr::filter(logFC > -2.5,)
p_adj <- p_adj %>% dplyr::filter(logFC > -4.0,)

d_scatter <- rbind(p_naive %>% dplyr::mutate(method = "no CN normalization"), 
                   p_aware %>% dplyr::mutate(method = "with CN normalization"))

scatter = ggplot(d_scatter, aes(x = log2(cnv_mean), y = logFC)) +
  geom_point(size=0.8) +
  #geom_smooth(method = "lm", color = "blue", se = T)+
  geom_smooth()+
  theme(legend.position = 'bottom') +
  labs(x ="log2(CN ratio)", y="log2FC") +
  facet_wrap(~method)+
  theme_classic()
scatter

#gridExtra::grid.arrange(scatter1, scatter2, nrow = 2)


# Log2FC comparison  - scatter plot #

lfc_naive <- p_naive %>% dplyr::select(logFC) %>% dplyr::rename(logFC=logFC)
lfc_aware <- p_adj %>% dplyr::select(logFC) %>% dplyr::rename(logFC_adj=logFC)
plot_lfc_edge <- merge(lfc_naive, lfc_aware, by = "row.names")
plot_lfc_pydeseq <- merge(lfc_naive, lfc_aware, by = "row.names")

#plot_lfc_edge <- subset(plot_lfc_edge, plot_lfc_edge$logFC < 2.0 & plot_lfc_edge$logFC > -2.0)

d_lfc <- rbind(plot_lfc_edge %>% dplyr::mutate(method = "edgeR"), 
                   plot_lfc_pydeseq %>% dplyr::mutate(method = "PyDESeq2"))


comparison_lfc <- ggplot(d_lfc, aes(x=logFC_adj, y=logFC)) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  xlab("log2FC with CN normalization") +
  ylab ("log2FC no CN normalization") +
  facet_wrap(~method)+
  theme_classic()+
  theme(legend.position="none")
comparison_lfc


# p-value comparison #
pval_naive <- p_naive %>% dplyr::select(FDR) %>% dplyr::rename(FDR=FDR)
pval_aware <- p_aware %>% dplyr::select(FDR) %>% dplyr::rename(FDR_adj = FDR)

plot_pval_edge <- merge(pval_naive, pval_aware, by = "row.names")
plot_pval_pydeseq <- merge(pval_naive, pval_aware, by = "row.names")

d_pval <- rbind(plot_pval_edge %>% dplyr::mutate(method = "edgeR"), 
               plot_pval_pydeseq %>% dplyr::mutate(method = "PyDESeq2"))

comparison_pval <- ggplot(d_pval, aes(x=-log10(FDR_adj), y=-log10(FDR))) + 
  geom_pointdensity(shape=20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed",
              color="darkred")+
  xlab("FDR with CN normalization") +
  ylab ("FDR no CN normalization") +
  facet_wrap(~method)+
  theme_classic()+
  theme(legend.position="none")
comparison_pval

gridExtra::grid.arrange(comparison_lfc, comparison_pval, nrow = 2)


### Volcano plot ###

res_naive <- read.csv("CN-aware-DGE/Python/results/res_CNnaive_test4.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/res_CNaware_test4.csv")
res_adj_shrink <- read.csv("CN-aware-DGE/Python/results/res_CNaware_sim_shrink.csv")

res_naive$diffexpressed <- "NO"
res_naive$diffexpressed[res_naive$logFC >= 1.0 & res_naive$FDR < 0.05] <- "UP"
res_naive$diffexpressed[res_naive$logFC <= -1.0 & res_naive$FDR < 0.05] <- "DOWN"
#res <- res %>% dplyr::filter(log2FoldChange < 2.0, log2FoldChange > -2.0,)

res_adj$diffexpressed <- "NO"
res_adj$diffexpressed[res_adj$logFC >= 1.0 & res_adj$FDR < 0.05] <- "UP"
res_adj$diffexpressed[res_adj$logFC <= -1.0 & res_adj$FDR < 0.05] <- "DOWN"
#res_adj <- res_adj %>% dplyr::filter(FDR > 2.180713e-65, logFC > -4.0,)

#res_adj_shrink$diffexpressed <- "NO"
#res_adj_shrink$diffexpressed[res_adj_shrink$log2FoldChange >= 1.0 & res_adj_shrink$padj < 0.05] <- "UP"
#res_adj_shrink$diffexpressed[res_adj_shrink$log2FoldChange <= -1.0 & res_adj_shrink$padj < 0.05] <- "DOWN"
#res_adj_shrink <- res_adj_shrink %>% dplyr::filter(log2FoldChange < 2.0, log2FoldChange > -1.5,)

d_volcano <- rbind(res_naive %>% dplyr::mutate(method = "no CN normalization"), 
                   res_adj %>% dplyr::mutate(method = "with CN normalization"))

p_volcano <- ggplot(data = d_volcano, aes(x = logFC, y = -log10(FDR), col = diffexpressed)) +
  geom_vline(xintercept = c(-1.0, 1.0), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("darkblue", "gray", "darkred"))+
  scale_x_continuous(breaks = seq(-8, 8, 1))+
  labs(x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 10, color = "black")+
  font("xlab", size = 10)+
  font("ylab", size = 10)+
  facet_wrap(~method)+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p_volcano