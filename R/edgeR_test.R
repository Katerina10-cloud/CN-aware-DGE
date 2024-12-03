## Test edgeR ##

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "edgeR")
sapply(pkgs, require, character.only = TRUE)

#rna <- read.csv("TCGA/lung/LUAD/rna_test_1.csv")
#cnv <- read.csv("TCGA/lung/LUAD/cnv_test_1.csv")
#metadata <- read.csv("TCGA/lung/LUAD/metadata_1.csv")

rna <- read.csv("CN-aware-DGE/simulations/data/replicates/10_rna_join_100_5000.csv")
cnv <- read.csv("CN-aware-DGE/simulations/data/replicates/10_cn_join_100_5000.csv") %>% as.data.frame()
metadata <- read.csv("CN-aware-DGE/simulations/data/replicates/10_metadata_100_5000.csv")

cnv <- cnv %>% remove_rownames %>% column_to_rownames(var="X")
rna <- rna %>% remove_rownames %>% column_to_rownames(var="X")
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="X")


# CN naive
cn_naive <- function(rna, metadata) {
  design <- model.matrix(~1+condition, data=metadata)
  edger.obj <- edgeR::DGEList(rna)
  edger.obj <- edgeR::calcNormFactors(edger.obj, method="TMM")
  edger.obj <- edgeR::estimateDisp(edger.obj, design)
  fit <- edgeR::glmFit(edger.obj, design)
  lrt <- edgeR::glmLRT(fit, coef=2)
  return(lrt)
}

lrt <- cn_naive(rna, metadata)
res_naive_edge <- edgeR::topTags(lrt, n=Inf)$table

# CN aware #
design <- model.matrix(~1+condition, data=metadata)
edger.obj <- edgeR::DGEList(rna, group = metadata$condition)
edger.obj <- edgeR::calcNormFactors(edger.obj)
offset <- outer(rep(1,nrow(edger.obj)), getOffset(edger.obj)) + log(cnv)
offset <- offset %>% filter_all(all_vars(!is.infinite(.))) %>% as.matrix()
rna <- rna[ rownames(rna) %in% rownames(offset),]
cnv <- cnv[ rownames(cnv) %in% rownames(rna),]

fit_adj <- edgeR::glmFit(y=rna, design=design, offset=offset, dispersion = lrt[["dispersion"]])
lrt_adj <- edgeR::glmLRT(fit_adj, coef=2)
res_aware_edge <- edgeR::topTags(lrt_adj, n=Inf)$table

rownames_idx <- match(rownames(res_naive_edge), rownames(res_aware_edge))
res_aware_edge <- res_aware_edge[rownames_idx,]

res_naive_pydeseq <- read.csv("CN-aware-DGE/simulations/results/replicates_pydeseq/cn_naive/1_res_CNnaive_100_5000.csv")
#res_aware_pydeseq <- read.csv("CN-aware-DGE/simulations/results/replicates_pydeseq_5000/cn_aware/1_res_CNaware_10_5000.csv")
res_naive_pydeseq <- res_naive_pydeseq %>% remove_rownames %>% column_to_rownames(var="X") 

res_naive_edge <- res_naive_edge %>% dplyr::rename(padj = FDR)
res_aware_edge <- res_aware_edge %>% dplyr::rename(padj = FDR)

rownames_idx <- match(rownames(res_naive_pydeseq), rownames(res_naive_edge))
res_aware_edge <- res_aware_edge[rownames_idx,] %>% na.omit()
res_naive_edge <- res_naive_edge[rownames_idx,] %>% na.omit()

saveRDS(res_naive_edge, file = "CN-aware-DGE/simulations/results/replicates_edgeR/cn_naive/10_res_CNnaive_100_5000.RDS")
saveRDS(res_aware_edge, file = "CN-aware-DGE/simulations/results/replicates_edgeR/cn_aware/10_res_CNaware_100_5000.RDS")


# Volcano plot #

lfc_cut <- 0.0
pval_cut <- .05
de_gene_colors <- c("Not significant" = "gray", "Down-reg" = "#8491B4B2", "Up-reg"="#E7969C")

res_aware_edge <- res_aware_edge %>%
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
colnames(res_aware_edge) <- colnames(res_naive_edge)

d_volcano <- rbind(res_naive_edge, res_aware_edge)

p_volcanos <-  res_naive_edge %>% 
  ggplot(mapping = aes(x=logFC, y=-log10(padj), col=DEtype)) +
  geom_point(size=.6) +
  theme_bw() +
  scale_color_manual(values = de_gene_colors) +
  ggh4x::facet_nested(factor(tool, levels = c("PyDESeq2", "edgeR"))~factor(method, levels = c("CN naive", "CN aware")), scales ="free", independent = "y")+
  #ggh4x::facet_nested(factor(n_genes, levels = c("1000 genes", "3000 genes"))~factor(method, levels = c("CN naive", "CN aware")), scales ="free", independent = "y")+
  #ggh4x::facet_nested(factor(method, levels = c("CN naive", "CN aware"))~factor(tumor_type, levels = c("LUAD", "BRCA", "LIHC", "HNSC", "COAD")), scales ="free", independent = "y")+
  #facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow=1, scale = "free")+
  ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col="") +
  ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  ggplot2::geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'bottom')
p_volcanos


