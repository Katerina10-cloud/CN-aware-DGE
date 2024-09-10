## Test edgeR ##

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "edgeR")
sapply(pkgs, require, character.only = TRUE)

rna <- read.csv("TCGA/lung_cancer/LUAD/rna_test_3.csv")
cnv <- read.csv("TCGA/lung_cancer/LUAD/cnv_test_3.csv")
metadata <- read.csv("TCGA/lung_cancer/LUAD/metadata_3.csv")

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
res_naive_edge <- edgeR::topTags(lrt, n=17667)$table

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
res_adj_edge <- edgeR::topTags(lrt_adj, n=3179)$table


rownames_idx <- match(rownames(res_naive_edge), rownames(res_adj_edge))
res_adj_edge <- res_adj_edge[rownames_idx,]




