## Test edgeR ##

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("dplyr", "ggplot2", "edgeR")
sapply(pkgs, require, character.only = TRUE)

rna <- read.csv("TCGA/lung_cancer/LUAD/rna_test_4.csv")
cnv <- read.csv("TCGA/lung_cancer/LUAD/cnv_test_4.csv")
metadata <- read.csv("TCGA/lung_cancer/LUAD/metadata_4.csv")

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
res <- edgeR::topTags(lrt, n=3179)$table

# CN aware #
cn_aware <- function(rna, cnv, metadata) {
  design <- model.matrix(~1+condition, data=metadata)
  edger.obj <- dgeR::DGEList(rna, group = metadata$condition)
  edger.obj <- edgeR::calcNormFactors(edger.obj, method="TMM")
  offset <- outer( rep(1,nrow(data_obj)), getOffset(data_obj)) + log(cnv)
  offset <- offset %>% filter_all(all_vars(!is.infinite(.))) %>% as.matrix()
  rna <- rna[ rownames(rna) %in% rownames(offset),]
  cnv <- cnv[ rownames(cnv) %in% rownames(rna),]
  fit_adj <- edgeR::glmFit(y=luad_rna, design=design, offset=offset, dispersion = edger.obj[["tagwise.dispersion"]])
  lrt_adj <- edgeR::glmLRT(fit_adj, coef=2)
  return(lrt_adj)
}

lrt_adj <- cn_aware(rna, cnv, metadata)
res_adj <- edgeR::topTags(lrt_adj, n=3179)$table

rownames_idx <- match(rownames(res), rownames(res_adj))
res_adj <- res_adj[rownames_idx,]





