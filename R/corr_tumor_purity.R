setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "ppcor", "ggplot2", "readxl", "DESeq2")
sapply(pkgs, require, character.only = TRUE)

source("CN-aware-DGE/R/utils.R")

# RNA and tumor purity
cnv_filt <- readRDS("TCGA/lung/cnv_filt.RDS")

data_path <- "TCGA/lung/LUAD/rna_counts.RDS"
dataset_name <- "LUAD_rna"

rna <- rna_processing(dataset_name, data_path, cnv_filt)
rna_norm <- rna[[1]]
rna_tum <- rna[[2]]

quantile_90 <- apply(rna_tum, 1, function(x) quantile(x, 0.9))
filt_genes <- rna_tum[quantile_90 >= 30, ]

#select the most variable genes
gene_variances <- apply(filt_genes, 1, var)
top_5000 <- order(gene_variances, decreasing = TRUE)[1:5000]
var_genes <- filt_genes[top_5000, ]


# Selection of frequently amplified/deleted genes # 
#cnv_filt <- readRDS("TCGA/lung/cnv_filt.RDS")

# Thresholds
deleted_threshold <- 0.25  # 25% threshold for deletion (CN = 1)
amplified_threshold_1 <- 0.50  # 25% threshold for CN = 3, 4
amplified_threshold_2 <- 0.50  # 5% threshold for CN > 4


# Filter genes as "Deleted" (CN = 1 in at least 25% of samples)
deleted_genes <- rownames(cnv_filt)[rowMeans(cnv_filt == 1) >= deleted_threshold]

# Filter genes as "Amplified" (CN = 3 or 4 in at least 25% of samples OR CN > 4 in at least 5% of samples)
amplified_genes <- rownames(cnv_filt)[
  rowMeans(cnv_filt >= 3.0 & cnv_filt <= 4.0) >= amplified_threshold_1 |  # CN = 3 or 4 in at least 25% of samples
    rowMeans(cnv_filt > 4) >= amplified_threshold_2  # CN > 4 in at least 5% of samples
]

all_genes <- c(deleted_genes, amplified_genes)
cnv_filt <- cnv_filt[all_genes,]

common_genes <- intersect(rownames(cnv_filt), rownames(var_genes))
cnv_filt <- cnv_filt[common_genes, ]
var_genes <- var_genes[common_genes, ]

#Purity
purity <- read_excel("CN-aware-DGE/R/TCGA_tumor_purity.xlsx")
purity <- purity[ (purity$`Cancer type` %in% c("LUAD")) , ]

purity$`Sample ID` <- substr(purity$`Sample ID`, 1, 12)
purity <- purity[ (purity$`Sample ID` %in% colnames(cnv_filt)), ]
purity_filt <- purity %>% dplyr::select(`Sample ID`,`Cancer type`, `IHC`) %>% replace(is.na(.), 0.499)
purity_filt <- purity_filt[!duplicated(purity_filt$`Sample ID`), ] %>% 
  remove_rownames %>% 
  column_to_rownames(var="Sample ID")

x <- rownames(purity_filt)
rownames(purity_filt) <- paste(x,"-01A")

x <- colnames(cnv_filt)
colnames(cnv_filt) <- paste(x,"-01A")

colnames_idx <- match(colnames(rna_tum), rownames(purity_filt))
purity_filt <- purity_filt[colnames_idx,]


# Regress expression of frequently amplified/deleted genes against CN while controlling for tumor purity 
# (use partial correlation methods)

# Partial Spearman correlation
filt_genes_normalized <- var_genes %>% as.matrix() %>% DESeq2::vst()

partial_spearman <- function(gene_ind, expr_mat, cn_mat, purity) {
  expression_values <- expr_mat[gene_ind, ]
  cn_values <- cn_mat[gene_ind, ]
  partial_corr_result <- pcor.test(expression_values, cn_values, purity, method = "spearman")
  return(c(rho = partial_corr_result$estimate, p_value = partial_corr_result$p.value))
}

n_genes <- nrow(filt_genes_normalized)

if (nrow(filt_genes) != nrow(cnv_filt)) {
  stop("Error: The number of genes in the expression matrix and CN matrix do not match.")
}

results <- t(sapply(1:n_genes, partial_spearman, expr_mat = as.matrix(filt_genes_normalized), 
                    cn_mat = as.matrix(cnv_filt), purity = purity_filt$IHC))

colnames(results) <- c("p_Spearman", "pval")
rownames(results) <- rownames(filt_genes_normalized)

#results <- as.data.frame(results)
#results$pval_adj <- p.adjust(results$pval, method = "BH")

ggplot(results, aes(x=p_Spearman)) + 
  geom_density(color="darkblue")+
  labs(x ="Partial Spearman (ρ)", y="Density")+
  ggplot2::geom_vline(xintercept = 0.4, linetype = 'dashed', colour="darkred")+
  theme_bw()



