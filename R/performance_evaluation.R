### Using a consensus gene list as a ground truth ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "edgeR", "VennDiagram", "pROC", "PRROC", "grid", "metaseqR2", "caret")
sapply(pkgs, require, character.only = TRUE)


# Results CN-Naive DGE Methods (PyDESeq2, edgeR)
res_naive_pydeseq <- read.csv("CN-aware-DGE/Python/results/HNSC/res_CNnaive_test.csv")
res_naive_edge <- readRDS("CN-aware-DGE/Python/results/HNSC/res_CNnaive_edge.RDS")

# Results CN-aware DGE Methods (PyDESeq2_CN, edgeR_CN)
res_aware_pydeseq <- read.csv("CN-aware-DGE/Python/results/HNSC/res_CNaware_test.csv")
res_aware_edge <- readRDS("CN-aware-DGE/Python/results/LUAD/res_CNaware_edge.RDS")

res_naive_pydeseq <- res_naive_pydeseq %>% dplyr::select(X,log2FoldChange, padj) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X") %>% 
  dplyr::rename(logFC = log2FoldChange)

res_aware_pydeseq <- res_aware_pydeseq %>% dplyr::select(X,log2FoldChange, padj) %>% 
  remove_rownames %>% 
  column_to_rownames(var="X") %>% 
  dplyr::rename(logFC = log2FoldChange)

res_naive_edge <- res_naive_edge %>% dplyr::select(logFC, FDR) %>% 
  dplyr::rename(padj = FDR)

res_aware_edge <- res_aware_edge %>% dplyr::select(logFC, FDR) %>% 
  dplyr::rename(padj = FDR)

common_genes <- intersect(rownames(res_aware_pydeseq), rownames(res_aware_edge))
res_aware_pydeseq <- res_aware_pydeseq[common_genes, ] %>% data.frame()
res_naive_pydeseq <- res_naive_pydeseq[common_genes, ] %>% data.frame()

rownames_idx <- match(rownames(res_aware_pydeseq), rownames(res_aware_edge))
res_aware_edge <- res_aware_edge[rownames_idx,] %>% na.omit()
res_naive_edge <- res_naive_edge[rownames_idx,] %>% na.omit()

# Extract DEG

lfc_cut <- 1.0
pval_cut <- .01

de_pydeseq_naive <- res_naive_pydeseq %>% dplyr::filter(padj < pval_cut ,)
de_pydeseq_naive_genes <- rownames(de_pydeseq_naive)

de_pydeseq_aware <-  res_aware_pydeseq %>% dplyr::filter(padj < pval_cut ,)
de_pydeseq_aware_genes <- rownames(de_pydeseq_aware)
  
de_edge_naive <- res_naive_edge %>% dplyr::filter(padj < pval_cut ,)
de_edge_naive_genes <- rownames(de_edge_naive)

de_edge_aware <- res_aware_edge %>% dplyr::filter(padj < pval_cut ,)
de_edge_aware_genes <- rownames(de_edge_aware)

  
# Find Consensus Across DGE Tools by identifying the overlap of DE genes across methods
de_results_list <- list(de_pydeseq_naive_genes, de_pydeseq_aware_genes, de_edge_naive_genes, de_edge_aware_genes)
consensus_genes <- Reduce(intersect, de_results_list)
#consensus_genes <- Reduce(union, list(de_pydeseq_naive_genes, de_pydeseq_aware_genes, de_edge_naive_genes, de_edge_aware_genes))

## Evaluate and compare methods ##

# Calculate overlap between CN-aware methods and consensus
edge_aware_overlap <- length(intersect(consensus_genes, de_edge_aware_genes)) / length(consensus_genes)
pydeseq_aware_overlap <- length(intersect(consensus_genes, de_pydeseq_aware_genes)) / length(consensus_genes)

# Calculate overlap between CN-naive methods (DESeq2 and edgeR) and consensus
pydeseq_naive_overlap <- length(intersect(consensus_genes, de_pydeseq_naive_genes)) / length(consensus_genes)
edge_naive_overlap <- length(intersect(consensus_genes, de_edge_naive_genes)) / length(consensus_genes)

# Calculate unique genes detected by CN-aware vs. CN-naive methods
unique_cn_aware <- setdiff(union(de_edge_aware_genes, de_pydeseq_aware_genes), union(de_edge_naive_genes, de_pydeseq_naive_genes))
unique_cn_naive <- setdiff(union(de_edge_naive_genes, de_pydeseq_naive_genes), union(de_edge_aware_genes, de_pydeseq_aware_genes))

# Report statistics
cat("Pydeseq CN-Aware Overlap with Consensus:", pydeseq_aware_overlap, "\n")
cat("Pydeseq naive Overlap with Consensus:", pydeseq_naive_overlap, "\n")
cat("edgeR CN-naive Overlap with Consensus:", edge_naive_overlap, "\n")
cat("edgeR CN-aware Overlap with Consensus:", edge_aware_overlap, "\n")
cat("Unique genes in CN-aware method:", length(unique_cn_aware), "\n")
cat("Unique genes in CN-naive methods:", length(unique_cn_naive), "\n")

# Create a Venn diagram for DE gene overlaps
venn.plot <- venn.diagram(
  x = list(
    "PyDESeq-CN-naive" = de_pydeseq_naive_genes,
    "PyDESeq-CN-aware" = de_pydeseq_aware_genes,
    "edgeR-CN-naive" = de_edge_naive_genes,
    "edgeR_CN-Aware" = de_edge_aware_genes
  ),
  category.names = c("PyDESeq-CN-naive", "PyDESeq-CN-aware", "edgeR-CN-naive", "edgeR-CN-aware"),
  fill = c("#FD744699", "#D5E4A299", "#D2AF8199", "#91D1C2B2"),  
  alpha = 0.5,  
  col = "black",  
  lty = "solid",  
  lwd = 2,  
  
  # Text style
  cex = 1.5,  
  fontface = "bold",
  fontfamily = "sans",  
  output = TRUE,
  filename = NULL
)
grid.draw(venn.plot)


## ROC & PR curves analysis ##

# ground_truth should be a binary vector, where 1 = DE gene, 0 = non-DE gene

proxy_ground_truth <- ifelse(rownames(res_aware_pydeseq) %in% consensus_genes, 1, 0) 

roc_pydeseq_naive <- roc(proxy_ground_truth, res_naive_pydeseq$padj)
roc_pydeseq_aware <- roc(proxy_ground_truth, res_aware_pydeseq$padj)
roc_edge_naive <- roc(proxy_ground_truth, res_naive_edge$padj)
roc_edge_aware <- roc(proxy_ground_truth, res_aware_edge$padj)

#plot(roc_edge_naive, main = "", col = "blue")

plot_roc_curves <- function(roc_list, method_names) {
  # Create an empty plot
  plot(roc_list[[1]], col = "#44AA99", main = "LUAD", lwd = 2, font.main = 1,
       xlab = "False Positive Rate", ylab = "True Positive Rate", cex.lab = 0.8)
  
  # Plot additional ROC curves
  colors <- c("#6699CC","#dcd300", "#CC6677")  
  
  for (i in 2:length(roc_list)) {
    plot(roc_list[[i]], col = colors[i - 1], add = TRUE, lwd = 2)
  }
  
  legend_text <- sapply(1:length(roc_list), function(i) {
    auc_value <- auc(roc_list[[i]]) 
    paste0(method_names[i], " (AUC = ", round(auc_value, 3), ")")
  })
  
  legend("bottomright", legend = legend_text, col = c("#44AA99", colors), lwd = 2, cex = 0.8)
}

roc_list <- list(roc_pydeseq_naive, roc_pydeseq_aware, roc_edge_naive, roc_edge_aware)
method_names <- c("PyDESeq-CN-naive", "PyDESeq-CN-aware", "edgeR-CN-naive", "edgeR-CN-aware")
plot_roc_curves(roc_list, method_names)


# Use metaseqR2 package #
p1 <- as.data.frame(res_naive_pydeseq$padj)
p2 <- as.data.frame(res_naive_edge$padj)
p3 <- as.data.frame(res_aware_pydeseq$padj)
p4 <- as.data.frame(res_aware_edge$padj)
p_values <- cbind(p1, p2, p3, p4)
colnames(p_values) <- c("PyDESeq-CN-naive", "EdgeR-CN-naive", "PyDESeq-CN-aware", "EdgeR-CN-aware")

roc <- diagplotRoc(truth = proxy_ground_truth, p = p_values, sig = 0.05, x = "fnr",
            y = "tnr", output = "x11", path = NULL,
            draw = TRUE)


# Calculate all metrics #

calculate_metrics <- function(predictions, ground_truth) {
  predictions <- as.numeric(predictions)
  ground_truth <- as.numeric(ground_truth)
  # Calculate confusion matrix components
  TP <- sum(predictions == 1 & ground_truth == 1)
  TN <- sum(predictions == 0 & ground_truth == 0)
  FP <- sum(predictions == 1 & ground_truth == 0)
  FN <- sum(predictions == 0 & ground_truth == 1)
  
  # Calculate metrics
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), NA)
  recall <- ifelse(TP + FN > 0, TP / (TP + FN), NA)  # Also called Sensitivity or True Positive Rate (TPR)
  specificity <- ifelse(TN + FP > 0, TN / (TN + FP), NA)  # Also called True Negative Rate (TNR)
  fpr <- ifelse(FP + TN > 0, FP / (FP + TN), NA)  # False Positive Rate
  f1_score <- ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), NA)
  fdr <- ifelse(TP + FP > 0, FP / (TP + FP), NA)  # False Discovery Rate (FDR)
  
  # MCC Calculation
  numerator <- (TP * TN) - (FP * FN)
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(denominator > 0, numerator / denominator, NA)
  
  # Return a list of all metrics
  return(list(
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    Specificity = specificity,
    FPR = fpr,
    FDR = fdr,
    F1_Score = f1_score,
    MCC = mcc
  ))
}

compare_methods_metrics <- function(methods_predictions, ground_truth, method_names) {
  results <- list()
  for (i in seq_along(methods_predictions)) {
    method_metrics <- calculate_metrics(methods_predictions[[i]], ground_truth)
    results[[method_names[i]]] <- method_metrics
  }
  return(results)
}

r_naive_pydeseq <- ifelse(rownames(res_naive_pydeseq) %in% rownames(de_pydeseq_naive), 1, 0) 
r_aware_pydeseq <- ifelse(rownames(res_aware_pydeseq) %in% rownames(de_pydeseq_aware), 1, 0) 
r_naive_edge <- ifelse(rownames(res_naive_edge) %in% rownames(de_edge_naive), 1, 0) 
r_aware_edge <- ifelse(rownames(res_aware_edge) %in% rownames(de_edge_aware), 1, 0) 

ground_truth <- proxy_ground_truth
methods_predictions <- list(r_naive_pydeseq, r_aware_pydeseq, r_naive_edge, r_aware_edge)
method_names <- c("Pydeseq-CN-naive", "Pydeseq-CN-aware", "EdgeR-CN-naive", "EdgeR-CN-aware")

metrics_results <- compare_methods_metrics(methods_predictions, ground_truth, method_names)
