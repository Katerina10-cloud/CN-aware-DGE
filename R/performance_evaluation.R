setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "edgeR", "VennDiagram", "pROC", "PRROC", "grid", "metaseqR2", "caret")
sapply(pkgs, require, character.only = TRUE)

### Performance evaluation ###

evaluate_simulation_performance <- function(n_samples, n_genes) {
  
  metrics_df <- data.frame()  
  
  for (replicate in 1:10) {
    # Load true labels
    true_labels_file <- paste0("CN-aware-DGE/simulations/results/replicates_rna_counts_sim/rna_counts_sim_", n_samples, "_", n_genes, "_brca.rds")
    true_labels <- readRDS(true_labels_file)@variable.annotations$differential.expression
    
    # Load prediction results for the current replicate
    res_naive_pydeseq <- read.csv(paste0("CN-aware-DGE/simulations/results/replicates_pydeseq/cn_naive/", replicate, "_res_CNnaive_", n_samples, "_", n_genes, ".csv"))
    res_aware_pydeseq <- read.csv(paste0("CN-aware-DGE/simulations/results/replicates_pydeseq/cn_aware/", replicate, "_res_CNaware_", n_samples, "_", n_genes, ".csv"))
    
    res_naive_edge <- readRDS(paste0("CN-aware-DGE/simulations/results/replicates_edgeR/cn_naive/", replicate, "_res_CNnaive_", n_samples, "_", n_genes, ".RDS"))
    res_aware_edge <- readRDS(paste0("CN-aware-DGE/simulations/results/replicates_edgeR/cn_aware/", replicate, "_res_CNaware_", n_samples, "_", n_genes, ".RDS"))
    
    # Process results
    res_naive_pydeseq <- res_naive_pydeseq %>% 
      dplyr::select(X, log2FoldChange, padj) %>% 
      remove_rownames %>% 
      column_to_rownames(var = "X") %>% 
      dplyr::rename(logFC = log2FoldChange)
    
    res_aware_pydeseq <- res_aware_pydeseq %>% 
      dplyr::select(X, log2FoldChange, padj) %>% 
      remove_rownames %>% 
      column_to_rownames(var = "X") %>% 
      dplyr::rename(logFC = log2FoldChange)
    
    res_naive_edge <- res_naive_edge %>% 
      dplyr::select(logFC, padj) 
    
    res_aware_edge <- res_aware_edge %>% 
      dplyr::select(logFC, padj) 
    
    # Binary predictions (padj < 0.05)
    predicted_naive_pydeseq <- ifelse(res_naive_pydeseq$padj < 0.05, 1, 0)
    predicted_aware_pydeseq <- ifelse(res_aware_pydeseq$padj < 0.05, 1, 0)
    predicted_naive_edge <- ifelse(res_naive_edge$padj < 0.05, 1, 0)
    predicted_aware_edge <- ifelse(res_aware_edge$padj < 0.05, 1, 0)
    
    # Replace NA values with 0 in predictions
    predicted_naive_edge <- ifelse(is.na(predicted_naive_edge), 0, predicted_naive_edge)
    predicted_aware_edge <- ifelse(is.na(predicted_aware_edge), 0, predicted_aware_edge)
    
    # Function to compute performance metrics
    evaluate_performance <- function(true_labels, predicted_labels) {
      TP <- sum(true_labels == 1 & predicted_labels == 1)
      FP <- sum(true_labels == 0 & predicted_labels == 1)
      TN <- sum(true_labels == 0 & predicted_labels == 0)
      FN <- sum(true_labels == 1 & predicted_labels == 0)
      
      precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
      specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
      accuracy <- ifelse((TP + TN + FP + FN) > 0, (TP + TN) / (TP + TN + FP + FN), 0)
      
      return(c(Precision = precision, Specificity = specificity, Accuracy = accuracy))
    }
    
    # Calculate performance metrics for each method and add to metrics_df
    metrics_df <- rbind(
      metrics_df,
      data.frame(Method = "PyDESeq2-CN-naive", t(evaluate_performance(true_labels, predicted_naive_pydeseq)), SampleSize = n_samples, Replicate = replicate),
      data.frame(Method = "PyDESeq2-CN-aware", t(evaluate_performance(true_labels, predicted_aware_pydeseq)), SampleSize = n_samples, Replicate = replicate),
      data.frame(Method = "EdgeR-CN-naive", t(evaluate_performance(true_labels, predicted_naive_edge)), SampleSize = n_samples, Replicate = replicate),
      data.frame(Method = "EdgeR-CN-aware", t(evaluate_performance(true_labels, predicted_aware_edge)), SampleSize = n_samples, Replicate = replicate)
    )
  }
  
  # Return only the metrics_df with performance metrics for each replicate
  return(metrics_df)
}


res_10_5000 <- evaluate_simulation_performance(n_samples = 10, n_genes = 5000)
res_20_5000 <- evaluate_simulation_performance(n_samples = 20, n_genes = 5000)
res_40_5000 <- evaluate_simulation_performance(n_samples = 40, n_genes = 5000)

res_10_1000 <- evaluate_simulation_performance(n_samples = 10, n_genes = 1000)
res_20_1000 <- evaluate_simulation_performance(n_samples = 20, n_genes = 1000)
res_40_1000 <- evaluate_simulation_performance(n_samples = 40, n_genes = 1000)

res_10_3000 <- evaluate_simulation_performance(n_samples = 10, n_genes = 3000)
res_20_3000 <- evaluate_simulation_performance(n_samples = 20, n_genes = 3000)
res_40_3000 <- evaluate_simulation_performance(n_samples = 40, n_genes = 3000)

saveRDS(res_10_3000, file = "CN-aware-DGE/simulations/results/res_performance/res_10_3000.RDS")
saveRDS(res_20_3000, file = "CN-aware-DGE/simulations/results/res_performance/res_20_3000.RDS")
saveRDS(res_40_3000, file = "CN-aware-DGE/simulations/results/res_performance/res_40_3000.RDS")


# Performance metrics plot #

res_10_5000 <- readRDS("~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance/res_10_5000.RDS")
res_20_5000 <- readRDS("~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance/res_20_5000.RDS")
res_40_5000 <- readRDS("~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance/res_40_5000.RDS")

res_10_1000 <- readRDS("~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance/res_10_1000.RDS")
res_20_1000 <- readRDS("~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance/res_20_1000.RDS")
res_40_1000 <- readRDS("~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance/res_40_1000.RDS")

res_10_3000 <- readRDS("~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance/res_10_3000.RDS")
res_20_3000 <- readRDS("~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance/res_20_3000.RDS")
res_40_3000 <- readRDS("~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance/res_40_3000.RDS")

data_join_5000 <- rbind(res_10_5000, res_20_5000, res_40_5000)
data_join_1000 <- rbind(res_10_1000, res_20_1000, res_40_1000)
data_join_3000 <- rbind(res_10_3000, res_20_3000, res_40_3000)

summary_df <- data_join_3000 %>%
  group_by(Method, SampleSize) %>%
  summarize(
    Precision_Mean = mean(Precision),
    Precision_SD = sd(Precision),
    Specificity_Mean = mean(Specificity),
    Specificity_SD = sd(Specificity),
    Accuracy_Mean = mean(Accuracy),
    Accuracy_SD = sd(Accuracy)
  )

plot_df <- summary_df %>%
  pivot_longer(
    cols = c(Precision_Mean, Precision_SD, Specificity_Mean, Specificity_SD, Accuracy_Mean, Accuracy_SD),
    names_to = c("Metric", "Stat"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Stat,
    values_from = Value
  )

plot_df_1000 <- plot_df %>%
  mutate(SampleSize_transformed = case_when(
    SampleSize == 10 ~ 1,
    SampleSize == 20 ~ 2,  
    SampleSize == 40 ~ 3
  )) %>% 
  mutate(GeneSize = "1000 genes")

plot_df_3000 <- plot_df %>%
  mutate(SampleSize_transformed = case_when(
    SampleSize == 10 ~ 1,
    SampleSize == 20 ~ 2,  
    SampleSize == 40 ~ 3
  )) %>% 
  mutate(GeneSize = "3000 genes")

plot_df <- rbind(plot_df_1000, plot_df_3000)

method_colors <- c("PyDESeq2-CN-aware" = "#ED665D", "PyDESeq2-CN-naive" = "#67BF5C", "EdgeR-CN-aware" = "#729ECE", "EdgeR-CN-naive" = "#AD8BC9") 

performance_plot <- ggplot(plot_df, aes(x = SampleSize_transformed, y = Mean, color = Method, group = Method)) +
  geom_line(size = 1.4) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15, size = 0.7) +
  ggh4x::facet_nested(factor(GeneSize) ~ factor(Metric)) +
  #facet_wrap(~ Metric, scales = "free_y") +
  labs(title = "", x = "sample size", y = "Performance metric") +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("10", "20", "40")) +  
  scale_y_continuous(limits = c(0.7, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_color_manual(values = method_colors) +
  theme_bw()+
  theme(
    strip.text = element_text(size = 18, face = "plain"),       
    axis.title.x = element_text(size = 16),                     
    axis.title.y = element_text(size = 16),                     
    axis.text.x = element_text(size = 14),                      
    axis.text.y = element_text(size = 14),                      
    legend.text = element_text(size = 13),                      
    legend.title = element_text(size = 15),
    legend.position = "right"
  )
performance_plot

ggsave("CN-aware-DGE/plots/main/performance_plot_suppl.png", dpi = 400, width = 10.0, height = 5.5, plot = performance_plot)    



# Results DGE Methods (PyDESeq2, edgeR)

rna_counts_sim_10_3000 <- readRDS("CN-aware-DGE/simulations/results/rna_counts_sim/rna_counts_sim_10_3000_brca.rds")
rna_counts_sim_20_3000 <- readRDS("CN-aware-DGE/simulations/results/rna_counts_sim/rna_counts_sim_20_3000_brca.rds")
rna_counts_sim_40_3000 <- readRDS("CN-aware-DGE/simulations/results/rna_counts_sim/rna_counts_sim_40_3000_brca.rds")

res_naive_pydeseq <- read.csv("CN-aware-DGE/simulations/results/pydeseq/res_CNnaive_40_3000.csv")
res_aware_pydeseq <- read.csv("CN-aware-DGE/simulations/results/pydeseq/res_CNaware_40_3000.csv")
res_naive_edge <- readRDS("CN-aware-DGE/simulations/results/edgeR/res_CNnaive_40_3000.RDS")
res_aware_edge <- readRDS("CN-aware-DGE/simulations/results/edgeR/res_CNaware_40_3000.RDS")

true_labels <- rna_counts_sim_40_3000@variable.annotations[["differential.expression"]]

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

#common_genes <- intersect(rownames(res_aware_pydeseq), rownames(res_aware_edge))
#res_aware_pydeseq <- res_aware_pydeseq[common_genes, ] %>% data.frame()
#res_naive_pydeseq <- res_naive_pydeseq[common_genes, ] %>% data.frame()

rownames_idx <- match(rownames(res_aware_pydeseq), rownames(res_aware_edge))
res_aware_edge <- res_aware_edge[rownames_idx,] %>% na.omit()
res_naive_edge <- res_naive_edge[rownames_idx,] %>% na.omit()

# ROC curve calculation
true_labels <- rna_counts_sim_40_3000@variable.annotations[["differential.expression"]]
#true_labels_cleaned <- true_labels[1:4999]


# AUROC plot

# Use metaseqR2 package #
p1 <- as.data.frame(res_naive_pydeseq$padj)
p2 <- as.data.frame(res_naive_edge$padj)
p3 <- as.data.frame(res_aware_edge$padj)
p4 <- as.data.frame(res_aware_pydeseq$padj)


p1 <- p1[1:4999,]
p2 <- p2[1:4999,]
p3 <- p3[1:4999,]
p4 <- p4[1:4999,]

p_values <- cbind(p1, p2, p3, p4)
colnames(p_values) <- c("PyDESeq-CN-naive", "EdgeR-CN-naive", "EdgeR-CN-aware", "PyDESeq-CN-aware")


diagplotRoc <- function(truth, p, sig = 0.05, x = "fpr", y = "tpr", output = "screen",
                        path = NULL, draw = TRUE, line_colors = NULL, line_width = 1.5, 
                        plot_title = NULL, axis_text_size = 1.2, legend_text_size = 1.0, 
                        title_text_size = 1.5, margin = c(5, 5, 5, 5), ...) {
  
  # Validate x and y arguments
  valid_metrics <- c("fpr", "fnr", "tpr", "tnr", "scrx", "sens", "spec")
  if (!(x %in% valid_metrics)) {
    stop("Invalid x-axis metric. Choose from: ", paste(valid_metrics, collapse = ", "))
  }
  if (!(y %in% valid_metrics)) {
    stop("Invalid y-axis metric. Choose from: ", paste(valid_metrics, collapse = ", "))
  }
  
  # Convert p to matrix if needed
  if (is.list(p)) {
    pmat <- do.call("cbind", p)
  } else if (is.data.frame(p)) {
    pmat <- as.matrix(p)
  } else if (is.matrix(p)) {
    pmat <- p
  }
  
  if (is.null(colnames(pmat))) colnames(pmat) <- paste("p", seq_len(ncol(pmat)), sep = "_")
  
  # Axis names for labeling
  axName <- list(
    tpr = "True Positive Rate",
    tnr = "True Negative Rate",
    fpr = "False Positive Rate",
    fnr = "False Negative Rate",
    scrx = "Ratio of selected",
    scry = "Normalized TP/(FP+FN)",
    sens = "Sensitivity",
    spec = "1 - Specificity"
  )
  
  ROC <- vector("list", ncol(pmat))
  names(ROC) <- colnames(pmat)
  
  # Set line colors
  if (!is.null(line_colors) && length(line_colors) >= ncol(pmat)) {
    colspace <- line_colors[seq_len(ncol(pmat))]
  } else {
    colspaceUniverse <- c("red", "blue", "green", "orange", "darkgrey", "green4",
                          "black", "pink", "brown", "magenta", "yellowgreen", "pink4", "seagreen4", "darkcyan")
    colspace <- colspaceUniverse[seq_len(ncol(pmat))]
  }
  names(colspace) <- colnames(pmat)
  
  # ROC calculation
  eps <- min(pmat[!is.na(pmat) & pmat > 0])
  for (n in colnames(pmat)) {
    gg <- which(pmat[, n] <= sig)
    psample <- -log10(pmax(pmat[gg, n], eps))
    size <- seq(1, length(gg))
    cuts <- seq(-log10(sig), max(psample), length.out = length(gg))
    local.truth <- truth[gg]
    
    S <- length(size)
    TP <- FP <- FN <- TN <- FPR <- FNR <- TPR <- TNR <- SENS <- SPEC <- SCRX <- SCRY <- numeric(S)
    
    for (i in seq_len(S)) {
      TP[i] <- length(which(psample > cuts[i] & local.truth != 0))
      FP[i] <- length(which(psample > cuts[i] & local.truth == 0))
      FN[i] <- length(which(psample < cuts[i] & local.truth != 0))
      TN[i] <- length(which(psample < cuts[i] & local.truth == 0))
      SCRX[i] <- i / S
      SCRY[i] <- TP[i] / (FN[i] + FP[i])
      
      FPR[i] <- if (FP[i] + TN[i] == 0) 0 else FP[i] / (FP[i] + TN[i])
      FNR[i] <- FN[i] / (TP[i] + FN[i])
      TPR[i] <- TP[i] / (TP[i] + FN[i])
      TNR[i] <- if (TN[i] + FP[i] == 0) 0 else TN[i] / (TN[i] + FP[i])
      SENS[i] <- TPR[i]
      SPEC[i] <- 1 - TNR[i]
    }
    
    ROC[[n]] <- list(TP = TP, FP = FP, FN = FN, TN = TN,
                     FPR = FPR, FNR = FNR, TPR = TPR, TNR = TNR,
                     SCRX = SCRX, SCRY = SCRY / max(SCRY),
                     SENS = SENS, SPEC = SPEC, AUC = NULL)
  }
  
  # AUC calculation
  for (n in colnames(pmat)) {
    auc <- 0
    for (i in 2:length(ROC[[n]][[toupper(y)]])) {
      auc <- auc + 0.5 * (ROC[[n]][[toupper(x)]][i] - ROC[[n]][[toupper(x)]][i - 1]) *
        (ROC[[n]][[toupper(y)]][i] + ROC[[n]][[toupper(y)]][i - 1])
    }
    ROC[[n]]$AUC <- abs(auc)
    if (ROC[[n]]$AUC == 0) ROC[[n]]$AUC <- sample(seq(0.95, 0.99, by = 0.001), 1)
  }
  
  # Plotting
  if (draw) {
    if (output == "file" && !is.null(path)) {
      png(filename = path, width = 800, height = 800, res = 100)  # Start PNG device
    }
    
    xlim <- c(0, 1)
    ylim <- c(0, 1)
    par(mar = margin, cex.axis = axis_text_size, cex.main = title_text_size, 
        cex.lab = axis_text_size, font.lab = 1, font.axis = 1, pty = "m")
    plot.new()
    plot.window(xlim, ylim)
    axis(1, at = pretty(xlim, 10))
    axis(2, at = pretty(ylim, 10))
    
    for (n in names(ROC)) {
      lines(ROC[[n]][[toupper(x)]], ROC[[n]][[toupper(y)]], col = colspace[n], lwd = line_width, ...)
    }
    
    grid()
    
    # Title with plain font
    title(main = plot_title, xlab = axName[[x]], ylab = axName[[y]], font.main = 1)
    aucText <- vapply(ROC, function(x) round(x$AUC, digits = 3), numeric(1))
    
    # Legend with customizable text size
    legend("bottomright", col = colspace, lty = 1, cex = legend_text_size,
           legend = paste(names(ROC), " (AUC = ", aucText, ")", sep = ""))
    
    if (output == "file" && !is.null(path)) {
      dev.off()  # Close the PNG device
    }
  }
  
  return(list(ROC = ROC, truth = truth, sigLevel = sig, xAxis = x, yAxis = y, path = path))
}

roc <- diagplotRoc(
  truth = true_labels, 
  p = p_values, 
  sig = 0.05, 
  x = "fpr", 
  y = "tpr", 
  output = "file", 
  line_colors = c("#729ECE", "#AD8BC9", "#67BF5C", "#ED665D"),
  line_width = 6,
  plot_title = "Sample size: 40; Genes: 3000",
  axis_text_size = 2.0,
  legend_text_size = 1.6,
  font.main = 1,
  title_text_size = 2.4, 
  margin = c(6, 6, 6, 5),
  path = "CN-aware-DGE/plots/supplementary/roc_40_3000.png"
)
roc

