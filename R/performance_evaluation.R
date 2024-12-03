setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "edgeR", "VennDiagram", "pROC", "PRROC", "grid", "metaseqR2", "caret")
sapply(pkgs, require, character.only = TRUE)
source("CN-aware-DGE/R/utils.R")

### Performance evaluation ###

n_genes_list <- c(1000, 3000, 5000)
n_samples_list <- c(10, 20, 40, 100)

for (n_genes in n_genes_list) {
  for (n_samples in n_samples_list) {
    res <- evaluate_simulation_performance(n_samples = n_samples, n_genes = n_genes)
    file_path <- sprintf("CN-aware-DGE/simulations/results/res_performance/res_%d_%d.RDS", n_samples, n_genes)
    saveRDS(res, file = file_path)
  }
}

load_results <- function(n_genes, sample_sizes, base_path) {
  files <- sprintf("%s/res_%d_%d.RDS", base_path, sample_sizes, n_genes)
  results <- lapply(files, readRDS)
  bind_rows(results)
}

# Load datasets
base_path <- "~/Documents/PhD_AI/CN-aware-DGE/simulations/results/res_performance"
sample_sizes <- c(10, 20, 40, 100)

data_1000 <- load_results(1000, sample_sizes, base_path)
data_5000 <- load_results(5000, sample_sizes, base_path)

rename_methods <- function(data) {
  data %>%
    mutate(Method = case_when(
      Method == "PyDESeq2-CN-naive" ~ "PyDESeq2",
      Method == "PyDESeq2-CN-aware" ~ "DeConveil",
      Method == "EdgeR-CN-naive" ~ "EdgeR",
      Method == "EdgeR-CN-aware" ~ "ABCD-DNA",
      TRUE ~ Method
    ))
}

data_1000 <- rename_methods(data_1000)
data_5000 <- rename_methods(data_5000)

# Combine and summarize the data
summarize_performance <- function(data) {
  data %>%
    group_by(Method, SampleSize) %>%
    summarize(
      Precision_Mean = mean(Precision),
      Precision_SD = sd(Precision),
      Specificity_Mean = mean(Specificity),
      Specificity_SD = sd(Specificity),
      Accuracy_Mean = mean(Accuracy),
      Accuracy_SD = sd(Accuracy),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(Precision_Mean, Precision_SD, Specificity_Mean, Specificity_SD, Accuracy_Mean, Accuracy_SD),
                 names_to = c("Metric", "Stat"),
                 names_sep = "_", values_to = "Value") %>%
    pivot_wider(names_from = Stat, values_from = Value)
}

summary_1000 <- summarize_performance(data_1000) %>% mutate(GeneSize = "1000 genes")
summary_5000 <- summarize_performance(data_5000) %>% mutate(GeneSize = "5000 genes")

plot_df <- bind_rows(summary_1000, summary_5000)


# Performance metrics plot #

method_colors <- c("DeConveil" = "#ED665D", "PyDESeq2" = "#67BF5C", "EdgeR" = "#729ECE", "ABCD-DNA" = "#AD8BC9")

performance_plot <- ggplot(plot_df, aes(x = SampleSize_transformed, y = Mean, color = Method, group = Method)) +
  geom_line(size = 1.4) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15, size = 0.7) +
  ggh4x::facet_nested(factor(GeneSize) ~ factor(Metric)) +
  #facet_wrap(~ Metric, scales = "free_y") +
  labs(title = "", x = "sample size", y = "Performance metric") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("10", "20", "40", "100")) +  
  scale_y_continuous(limits = c(0.7, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_color_manual(values = method_colors) +
  theme_bw()+
  theme(
    strip.text = element_text(size = 18, face = "plain", color = "black"),       
    axis.title.x = element_text(size = 16, color = "black"),                     
    axis.title.y = element_text(size = 16, color = "black"),                     
    axis.text.x = element_text(size = 14, color = "black"),                      
    axis.text.y = element_text(size = 14, color = "black"),                      
    legend.text = element_text(size = 13, color = "black"),                      
    legend.title = element_text(size = 15, color = "black"),
    legend.position = "right"
  )
performance_plot

ggsave("CN-aware-DGE/plots/main/performance_plot.png", dpi = 400, width = 10.0, height = 6.0, plot = performance_plot)    


# Results CN-aware & CN-naive methods

sample_sizes <- c(10, 20, 40, 100)
gene_counts <- c(1000, 3000, 5000)

load_rna_counts <- function(sample, gene) {
  file_path <- paste0("CN-aware-DGE/simulations/results/rna_counts_sim/rna_counts_sim_", sample, "_", gene, "_brca.rds")
  readRDS(file_path)
}

load_results <- function(method, sample, gene, extension = "csv") {
  file_dir <- ifelse(grepl("edge", method), "edgeR", "pydeseq")
  file_name <- paste0("res_", toupper(method), "_", sample, "_", gene, ".", extension)
  file_path <- paste0("CN-aware-DGE/simulations/results/", file_dir, "/", file_name)
  
  if (extension == "csv") {
    read.csv(file_path)
  } else {
    readRDS(file_path)
  }
}

rna_counts_list <- list()
results_list <- list()

# Load all combinations of RNA counts and results
for (gene in gene_counts) {
  for (sample in sample_sizes) {
    rna_counts_list[[paste0("rna_", sample, "_", gene)]] <- load_rna_counts(sample, gene)
    
    # Load PyDESeq2 (Naive & Aware) results
    results_list[[paste0("pydeseq_naive_", sample, "_", gene)]] <- load_results("CNnaive", sample, gene, "csv")
    results_list[[paste0("pydeseq_aware_", sample, "_", gene)]] <- load_results("CNaware", sample, gene, "csv")
    
    # Load EdgeR (Naive & Aware) results
    results_list[[paste0("edge_naive_", sample, "_", gene)]] <- load_results("CNnaive", sample, gene, "RDS")
    results_list[[paste0("edge_aware_", sample, "_", gene)]] <- load_results("CNaware", sample, gene, "RDS")
  }
}


process_results_pydeseq <- function(df) {
  df %>%
    dplyr::select(X, log2FoldChange, padj) %>%
    remove_rownames() %>%
    column_to_rownames(var = "X") %>%
    dplyr::rename(logFC = log2FoldChange) %>%
    na.omit()
}

process_results_edge <- function(df) {
  df %>%
    dplyr::select(logFC, FDR) %>%
    dplyr::rename(padj = FDR) %>%
    na.omit()
}

res_naive_pydeseq <- process_results_pydeseq(results_list[["pydeseq_naive_10_1000"]])
res_naive_edge <- process_results_edge(results_list[["edge_naive_10_1000"]])

res_aware_pydeseq <- process_results_pydeseq(results_list[["pydeseq_aware_10_1000"]])
res_aware_edge <- process_results_edge(results_list[["edge_naive_10_1000"]])

true_labels <- rna_counts_list[["rna_10_1000"]]@variable.annotations[["differential.expression"]]

rownames_idx <- match(rownames(res_aware_pydeseq), rownames(res_aware_edge))
res_aware_edge <- res_aware_edge[rownames_idx,] %>% na.omit()
res_naive_edge <- res_naive_edge[rownames_idx,] %>% na.omit()


# ROC curve calculation
true_labels <- rna_counts_list[["rna_10_1000"]]@variable.annotations[["differential.expression"]]

# AUROC plot

p1 <- as.data.frame(res_naive_pydeseq$padj)
p2 <- as.data.frame(res_naive_edge$padj)
p3 <- as.data.frame(res_aware_edge$padj)
p4 <- as.data.frame(res_aware_pydeseq$padj)


p1 <- p1[1:4993,]
p2 <- p2[1:4993,]
p3 <- p3[1:4993,]
p4 <- p4[1:4993,]

p_values <- cbind(p1, p2, p3, p4)
colnames(p_values) <- c("PyDESeq2", "EdgeR", "ABCD-DNA", "DeConveil")


roc <- diagplotRoc(
  truth = true_labels, 
  p = p_values, 
  sig = 0.05, 
  x = "fpr", 
  y = "tpr", 
  output = "file", 
  line_colors = c("#67BF5C", "#AD8BC9", "#729ECE", "#ED665D"),
  line_width = 6,
  plot_title = "Sample size: 100; Genes: 5000",
  axis_text_size = 2.0,
  legend_text_size = 1.6,
  font.main = 1,
  title_text_size = 2.4, 
  margin = c(6, 6, 6, 5),
  path = "CN-aware-DGE/plots/supplementary/roc_100_5000.png"
)
roc


# Radar chart - plot AUC values across methods #

library(fmsb)
library(scales)  # For alpha transparency in colors

data_1000 <- data.frame(
  Sample_Size = c("10 samples", "20 samples", "40 samples", "100 samples"),
  PyDESeq2 = c(0.779, 0.774, 0.78, 0.916),
  DeCNVeil = c(0.836, 0.904, 0.98, 0.988),
  EdgeR = c(0.735, 0.796, 0.80, 0.915),
  EdgeR_CN_aware = c(0.848, 0.904, 0.979, 0.982)
  )

data_3000 <- data.frame(
  Sample_Size = c("10 samples", "20 samples", "40 samples", "100 samples"),
  PyDESeq2 = c(0.774, 0.812, 0.819, 0.902),
  DeCNVeil = c(0.867, 0.936, 0.971, 0.991),
  EdgeR = c(0.762, 0.822, 0.835, 0.901),
  EdgeR_CN_aware = c(0.855, 0.929, 0.966, 0.989)
)

data_5000 <- data.frame(
  Sample_Size = c("10 samples", "20 samples", "40 samples", "100 samples"),
  PyDESeq2 = c(0.76, 0.809, 0.824, 0.897),
  DeCNVeil = c(0.87, 0.935, 0.97, 0.993),
  EdgeR = c(0.767, 0.80, 0.834, 0.908),
  EdgeR_CN_aware = c(0.871, 0.935, 0.967, 0.99)
)

prepare_radar_data <- function(data, min_val = 0.6, max_val = 1.0) {
  transposed <- as.data.frame(t(data[,-1]))  
  colnames(transposed) <- data$Sample_Size  
  
  # Create radar chart data with min and max limits
  radar_data <- as.data.frame(rbind(
    Max = rep(max_val, ncol(transposed)),  
    Min = rep(min_val, ncol(transposed)),
    transposed
  ))
  return(radar_data)
}

# Apply the function to all datasets
radar_data_1000 <- prepare_radar_data(data_1000)
radar_data_3000 <- prepare_radar_data(data_3000)
radar_data_5000 <- prepare_radar_data(data_5000)

colors <- c("#0072B2", "#6666FF", "#67BF5C", "#ED665D")
method_labels <- c("PyDESeq2", "EdgeR-CN-aware", "EdgeR", "DeCNVeil")

create_radarchart <- function(data, color = colors, 
                              vlabels = colnames(data), vlcex = 1.4,
                              caxislabels = c(0.6, 0.7, 0.8, 0.9, 1.0), 
                              title = "AUC Plot") {
  radarchart(
    data, axistype = 1,
    pcol = color,                      
    pfcol = alpha(color, 0.1),         
    plwd = 3,                          
    plty = 1,                          
    cglcol = "darkgray", cglty = 1,    
    cglwd = 1.0,                       
    axislabcol = "black",              
    vlcex = vlcex,                     
    vlabels = vlabels,                 
    caxislabels = caxislabels                           
  )
  title(main = title, font.main = 1, cex.main = 1.4)
}

# Sample size labels
vlabels <- c("n=10", "n=20", "n=40", "n=100")

# Plot for each radar dataset
create_radarchart(radar_data_1000, vlabels = vlabels, title = "AUC - 1000 genes")
create_radarchart(radar_data_3000, vlabels = vlabels, title = "AUC - 3000 genes")
create_radarchart(radar_data_5000, vlabels = vlabels, title = "AUC - 5000 genes")


legend("bottomleft", legend = method_labels, col = colors, 
       lty = 1, lwd = 2, bty = "n", cex = 1.2, 
       title = "Methods")
