clustering_patients_cnv <- function(dataset_name, data_path) {
  if (dataset_name == "LUAD_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    pca_result <- prcomp(cnv_data_normalized, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    #set.seed(123)  
    kmeans_result <- kmeans(pca_data, centers = 3)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "LUSC_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    pca_result <- prcomp(cnv_data_normalized, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    #set.seed(123)  
    kmeans_result <- kmeans(pca_data, centers = 3)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "BRCA_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    pca_result <- prcomp(cnv_data_normalized, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    #set.seed(123)  
    kmeans_result <- kmeans(pca_data, centers = 5)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "LIHC_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    pca_result <- prcomp(cnv_data_normalized, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    #set.seed(123)  
    kmeans_result <- kmeans(pca_data, centers = 2)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "HNSC_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    non_constant_columns <- apply(cnv_data_normalized, 2, function(x) var(x) != 0)
    cnv_data_filtered <- cnv_data_normalized[, non_constant_columns]
    pca_result <- prcomp(cnv_data_filtered, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    kmeans_result <- kmeans(pca_data, centers = 2)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  } else if (dataset_name == "LUSC_cnv") {
    cnv_tumor <- readRDS(data_path)
    cnv_tumor <- as.matrix(t(cnv_tumor)) 
    cnv_data_normalized <- scale(cnv_tumor)
    non_constant_columns <- apply(cnv_data_normalized, 2, function(x) var(x) != 0)
    cnv_data_filtered <- cnv_data_normalized[, non_constant_columns]
    pca_result <- prcomp(cnv_data_filtered, scale. = TRUE)
    pca_data <- as.data.frame(pca_result$x[, 1:10])
    kmeans_result <- kmeans(pca_data, centers = 3)
    cnv_tumor <- as.data.frame(cnv_tumor)
    cnv_tumor$Cluster <- kmeans_result$cluster
    plot_clusters <- fviz_cluster(kmeans_result, data = pca_data, geom = "point", stand = FALSE) + 
      ggtitle("PCA of CNV Profiles")+
      theme_classic()
  }  else {
    stop("Dataset name not recognized")
  }
  return(list(cnv_tumor, plot_clusters))
}


rna_processing <- function(dataset_name, data_path, cnv_filt) {
  if (dataset_name == "LUAD_rna") {
    rna <- readRDS(data_path)
    colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% dplyr::select(1:45)
    rna_tum <- rna %>% dplyr::select(46:90)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_filt <- as.data.frame(cnv_filt)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_filt)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_filt)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")
  } else if (dataset_name == "LUSC_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:51)
    rna_tum <- rna %>% select(52:102)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_filt <- as.data.frame(cnv_filt)
    colnames(cnv_filt) <- stringr::str_sub(colnames(cnv_filt),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_filt)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_filt)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")
  } else if (dataset_name == "BRCA_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:110)
    rna_tum <- rna %>% select(111:220)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_filt <- as.data.frame(cnv_filt)
    colnames(cnv_filt) <- stringr::str_sub(colnames(cnv_filt),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_filt)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_filt)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")
  } else if (dataset_name == "LIHC_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:50)
    rna_tum <- rna %>% select(51:100)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_filt <- as.data.frame(cnv_filt)
    colnames(cnv_filt) <- stringr::str_sub(colnames(cnv_filt),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_filt)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_filt)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")  
  } else if (dataset_name == "HNSC_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:21)
    rna_tum <- rna %>% select(22:42)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_tumor <- as.data.frame(cnv_tumor)
    colnames(cnv_tumor) <- stringr::str_sub(colnames(cnv_tumor),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_tumor)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_tumor)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")  
  } else if (dataset_name == "COLON_rna") {
    rna <- readRDS(data_path)
    #colnames(rna) <- gsub(pattern = "\\.", replacement = "-", colnames(rna))
    rna_norm <- rna %>% select(1:12)
    rna_tum <- rna %>% select(13:24)
    colnames(rna_norm) <- stringr::str_sub(colnames(rna_norm),1,12)
    colnames(rna_tum) <- stringr::str_sub(colnames(rna_tum),1,12)
    cnv_tumor <- as.data.frame(cnv_tumor)
    colnames(cnv_tumor) <- stringr::str_sub(colnames(cnv_tumor),1,12)
    rna_tum <- rna_tum[,colnames(rna_tum) %in% colnames(cnv_tumor)]
    rna_norm <- rna_norm[,colnames(rna_norm) %in% colnames(cnv_tumor)]
    x <- colnames(rna_norm)
    names(rna_norm) <- paste(x,"-11A")
    x <- colnames(rna_tum)
    names(rna_tum) <- paste(x,"-01A")  
  }  else {
    stop("Dataset name not recognized")
  }
  return(list(rna_norm, rna_tum))
}


calculate_correlation <- function(cnv_values, expression_values) {
  cnv_values <- as.numeric(as.character(cnv_values))
  expression_values <- as.numeric(as.character(expression_values))
  valid_indices <- !is.na(cnv_values) & !is.na(expression_values)
  # Check if there are enough valid values to calculate correlation
  if (sum(valid_indices) > 1) {
    return(cor.test(cnv_values[valid_indices], expression_values[valid_indices], method = "pearson")$estimate)
  } else {
    return(NA)  # Return NA if not enough valid data
  }
}

calculate_p_value <- function(cnv_values, expression_values) {
  cnv_values <- as.numeric(as.character(cnv_values))
  expression_values <- as.numeric(as.character(expression_values))
  valid_indices <- !is.na(cnv_values) & !is.na(expression_values)
  # Check if there are enough valid values to calculate correlation
  if (sum(valid_indices) > 1) {
    return(cor.test(cnv_values[valid_indices], expression_values[valid_indices], method = "pearson")$p.value)
  } else {
    return(NA)  # Return NA if not enough valid data
  }
  
}

fit_linear_model <- function(cnv_values, expression_values) {
  cnv_values <- as.numeric(cnv_values)
  expression_values <- as.numeric(expression_values)
  
  valid_indices <- !is.na(cnv_values) & !is.na(expression_values)
  
  if (sum(valid_indices) > 1) {
    model <- lm(expression_values[valid_indices] ~ cnv_values[valid_indices])
    return(c(slope = coef(model)[2], p.value = summary(model)$coefficients[2, 4]))
  } else {
    return(c(slope = NA, p.value = NA))
  }
}

