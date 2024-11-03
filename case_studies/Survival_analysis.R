### Survival analysis of TCGA-BRCA dataset ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "survival", "glmnet", "survminer", "survcomp", "DESeq2", "forestplot", "caret", "randomForestSRC")
sapply(pkgs, require, character.only = TRUE)

# Prepare data #
d_compensated <- readRDS("TCGA/brca/case_study/d_compensated_genes.RDS")
clinical_data <- readRDS("TCGA/brca/clinical_full.RDS")
rna <- read.csv("TCGA/brca/case_study/rna.csv")
rna_tumor <- readRDS("TCGA/brca/rna_tumor.RDS")
cnv_tumor <- readRDS("TCGA/brca/cnv_tumor.RDS")

colnames(rna_tumor) <- substr(colnames(rna_tumor), 1, 12)
colnames(cnv_tumor) <- substr(colnames(rna_tumor), 1, 12)
rna_tumor <- rna_tumor[(rownames(rna_tumor) %in% rna$X),]
cnv_tumor <- cnv_tumor[(rownames(cnv_tumor) %in% rna$X),]
clinical_data <- clinical_data[clinical_data$bcr_patient_barcode %in% colnames(rna_tumor),]
rna_tumor <- rna_tumor[,(colnames(rna_tumor) %in% clinical_data$bcr_patient_barcode)]

clinical_data <- clinical_data %>% dplyr::select(bcr_patient_barcode, days_to_death, days_to_last_followup)
clinical_data$event <- ifelse(!is.na(clinical_data$days_to_death), 1, 0)
clinical_data$time <- ifelse(!is.na(clinical_data$days_to_death), 
                            clinical_data$days_to_death, 
                            clinical_data$days_to_last_followup)

clinical_data <- clinical_data %>%
  dplyr::select(bcr_patient_barcode, time, event) %>%   
  drop_na() 
colnames(clinical_data) <- c("patientID", "time", "event")

rna_d_compensated_genes <- rna_tumor[(rownames(rna_tumor) %in% d_compensated$geneID),]

# GE normalization 
rna_log_norm<- rna_d_compensated_genes %>% as.matrix() %>% DESeq2::rlog()
cnv_tumor <- cnv_tumor[(rownames(cnv_tumor) %in% rownames(rna_log_norm)) ,]
cnv_tumor <- cnv_tumor/2
rna_log_norm <- rna_log_norm * cnv_tumor

# Reorder patients indexing
idx <- match(clinical_data$patientID, colnames(rna_log_norm))
rna_log_norm <- rna_log_norm[,idx]
rna_log_norm <- t(rna_log_norm)
clinical_data <- clinical_data %>% remove_rownames %>% column_to_rownames(var="patientID") 
data <- cbind(clinical_data, rna_log_norm)

# Split data into training and test set 
#set.seed(42)
# Create stratified split based on the "event" column (ensures balanced distribution)
train_index <- createDataPartition(data$event, p = 0.7, list = FALSE)
# Split into training and test sets
train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

# Check the distribution of events in both sets
table(train_data$event)
table(test_data$event)

# Extract the gene expression and clinical data for training and test sets
rna_train <- train_data %>% select(-time, -event)
clinical_train <- train_data %>% select(time, event)

rna_test <- test_data %>% select(-time, -event)
clinical_test <- test_data %>% select(time, event)

# Univariate Cox regression
surv_object <- survival::Surv(time = clinical_train$time, event = clinical_train$event)
cox_results <- data.frame(Gene = character(), p.value = numeric(), HR = numeric(), CI_lower = numeric(), CI_upper = numeric())

for (gene in colnames(rna_train)) {
  cox_model <- coxph(surv_object ~ rna_train[, gene], data = train_data)
  summary_cox <- summary(cox_model)
  
  cox_results <- rbind(cox_results, data.frame(
    Gene = gene,
    p.value = summary_cox$coefficients[,"Pr(>|z|)"],
    HR = summary_cox$coefficients[,"exp(coef)"],
    CI_lower = summary_cox$conf.int[,"lower .95"],
    CI_upper = summary_cox$conf.int[,"upper .95"]
  ))
}

significant_genes <- cox_results %>% filter(p.value < 0.05)

saveRDS(significant_genes, file = "TCGA/brca/case_study/significant_genes_cox.RDS")

# The Hazard Ratio (HR) and confidence intervals (CI) give an indication of the prognostic impact of each gene. 
# Genes with HR > 1 may be associated with poor prognosis, while HR < 1 suggests good prognosis.

# Feature Selection Using LASSO (Penalized Cox Model) #
# Build Gene Prognostic Signature Using Lasso Regression (selects the most important genes and avoids overfitting) 

# Prepare the matrix of gene expression (X) and the survival outcome (y)
X <- as.matrix(rna_train[, significant_genes$Gene])  
y <- surv_object  

# Perform Lasso Cox regression
lasso_model <- glmnet(X, y, family = "cox", alpha = 1)

# Use cross-validation to choose the best lambda
cv_model <- cv.glmnet(X, y, family = "cox", type.measure="C", alpha = 0.95, )
optimal_lambda <- cv_model$lambda.min

#plot(cv_model, 
     #xlab = "Log(λ)", ylab = "Partial Likelihood Deviance", 
     #col = "red", lwd = 2)
plot(cv_model)


# Get coefficients for the best lambda (selected genes)
lasso_coefs <- coef(lasso_model, s = optimal_lambda)
lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs))

# Extract genes with non-zero coefficients (i.e., selected genes)
selected_genes <- rownames(lasso_coefs_df)[lasso_coefs_df$`1` != 0]

## Selected genes form our prognostic signature
selected_genes

saveRDS(selected_genes, file = "TCGA/brca/case_study/prognostic_signature_aware.RDS")

# Train the Cox model using only the selected genes
cox_selected <- coxph(Surv(clinical_train$time, clinical_train$event) ~ ., data = rna_train[, selected_genes])


# Calculate the prognostic score for each patient
prognostic_score <- X[, selected_genes] %*% lasso_coefs[selected_genes,]
colnames(prognostic_score) <- c("progn_score")
clinical_train <- cbind(clinical_train, prognostic_score)

# Random Survival Forest (to account for nonlinear relationships between features and survival)
rsf_model <- rfsrc(Surv(time, event) ~ ., data = train_data, ntree = 100)
#var_importance <- vimp(rsf_model)

#png("var_importance.png", width = 800, height = 600)
#plot(var_importance)
#dev.off()

# Predict survival on the test set
rsf_predictions <- predict(rsf_model, newdata = test_data)


# Model validation #

# Cox model evaluation on test data
cox_model <- coxph(Surv(clinical_train$time, clinical_train$event) ~ ., data = rna_train)
cox_predictions <- predict(cox_model, newdata = rna_test, type = "lp")
c_index <- concordance(Surv(test_data$time, test_data$event) ~ cox_predictions)$concordance
print(paste("Concordance Index (Cox):", c_index))

#cindex <- concordance.index(clinical_data$progn_score, surv.time = clinical_data$time, surv.event = clinical_data$event)
#cindex$c.index  # C-index close to 1 indicates better predictive ability

# Make predictions on the test set using the selected genes
cox_predictions <- predict(cox_selected, newdata = rna_test[, selected_genes], type = "lp")
c_index <- concordance(Surv(clinical_test$time, clinical_test$event) ~ cox_predictions)$concordance
cat("Concordance Index with LASSO-selected genes: ", c_index, "\n")

# RSF model evaluation
rsf_c_index <- rsf_model$err.rate
print(paste("Concordance Index (RSF):", 1 - rsf_c_index))


#Stratify Patients and Perform Survival Analysis

#To assess the prognostic significance of the gene signature, 
#stratify patients into high-risk and low-risk groups based on the median prognostic score. 
#Then, use Kaplan-Meier survival curves to compare survival between these groups.

# Training data
clinical_train$risk_group <- ifelse(clinical_train$progn_score > median(clinical_train$progn_score), 
                                    "High risk", "Low risk")
                                    
# Test data
clinical_test$risk_group <- ifelse(cox_predictions>0,"High risk","Low risk")

surv_fit <- survfit(Surv(time, event) ~ risk_group, data = clinical_train)


ggsurvplot(surv_fit, data = clinical_train, pval = TRUE, 
           conf.int = TRUE, risk.table = TRUE, palette = c("#AD002AB2", "#00468BB2"),  
           title = "",  xlab = "Time (days)",  ylab = "Survival Probability",
           font.main = c(14, "bold", "black"),  
           font.x = c(12, "plain"),  font.y = c(12, "plain"),  
           font.tickslab = c(12, "plain"),  
           legend = "bottom",  
           legend.title = "Risk group",  
           legend.labs = c("High risk", "Low risk"),  
           risk.table.height = 0.25,  risk.table.y.text = TRUE,  
           risk.table.fontsize = 4,  risk.table.title = "Number at risk",  
           risk.table.col = "strata",  
           pval.coord = c(1000, 0.2),  pval.size = 5,  
           ggtheme = theme_classic2(),  
           risk.table.title.fontface = "bold"  
)



# Forest plot of univariate Cox proportional hazards regression analysis 
sel_genes_data <- significant_genes %>% dplyr::filter(Gene %in% c(selected_genes))

table_text <- cbind(
  c("", as.character(sel_genes_data$Gene)),  # Gene names (first column)
  c("Hazard Ratio", sprintf("%.2f", sel_genes_data$HR)),  # Hazard Ratio without CI (second column)
  c("p-value", sprintf("%.3f", sel_genes_data$p.value))  # p-value (third column)
)
colnames(table_text) <- c("", "Hazard ratio", "p_value")
table_text <- table_text[2:10,]

box_colors <- ifelse(sel_genes_data$HR > 1, "darkred", "darkblue")

forestplot(
  labeltext = table_text, 
  mean = c(NA, sel_genes_data$HR),
  lower = c(NA, sel_genes_data$CI_lower),
  upper = c(NA, sel_genes_data$CI_upper),
  new_page = TRUE,
  is.summary = c(TRUE, rep(FALSE, nrow(sel_genes_data))),
  #is.summary = rep(FALSE, nrow(table_text)),
  xlog = FALSE, 
  col = fpColors(box = "royalblue", lines = "darkblue", zero = "black"),
  xlab = "Hazard Ratio",
  txt_gp = fpTxtGp(
    label = gpar(
      fontface = c("plain"),  
      cex = 0.9
    ),
    ticks = gpar(cex = 0.8),
    xlab = gpar(cex = 1)
  ),
  boxsize = 0.1,
  zero = 1,
  line.margin = 0.5,
  graphwidth = unit(3, "inches"),
  ci.vertices = TRUE
)


# Boxplot - Concordance Index #
ci_cox_aware <- data.frame(CI = c(0.72, 0.56, 0.61, 0.65, 0.63, 0.55, 0.60, 0.67, 0.58, 0.77, 0.60, 0.55, 0.77, 0.68,
                                  0.70))
ci_cox_aware <- ci_cox_aware %>% dplyr::mutate(method = "CN-aware") %>% 
  dplyr::mutate(model = "Cox")

ci_cox_naive <- data.frame(CI = c(0.47, 0.38, 0.56, 0.38, 0.43, 0.62, 0.63, 0.49, 0.40, 0.48, 0.42, 0.51, 0.42, 0.41, 
                                  0.39))  
ci_cox_naive <- ci_cox_naive %>% dplyr::mutate(method = "CN-naive") %>% 
  dplyr::mutate(model = "Cox")

ci_rsf_aware <- data.frame(CI = c(0.54, 0.48, 0.62, 0.54, 0.46, 0.45, 0.61, 0.55, 0.52, 0.53, 0.56, 0.59, 0.56, 0.55,
                                  0.45))
ci_rsf_aware <- ci_rsf_aware %>% dplyr::mutate(method = "CN-aware") %>% 
  dplyr::mutate(model = "RSF")

ci_rsf_naive <- data.frame(CI = c(0.57, 0.54, 0.59, 0.50, 0.58, 0.57, 0.53, 0.65, 0.62, 0.49, 0.43, 0.63, 0.58, 0.64,
                                  0.60))
ci_rsf_naive <- ci_rsf_naive %>% dplyr::mutate(method = "CN-naive") %>% 
  dplyr::mutate(model = "RSF")

ci_data <- rbind(ci_cox_aware, ci_cox_naive, ci_rsf_aware, ci_rsf_naive)

method_colors <- c("#ADB17DFF", "#D49464FF")

bxp <- ggplot(ci_data, aes(x = model, y = CI, fill = method)) + 
  geom_boxplot(position = position_dodge())+
  labs(x="model", y = "Concordance Index", title = "")+
  labs(fill = "method")+
  theme_bw()+
  ggplot2::theme(legend.position = 'right')+
  #facet_wrap(~cancer_type)+
  scale_fill_manual(values = method_colors) +
  font("xy.text", size = 14, color = "black", face = "plain")+
  font("title", size = 16, color = "black")+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  #ggtitle("LUSC")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 14))
bxp
