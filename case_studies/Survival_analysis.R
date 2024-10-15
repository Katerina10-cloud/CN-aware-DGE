setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "survival", "glmnet", "survminer", "survcomp", "DESeq2", "forestplot")
sapply(pkgs, require, character.only = TRUE)

# Prepare data #
d_compensated <- readRDS("TCGA/brca/case_study/d_compensated_genes.RDS")
clinical_data <- readRDS("TCGA/brca/clinical_full.RDS")
rna <- read.csv("TCGA/brca/case_study/rna.csv")
rna_tumor <- readRDS("TCGA/brca/rna_tumor.RDS")

colnames(rna_tumor) <- substr(colnames(rna_tumor), 1, 12)
rna_tumor <- rna_tumor[(rownames(rna_tumor) %in% rna$X),]
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

# Reorder patients indexing
idx <- match(clinical_data$patientID, colnames(rna_log_norm))
rna_log_norm <- rna_log_norm[,idx]
rna_log_norm <- t(rna_log_norm)
clinical_data <- clinical_data %>% remove_rownames %>% column_to_rownames(var="patientID") 
data <- cbind(clinical_data, rna_log_norm)

# Univariate Cox regression
surv_object <- survival::Surv(time = clinical_data$time, event = clinical_data$event)
cox_results <- data.frame(Gene = character(), p.value = numeric(), HR = numeric(), CI_lower = numeric(), CI_upper = numeric())

for (gene in colnames(rna_log_norm)) {
  cox_model <- coxph(surv_object ~ rna_log_norm[, gene], data = data)
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


# Build Gene Prognostic Signature Using Lasso Regression (selects the most important genes and avoids overfitting) #

# Prepare the matrix of gene expression (X) and the survival outcome (y)
X <- as.matrix(rna_log_norm[, significant_genes$Gene])  
y <- surv_object  

# Perform Lasso Cox regression
lasso_model <- glmnet(X, y, family = "cox", alpha = 1)

# Use cross-validation to choose the best lambda
cv_model <- cv.glmnet(X, y, family = "cox", alpha = 1)
optimal_lambda <- cv_model$lambda.min

plot(cv_model, 
     xlab = "Log(λ)", ylab = "Partial Likelihood Deviance", 
     col = "red", lwd = 2)

# Get coefficients for the best lambda (selected genes)
lasso_coefs <- coef(lasso_model, s = optimal_lambda)
lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs))

# Extract genes with non-zero coefficients (i.e., selected genes)
selected_genes <- rownames(lasso_coefs_df)[lasso_coefs_df$`1` != 0]

## Selected genes form our prognostic signature
selected_genes

saveRDS(selected_genes, file = "TCGA/brca/case_study/prognostic_signature.RDS")

# Calculate the prognostic score for each patient
prognostic_score <- X[, selected_genes] %*% lasso_coefs[selected_genes,]
colnames(prognostic_score) <- c("progn_score")
clinical_data <- cbind(clinical_data, prognostic_score)

#Stratify Patients and Perform Survival Analysis

#To assess the prognostic significance of the gene signature, 
#stratify patients into high-risk and low-risk groups based on the median prognostic score. 
#Then, use Kaplan-Meier survival curves to compare survival between these groups.

clinical_data$risk_group <- ifelse(clinical_data$progn_score > median(clinical_data$progn_score), 
                                   "High risk", "Low risk")
surv_fit <- survfit(Surv(time, event) ~ risk_group, data = clinical_data)


ggsurvplot(surv_fit, data = clinical_data, pval = TRUE, 
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


# Concordance index (C-index) for the model
cindex <- concordance.index(clinical_data$progn_score, surv.time = clinical_data$time, surv.event = clinical_data$event)
cindex$c.index  # C-index close to 1 indicates better predictive ability


# Forest plot of univariate Cox proportional hazards regression analysis 
sel_genes_data <- significant_genes %>% dplyr::filter(Gene %in% c(selected_genes))

table_text <- cbind(
  c("", as.character(sel_genes_data$Gene)),  # Gene names (first column)
  c("Hazard Ratio", sprintf("%.2f", sel_genes_data$HR)),  # Hazard Ratio without CI (second column)
  c("p-value", sprintf("%.3f", sel_genes_data$p.value))  # p-value (third column)
)
colnames(table_text) <- c("", "Hazard ratio", "p_value")
table_text <- table_text[2:9,]

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
  line.margin = 0.1,
  graphwidth = unit(2, "inches"),
  ci.vertices = TRUE
)




