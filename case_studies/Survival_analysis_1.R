### Survival analysis of TCGA-BRCA dataset ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "survival", "glmnet", "survminer", "survcomp", "DESeq2", "forestplot", "caret", "randomForestSRC", "gridExtra")
sapply(pkgs, require, character.only = TRUE)

# Prepare data #
d_compensated <- readRDS("TCGA/brca/case_study/d_compensated_genes_2.RDS")
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

rna_d_compensated_genes <- rna_tumor[(rownames(rna_tumor) %in% rownames(d_compensated)),]

# GE normalization 
rna_log_norm<- rna_d_compensated_genes %>% as.matrix() %>% DESeq2::varianceStabilizingTransformation()
cnv_tumor <- cnv_tumor[(rownames(cnv_tumor) %in% rownames(rna_log_norm)) ,]
cnv_tumor <- cnv_tumor[,(colnames(cnv_tumor) %in% colnames(rna_log_norm)) ]
cnv_tumor <- apply(cnv_tumor, 2, function(x) ifelse(x > 10, 10, x))
cnv_tumor <- cnv_tumor/2
rna_log_norm <- rna_log_norm * cnv_tumor

# Reorder patients indexing
idx <- match(clinical_data$patientID, colnames(rna_log_norm))
rna_log_norm <- rna_log_norm[,idx]
rna_log_norm <- t(rna_log_norm)
clinical_data <- clinical_data %>% remove_rownames %>% column_to_rownames(var="patientID") 
data <- cbind(clinical_data, rna_log_norm)

# Split the data into training (70%) and test (30%) sets
#train_index <- createDataPartition(data$event, p = 0.7, list = FALSE)
#train_data <- data[train_index, ]
#test_data  <- data[-train_index, ]

#train_data <- data[train_index, ]
#test_data  <- data[-train_index, ]

#rna_train <- train_data %>% select(-time, -event)
#clinical_train <- train_data %>% select(time, event)

#rna_test <- test_data %>% select(-time, -event)
#clinical_test <- test_data %>% select(time, event)

#table(train_data$event)
#table(test_data$event)

# Extract the gene expression and clinical data for training and test sets
#rna_train <- train_data %>% select(-time, -event)
#clinical_train <- train_data %>% select(time, event)

#rna_test <- test_data %>% select(-time, -event)
#clinical_test <- test_data %>% select(time, event)

rna <- data %>% select(-time, -event)
clinical <- data %>% select(time, event)


## Cox model ##

# Initial Cox model on Training set
surv_object <- survival::Surv(time = clinical$time, event = clinical$event)
cox_results <- data.frame(Gene = character(), p.value = numeric(), HR = numeric(), CI_lower = numeric(), CI_upper = numeric())

for (gene in colnames(rna)) {
  cox_model <- coxph(surv_object ~ rna[, gene], data = data)
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

# LASSO for further gene selection
X <- as.matrix(rna[, significant_genes$Gene])  
y <- surv_object  

lasso_model <- glmnet(X, y, family = "cox", alpha = 1)
cv_model <- cv.glmnet(X, y, family = "cox", type.measure="C", alpha = 1, )
optimal_lambda <- cv_model$lambda.min

plot(cv_model)

lasso_coefs <- coef(lasso_model, s = optimal_lambda)
lasso_coefs_df <- as.data.frame(as.matrix(lasso_coefs))
selected_genes <- rownames(lasso_coefs_df)[lasso_coefs_df$`1` != 0]
selected_genes

saveRDS(selected_genes, file = "TCGA/brca/case_study/prognostic_signature_aware.RDS")

# Train the Cox model using only the selected genes
cox_selected <- coxph(Surv(clinical$time, clinical$event) ~ ., data = rna[, selected_genes])

# Calculate the prognostic score for each patient
prognostic_score <- X[, selected_genes] %*% lasso_coefs[selected_genes,]
colnames(prognostic_score) <- c("progn_score")
clinical <- cbind(clinical, prognostic_score)


#Stratify Patients and Perform Survival Analysis

#To assess the prognostic significance of the gene signature, 
#stratify patients into high-risk and low-risk groups based on the median prognostic score. 
#Then, use Kaplan-Meier survival curves to compare survival between these groups.

# Training data
clinical$risk_group <- ifelse(clinical$progn_score > median(clinical$progn_score), 
                                    "High risk", "Low risk")
                                    
# Test data
#clinical_train$risk_group <- ifelse(cox_predictions>0,"High risk","Low risk")

surv_fit <- survfit(Surv(time, event) ~ risk_group, data = clinical)


surv_plot <- ggsurvplot(surv_fit, data = clinical, pval = TRUE, 
           conf.int = F, risk.table = TRUE, palette = c("#DD2A7B", "#515BD4"),  
           title = "",  xlab = "Time (days)",  ylab = "Survival Probability",
           font.main = c(16, "bold", "black"),  
           font.x = c(16, "plain"),  
           font.y = c(16, "plain"),  
           font.tickslab = c(16, "plain"),  
           legend = "bottom",  
           legend.title = "Risk group",  
           legend.labs = c("High risk", "Low risk"),  
           font.legend = c(16, "plain"),
           risk.table.height = 0.25,  
           risk.table.y.text = TRUE,  
           risk.table.fontsize = 6,  
           risk.table.title = "Number at risk",  
           risk.table.col = "strata",  
           pval.coord = c(1000, 0.2),  pval.size = 5,  
           ggtheme = theme_classic2()+
             theme(
             strip.text = element_text(size = 14, face = "plain") 
             ),
           risk.table.title.fontface = "bold"  
)

surv_plot$plot <- surv_plot$plot + 
  theme(
    legend.text = element_text(size = 16)) 
surv_plot$table <- surv_plot$table + 
  theme(
    text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "plain"))
surv_plot

#ggsave("CN-aware-DGE/case_studies/plots/brca/KM_survival.png", dpi = 400, width = 5.0, height = 10.0, plot = surv_plot)

# Forest plot of univariate Cox proportional hazards regression analysis 
sel_genes_data <- significant_genes %>% dplyr::filter(Gene %in% c(selected_genes))

table_text <- cbind(
  c("", as.character(sel_genes_data$Gene)),  # Gene names (first column)
  c("Hazard Ratio", sprintf("%.2f", sel_genes_data$HR)),  # Hazard Ratio without CI (second column)
  c("p-value", sprintf("%.3f", sel_genes_data$p.value))  # p-value (third column)
)
colnames(table_text) <- c("", "Hazard ratio", "p_value")
table_text <- table_text[2:10,]

box_colors <- ifelse(sel_genes_data$HR > 1, "deeppink4", "blue3")

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

# Using ggplot

sel_genes_data <- significant_genes %>% dplyr::filter(Gene %in% c(selected_genes))

plot_data <- sel_genes_data %>%
  mutate(
    gene = factor(Gene, levels = rev(Gene)),  # Reversing for top-to-bottom order
    is_summary = FALSE, # Summary flag for first row
    box_color = ifelse(HR > 1, "#D60C00FF", "blue3")  # Color based on HR values
  )

plot_data$is_summary[1] <- TRUE

library(RColorBrewer)
custom_colors <- c("category1" = "#D60C00FF", "category2" = "blue3")

forest_plot <- ggplot(plot_data, aes(x = HR, y = gene)) +
  geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper, color = box_color), 
                  size = 0.9, show.legend = TRUE) +
  #geom_point(data = subset(plot_data, is_summary == TRUE), aes(x = HR), 
             #shape = 23, fill = "royalblue", color = "royalblue", size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Hazard Ratio", y = NULL) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 15, face = "plain"),
    panel.grid.minor = element_blank(),  
    panel.grid.major.y = element_blank()) +
  scale_color_identity() +
  #scale_color_manual(values = custom_colors)+
  coord_cartesian(xlim = c(min(plot_data$CI_lower), max(plot_data$CI_upper))) +
  scale_x_continuous(trans = "identity")
forest_plot

# Create text table for gene labels, HR, and p-values (mimicking table_text)
table_text <- plot_data %>%
  transmute(
    Gene = as.character(Gene),
    `Hazard Ratio` = sprintf("%.2f", HR),
    `p-value` = sprintf("%.3f", p.value)
  )

# Create a table plot for the labels using gridExtra
table_plot <- tableGrob(table_text, rows = NULL, 
                        theme = ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0.1))))

# Combine forest plot and table plot
grid.arrange(table_plot, forest_plot, ncol = 2, widths = c(1, 2))

forest_plot
