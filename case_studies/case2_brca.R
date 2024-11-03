setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "ggrepel")
sapply(pkgs, require, character.only = TRUE)

# BRCA cancer #

clinical <- readRDS("TCGA/brca/clinical_filt.RDS")
rna_normal <- readRDS("TCGA/brca/rna_normal.RDS")
rna_tumor <- readRDS("TCGA/brca/rna_tumor.RDS")
cnv_tumor <- readRDS("TCGA/brca/cnv_tumor.RDS")


colnames(cnv_tumor) <- substr(colnames(cnv_tumor), 1, 12)
colnames(rna_normal) <- substr(colnames(rna_normal), 1, 12)
colnames(rna_tumor) <- substr(colnames(rna_tumor), 1, 12)

colnames(clinical) <- c("patientID", "gender", "age", "stage")
stage_1 <- clinical %>% dplyr::filter(stage %in% c("Stage I", "Stage IA", "Stage IIA"))
stage_2 <- clinical %>% dplyr::filter(stage %in% c("Stage IIB", "Stage IIIA", "Stage IIIB"))

#Stage_1
rna_tum_stage1 <- rna_tumor[,(colnames(rna_tumor) %in% stage_1$patientID)]
rna_norm_stage1 <- rna_normal[,(colnames(rna_normal) %in% stage_1$patientID)]
cnv_tumor_stage1 <- cnv_tumor[,(colnames(cnv_tumor) %in% stage_1$patientID)]

cnv_tumor_stage1 <- apply(cnv_tumor_stage1, 2, function(x) ifelse(x > 10, 10, x))
cnv_tumor_stage1 <- as.data.frame(cnv_tumor_stage1)

# Exclude genes with low expression in normal tissue #
low_expression_threshold <- 100
expression_summary <- data.frame(
  Gene = rownames(rna_normal),
  MeanExpression = rowMeans(rna_normal)
)

filtered_genes <- expression_summary %>%
  filter(MeanExpression > low_expression_threshold)

rna_normal <- rna_normal[filtered_genes$Gene, ]
rna_tumor <- rna_tumor[filtered_genes$Gene, ]
cnv_tumor <- as.data.frame(cnv_tumor)
cnv_tumor <- cnv_tumor[filtered_genes$Gene, ]

x <- colnames(rna_normal)
names(rna_normal) <- paste(x,"-11A")

x <- colnames(rna_tumor)
names(rna_tumor) <- paste(x,"-01A")

x <- colnames(cnv_tumor)
names(cnv_tumor) <- paste(x,"-01A")

rna <- cbind(rna_normal, rna_tumor)

hist(rowMeans(cnv_tumor),
     main = "GLM", 
     xlab = "CN state",
     ylab = "Proportion",
     col = "#E1DEFC",
     prob = TRUE,
     breaks = 15)

cnv_mean <- cnv_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_tumor)) %>% 
  dplyr::select(cnv_mean)

#Generate CN normal
cnv_normal <- matrix(2, nrow(cnv_tumor), 110)
rownames(cnv_normal) <- rownames(cnv_tumor)
cnv <- cbind(cnv_normal, cnv_tumor)
colnames(cnv) <- colnames(rna)
cnv <- cnv/2
colnames(rna) <- paste0("sample", 1:(ncol(rna)))
colnames(cnv) <- colnames(rna)

#Generate metadata#
metadata <- data.frame(patID = colnames(rna),
                       condition = rep(c("A", "B"), each = 110))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 
metadata$condition <- as.factor(metadata$condition)

write.csv(cnv, file = "TCGA/brca/case_study/cnv.csv", row.names = T)
write.csv(rna, file = "TCGA/brca/case_study/rna.csv", row.names = T)
write.csv(metadata, file = "TCGA/brca/case_study/metadata.csv", row.names = T)


### Downstaream analysis ###
res_naive <- read.csv("CN-aware-DGE/Python/results/case_studies/BRCA/res_CNnaive.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/case_studies/BRCA/res_CNaware.csv")
cnv <- read.csv("TCGA/brca/case_study/cnv.csv")
cancer_genes <- read.delim("TCGA/cancerGeneList.tsv")

# Enrichment analysis #

# Hallmark pathways
H_path_sensitive <- c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_E2F_TARGETS", "HALLMARK_ESTROGEN_RESPONSE_LATE",
                      "HALLMARK_MTORC1_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE")

H_path_compensated <- c("HALLMARK_MYOGENESIS", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_UV_RESPONSE_DN", 
                        "HALLMARK_ANGIOGENESIS", "HALLMARK_P53_PATHWAY", "HALLMARK_WNT_BETA_CATENIN_SIGNALING")

H_path_insensitive <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_KRAS_SIGNALING_UP",
                        "HALLMARK_ADIPOGENESIS", "HALLMARK_APICAL_JUNCTION")

# KEGG pathways #
K_path_sensitive <- c("Cell cycle", "DNA replication", "Biosynthesis of nucleotide sugars", "Biosynthesis of amino acids",
                      "Hippo signaling pathway")

K_path_compensated <- c("Wnt signaling pathway", "MAPK signaling pathway", "PPAR signaling pathway", 
                        "GnRH signaling pathway", "Parathyroid hormone synthesis, secretion and action")

K_path_insensitive <- c("Focal adhesion", "Calcium signaling pathway", "Motor proteins", "Regulation of lipolysis in adipocytes", 
                        "Oxytocin signaling pathway")


# GO pathways of prognostic genes
GO_path_compensated <- c("regulation of small GTPase mediated signal transduction", "Ras protein signal transduction",
                         "substrate adhesion-dependent cell spreading", "cell chemotaxis", "cellular response to salt")



