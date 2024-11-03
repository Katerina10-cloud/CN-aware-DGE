setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")

pkgs <- c("tidyverse", "ggplot2", "ggrepel")
sapply(pkgs, require, character.only = TRUE)

# LUSC cancer #

res_naive <- read.csv("CN-aware-DGE/Python/results/case_studies/LUSC/res_CNnaive.csv")
res_adj <- read.csv("CN-aware-DGE/Python/results/case_studies/LUSC/res_CNaware.csv")
cnv <- read.csv("TCGA/lung/LUSC/cnv.csv")
cancer_genes <- read.delim("TCGA/cancerGeneList.tsv")
#hk_genes <- readRDS("TCGA/lung/housekeeping_genes_lung.RDS")

#hk_genes <- hk_genes %>% dplyr::select(Gene.Symbol) %>% 
  #rename(geneID = Gene.Symbol)

# Prepare data #

oncogenes <- cancer_genes %>% dplyr::filter(Is.Oncogene=="Yes") %>% 
  dplyr::select(Hugo.Symbol) %>% 
  dplyr::rename(geneID = Hugo.Symbol) %>% 
  dplyr::mutate(gene_type = "Oncogene")

tsg <- cancer_genes %>% dplyr::filter(Is.Tumor.Suppressor.Gene=="Yes") %>% 
  dplyr::select(Hugo.Symbol) %>% 
  dplyr::rename(geneID=Hugo.Symbol) %>% 
  dplyr::mutate(gene_type = "TSG")

cancer_genes <- rbind(oncogenes, tsg)

cnv <- cnv %>% 
  remove_rownames %>% 
  column_to_rownames(var="X") %>%  
  dplyr::select(111:220,)

cnv <- cnv * 2

cnv_mean <- cnv %>% 
  as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv)) %>% 
  dplyr::mutate(geneID = rownames(cnv)) %>% 
  dplyr::select(geneID, cnv_mean) 

res_naive <- res_naive %>% dplyr::select(X,log2FoldChange, padj) 
res_adj <- res_adj %>% dplyr::select(X,log2FoldChange, padj) 

lfc_cut <- 1.0
pval_cut <- .05

#de_gene_colors <- c("Not significant" = "gray", "Down-reg" = "royalblue3", "Up-reg"="hotpink3")  

res_adj <- res_adj %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype) 
  
res_naive <- res_naive %>%
  dplyr::mutate(isDE = (abs(log2FoldChange) >= lfc_cut) & (padj <= pval_cut)) %>%
  dplyr::mutate(DEtype = if_else(!isDE, "Not significant", if_else(log2FoldChange > 0, "Up-reg", "Down-reg"))) %>%
  dplyr::select(X,log2FoldChange, padj, isDE, DEtype)
 
colnames(res_naive) <- c("geneID", "log2FC_naive", "padj_naive", "isDE_naive", "DE_type_naive")
colnames(res_adj) <- c("geneID", "log2FC_adj", "padj_adj", "isDE_adj", "DE_type_adj")

res_joint <- dplyr::left_join(res_naive, res_adj, by = "geneID")


### Dosage-sensitive | Dosage-insensitive | Dosage-compensated ###

d_sensitive <- res_joint %>%
  dplyr::filter(DE_type_naive == "Up-reg" & DE_type_adj == "Not significant" | 
                  DE_type_naive == "Down-reg" & DE_type_adj == "Not significant") 

d_sensitive <- d_sensitive %>% left_join(cnv_mean, by = "geneID")

d_insensitive <- res_joint %>% 
  dplyr::filter(DE_type_naive == "Down-reg" & DE_type_adj == "Down-reg" |
                  DE_type_naive == "Up-reg" & DE_type_adj == "Up-reg") 

d_insensitive <- d_insensitive %>% left_join(cnv_mean, by = "geneID")


d_compensated <- res_joint %>% 
  dplyr::filter(DE_type_naive == "Not significant" & DE_type_adj == "Down-reg" | 
                  DE_type_naive == "Not significant" & DE_type_adj == "Up-reg") 

d_compensated <- d_compensated %>% left_join(cnv_mean, by = "geneID")

saveRDS(d_compensated, file = "TCGA/brca/case_study/d_compensated_genes.RDS")


d_sensitive_cancer_g <- d_sensitive[d_sensitive$geneID %in% cancer_genes$geneID ,]
d_sensitive_cancer_g <- d_sensitive_cancer_g[!duplicated(d_sensitive_cancer_g$geneID), ]

d_compensated_cancer_g <- d_compensated[d_compensated$geneID %in% cancer_genes$geneID ,]
d_compensated_cancer_g <- d_compensated_cancer_g[!duplicated(d_compensated_cancer_g$geneID), ]

d_insensitive_cancer_g <- d_insensitive[d_insensitive$geneID %in% cancer_genes$geneID ,]
d_insensitive_cancer_g <- d_insensitive_cancer_g[!duplicated(d_insensitive_cancer_g$geneID), ]


# CN-naive #
cn_naive_d_sensitive <- d_sensitive %>% dplyr::select(geneID, log2FC_naive, padj_naive, isDE_naive, DE_type_naive, cnv_mean) %>% 
  dplyr::mutate(method = "CN naive") %>%
  dplyr::mutate(gene_group = "Dosage-sensitive")
colnames(cn_naive_d_sensitive) <- c("geneID", "log2FC", "padj", "isDE", "DE_type", "cnv_mean", "method", "gene_group")

cn_naive_d_insensitive <- d_insensitive %>% dplyr::select(geneID, log2FC_naive, padj_naive, isDE_naive, DE_type_naive, cnv_mean) %>% 
  dplyr::mutate(method = "CN naive") %>%
  dplyr::mutate(gene_group = "Dosage-insensitive")
colnames(cn_naive_d_insensitive) <- c("geneID", "log2FC", "padj", "isDE", "DE_type", "cnv_mean", "method", "gene_group")

cn_naive_d_compensated <- d_compensated %>% dplyr::select(geneID, log2FC_naive, padj_naive, isDE_naive, DE_type_naive, cnv_mean) %>% 
  dplyr::mutate(method = "CN naive") %>%
  dplyr::mutate(gene_group = "Dosage-compensated")
colnames(cn_naive_d_compensated) <- c("geneID", "log2FC", "padj", "isDE", "DE_type", "cnv_mean", "method", "gene_group")

# CN-aware #
cn_aware_d_sensitive <- d_sensitive %>% dplyr::select(geneID, log2FC_adj, padj_adj, isDE_adj, DE_type_adj, cnv_mean) %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::mutate(gene_group = "Dosage-sensitive")
colnames(cn_aware_d_sensitive) <- c("geneID", "log2FC", "padj", "isDE", "DE_type", "cnv_mean","method", "gene_group")

cn_aware_d_insensitive <- d_insensitive %>% dplyr::select(geneID, log2FC_adj, padj_adj, isDE_adj, DE_type_adj, cnv_mean) %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::mutate(gene_group = "Dosage-insensitive")
colnames(cn_aware_d_insensitive) <- c("geneID", "log2FC", "padj", "isDE", "DE_type", "cnv_mean", "method", "gene_group")

cn_aware_d_compensated <- d_compensated %>% dplyr::select(geneID, log2FC_adj, padj_adj, isDE_adj, DE_type_adj, cnv_mean) %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::mutate(gene_group = "Dosage-compensated")
colnames(cn_aware_d_compensated) <- c("geneID", "log2FC", "padj", "isDE", "DE_type", "cnv_mean", "method", "gene_group")


# Volcano plot #
gene_group_colors <- c("Dosage-insensitive" = "#79AF9799", "Dosage-sensitive" = "#B2474599", "Dosage-compensated"="#3C5488B2")  
cnv_colors <- c("loss" = "#0073C299", "neutral" = "#86868699", "gain" = "#CFA127", "amplification" = "#DC0000B2")
cancer_g_compens <- d_compensated_cancer_g$geneID

v_plot_data <- rbind(cn_naive_d_sensitive, cn_naive_d_insensitive, cn_naive_d_compensated,
                    cn_aware_d_sensitive, cn_aware_d_insensitive, cn_aware_d_compensated)

v_plot_data <- v_plot_data %>% dplyr::filter(padj > 3.246376e-220 ,)

p_volcanos <-  v_plot_data %>% 
  ggplot(mapping = aes(x=log2FC, y=-log10(padj), col=gene_group)) +
  geom_point(size=1.0) +
  theme_bw() +
  scale_color_manual(values = gene_group_colors) +
  #ggh4x::facet_nested(factor(method, levels = c("CN naive", "CN aware"))~factor(tumor_type, levels = c("LUAD", "BRCA", "LIHC")), scales ="free", independent = "y")+
  facet_wrap(~factor(method, levels = c("CN naive", "CN aware")), nrow=1, scale = "free")+
  ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col="") +
  ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 'dashed') +
  ggplot2::geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  ggplot2::theme(legend.position = 'bottom')
p_volcanos


# Barplot #
cancer_g <- c(cancer_genes$geneID)
#hk_genes <- c(hk_genes$geneID)

label_genes <- function(gene_list, cancer_genes, hk_genes) {
  ifelse(gene_list %in% cancer_genes, "Cancer genes",
         ifelse(gene_list %in% hk_genes, "Housekeeping genes", "Other genes"))
}

#cn_aware_d_sensitive$gene_subcategory <- label_genes(cn_aware_d_sensitive$geneID, cancer_genes, hk_genes)
#cn_aware_d_insensitive$gene_subcategory <- label_genes(cn_aware_d_insensitive$geneID, cancer_genes, hk_genes)
#cn_aware_d_compensated$gene_subcategory <- label_genes(cn_aware_d_compensated$geneID, cancer_genes, hk_genes)

combined_data <- rbind(cn_aware_d_sensitive, cn_aware_d_insensitive, cn_aware_d_compensated)

combined_data <- combined_data %>% 
  dplyr::mutate(cnv_group = case_when(
    cnv_mean > 0.5 & cnv_mean <= 1.7  ~ "loss",
    cnv_mean > 1.7 & cnv_mean <= 2.4  ~ "neutral",
    cnv_mean > 2.4 & cnv_mean <=   4.3 ~ "gain",
    cnv_mean > 4.3 ~ "amplification"))

barplot_data <- combined_data %>%
  group_by(gene_group) %>%
  summarise(Count = n()) %>%
  mutate(total = sum(Count)) %>%
  mutate(percentage = (Count / total) * 100) %>%
  ungroup()

ggplot2::ggplot(barplot_data, aes(x = gene_group, y = percentage, fill = gene_group)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) + 
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4) +  
  scale_fill_manual(values = gene_group_colors) +  
  theme_classic() +  
  labs(y = "percentage of genes", x = "", title = "") +  
  theme(axis.text.x = element_text(size = 10),  
        axis.ticks.x = element_blank())+
  ggplot2::theme(legend.position = '')

ggplot2::ggplot(combined_data, aes(x = gene_group, fill = cnv_group)) +
  geom_bar(position = "stack", width = 0.6) + 
  #geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4) +  
  scale_fill_manual(values = cnv_colors) +  
  theme_classic() +  
  labs(y = "gene counts", x = "", title = "") +  
  theme(axis.text.x = element_text(size = 10),  
        axis.ticks.x = element_blank())+
  ggplot2::theme(legend.position = 'bottom')


#rm(barplot_data)


# Plot cancer genes | Dosage-sensitive | Dosage-compensated | CN-naive | CN-aware

d_sens_naive <- d_sensitive_cancer_g %>% dplyr::select(geneID, log2FC_naive) %>% 
  dplyr::mutate(method = "CN naive") %>%
  dplyr::mutate(gene_group = "Dosage-sensitive") %>% 
  dplyr::rename(log2FC = log2FC_naive)

d_insens_naive <- d_insensitive_cancer_g %>% dplyr::select(geneID, log2FC_naive) %>% 
  dplyr::mutate(method = "CN naive") %>%
  dplyr::mutate(gene_group = "Dosage-insensitive") %>% 
  dplyr::rename(log2FC = log2FC_naive)

d_comp_naive <- d_compensated_cancer_g %>% dplyr::select(geneID, log2FC_naive) %>% 
  dplyr::mutate(method = "CN naive") %>%
  dplyr::mutate(gene_group = "Dosage-compensated") %>% 
  dplyr::rename(log2FC = log2FC_naive)

d_sens_naive <- d_sens_naive %>% left_join(cancer_genes, by = "geneID") %>% 
  arrange(log2FC)
d_insens_naive <- d_insens_naive %>% left_join(cancer_genes, by = "geneID")
d_comp_naive <- d_comp_naive %>% left_join(cancer_genes, by = "geneID") %>% 
  arrange(log2FC)

d_insens_naive_30 <- d_insens_naive %>%
  dplyr::filter(log2FC < -2.13 | log2FC > 2.20) %>% 
  arrange(log2FC)
d_insens_naive_30 <- d_insens_naive_30[!duplicated(d_insens_naive_30$geneID), ]


d_sens_aware <- d_sensitive_cancer_g %>% dplyr::select(geneID, log2FC_adj) %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::mutate(gene_group = "Dosage-sensitive") %>% 
  dplyr::rename(log2FC = log2FC_adj)

d_insens_aware <- d_insensitive_cancer_g %>% dplyr::select(geneID, log2FC_adj) %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::mutate(gene_group = "Dosage-insensitive") %>% 
  dplyr::rename(log2FC = log2FC_adj)

d_comp_aware <- d_compensated_cancer_g %>% dplyr::select(geneID, log2FC_adj) %>% 
  dplyr::mutate(method = "CN aware") %>%
  dplyr::mutate(gene_group = "Dosage-compensated") %>% 
  dplyr::rename(log2FC = log2FC_adj)

d_insens_aware_30 <- d_insens_aware[d_insens_aware$geneID %in% d_insens_naive_30$geneID , ]

d_sens_aware <- d_sens_aware %>% left_join(cancer_genes, by = "geneID") %>% 
  arrange(log2FC)

d_insens_aware_30 <- d_insens_aware_30 %>% left_join(cancer_genes, by = "geneID") %>% 
  arrange(log2FC)
d_insens_aware_30 <- d_insens_aware_30[!duplicated(d_insens_aware_30$geneID), ]

d_comp_aware <- d_comp_aware %>% left_join(cancer_genes, by = "geneID") %>% 
  arrange(log2FC)

d_sens_aware <- d_sens_aware[!duplicated(d_sens_aware$geneID), ]
d_comp_aware <- d_comp_aware[!duplicated(d_comp_aware$geneID), ]
d_sens_naive <- d_sens_naive[!duplicated(d_sens_naive$geneID), ]
d_comp_naive <- d_comp_naive[!duplicated(d_comp_naive$geneID), ]

d_insens_aware_30 <- d_insens_aware[d_insens_aware$geneID %in% d_insens_naive_30$geneID , ]
d_insens_aware_30 <- d_insens_aware_30[!duplicated(d_insens_aware_30$geneID), ]



bargraph_data <- rbind(d_sens_naive, d_insens_naive_30, d_comp_naive, d_sens_aware, d_insens_aware_30,
                       d_comp_aware)

bargraph_data <- bargraph_data %>% 
  group_by(gene_group, method) %>%      
  arrange(log2FC, .by_group = TRUE) 

# Bar graph 
ggplot(bargraph_data, aes(x = reorder(geneID, log2FC), y = log2FC, fill = gene_type)) +
  geom_bar(stat = "identity") +  
  geom_hline(yintercept = 0, linetype = "dashed") +  
  theme_bw() +  
  labs(x = "Gene symbol", y = "log2FC") +  
  scale_fill_manual(values = c("Oncogene" = "#D43F3A99", "TSG" = "#357EBD99"))+
  ggh4x::facet_nested(factor(method, levels = c("CN naive", "CN aware"))~factor(gene_group, levels = c("Dosage-sensitive", "Dosage-compensated", "Dosage-insensitive")),
                      scales = "free", independent = "y")+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))


### Functional Enrichment analysis ###

cn_aware_d_sensitive <- cn_aware_d_sensitive %>% dplyr::select(geneID, log2FC)
cn_aware_d_compensated <- cn_aware_d_compensated %>% dplyr::select(geneID, log2FC)
cn_aware_d_insensitive <- cn_aware_d_insensitive %>% dplyr::select(geneID, log2FC)

pkgs <- c("ggplot2", "dplyr","tidyr","reactome.db", "fgsea", "org.Hs.eg.db", "data.table", "clusterProfiler", "enrichplot", "ggpubr", "msigdbr")
sapply(pkgs, require, character.only = TRUE)

# Prepare data
hs <- org.Hs.eg.db
my_symbols <- cn_aware_d_compensated$geneID
gene_list <- AnnotationDbi::select(hs,
                                   keys = my_symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")

# Overrapresentation analysis #
gene_l <- as.vector(gene_list$ENTREZID)

# MSigDb 
m_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") 
msig_H <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
msig_H <- enricher(gene_l, minGSSize = 10, maxGSSize = 500,
                                     pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = msig_H)

res_ora_H_sensitive <- msig_H@result %>% mutate(gene_group = "Dosage-sensitive")
res_ora_H_compensated <- msig_H@result %>% mutate(gene_group = "Dosage-compensated")
res_ora_H_insensitive <- msig_H@result %>% mutate(gene_group = "Dosage-insensitive")

# GO
oraGO <- enrichGO(gene = gene_l, ont = "BP", OrgDb = org.Hs.eg.db, 
                  keyType = "ENTREZID", pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 350)

res_ora_GO_sensitive <- oraGO@result %>% mutate(gene_group = "Dosage-sensitive")
res_ora_GO_compensated <- oraGO@result %>% mutate(gene_group = "Dosage-compensated")
res_ora_GO_insensitive <- oraGO@result %>% mutate(gene_group = "Dosage-insensitive")

# KEGG
kegg_organism = "hsa"
oraKEGG <- enrichKEGG(gene = gene_l, organism = kegg_organism,
                           minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")

res_ora_KEGG_sensitive <- oraKEGG@result %>% mutate(gene_group = "Dosage-sensitive")
res_ora_GO_compensated <- oraGO@result %>% mutate(gene_group = "Dosage-compensated")
res_ora_KEGG_insensitive <- oraKEGG@result %>% mutate(gene_group = "Dosage-insensitive")


# Visualize Gene - Biological term network

GO_path_compensated <- c("regulation of small GTPase mediated signal transduction", "Ras protein signal transduction",
                         "substrate adhesion-dependent cell spreading", "cell chemotaxis", "cellular response to salt")

prognostic_genes <- c("EPS8", "MYADM", "PRKCQ", "RASGRP2")

#ora_results <- ora_results %>% filter(pathway %in% GO_path_compensated & gene %in% prognostic_genes)
ora_results <- ora_results %>% filter(pathway %in% GO_path_compensated)

data <- res_GO_compensated %>%
  dplyr::select(Description, geneID_symbol) %>%   
  tidyr::separate_rows(geneID_symbol, sep = "/") %>% 
  dplyr::rename(term = Description, gene = geneID_symbol)

data$type <- ifelse(data$gene %in% data$gene, "gene", "term")

g <- graph_from_data_frame(data, directed = FALSE)
V(g)$type <- ifelse(V(g)$name %in% data$gene, "gene", "term")

# Define a custom label: label only biological terms and specific genes of interest
V(g)$label <- ifelse(V(g)$type == "term" | V(g)$name %in% prognostic_genes, V(g)$name, NA)


ggraph(g, layout = "fr") +  
  geom_edge_link(aes(edge_alpha = 0.2), color = "gray") +  
  geom_node_point(aes(color = type), size = 4) +  
  geom_node_text(aes(label = label, 
                     fontface = ifelse(type == "gene", "bold", "plain")),  
                 repel = TRUE, 
                 size = 4) +  
  labs(title = "") +
  theme_void() +
  scale_color_manual(values = c("gene" = "#7AA6DC99", "term" = "#E64B35B2"))+
  ggplot2::theme(legend.position = 'bottom')




# Hallmark pathways
H_path_sensitive <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_E2F_TARGETS", "HALLMARK_UV_RESPONSE_UP",
                      "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_MTORC1_SIGNALING")

H_path_compensated <- c("HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_ALLOGRAFT_REJECTION",
                        "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_MYOGENESIS",
                        "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

H_path_insensitive <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_INFLAMMATORY_RESPONSE",
                        "HALLMARK_ANGIOGENESIS", "HALLMARK_ESTROGEN_RESPONSE_LATE", "HALLMARK_HYPOXIA", "HALLMARK_IL2_STAT5_SIGNALING")

res_H_sensitive <- res_ora_H_sensitive %>% filter(Description %in% H_path_sensitive)
res_H_compensated <- res_ora_H_compensated %>% filter(Description %in% H_path_compensated)
res_H_insensitive <- res_ora_H_insensitive%>% filter(Description %in% H_path_insensitive)

res_H_sensitive <- res_H_sensitive %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_H_compensated <- res_H_compensated %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_H_insensitive <- res_H_insensitive %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

#res_H_sensitive$log_padjust <- -log10(res_H_sensitive$p.adjust)
#res_H_insensitive$log_padjust <- -log10(res_H_insensitive$p.adjust)
#res_H_compensated$log_padjust <- -log10(res_H_compensated$p.adjust)

p_data <- rbind(res_H_sensitive, res_H_compensated, res_H_insensitive)
p_data$Description <- fct_reorder(p_data$Description, p_data$gene_group)

#KEGG pathways
K_path_sensitive <- c("Ribosome biogenesis in eukaryotes", "Base excision repair", "Spliceosome",
                      "RNA degradation", "Biosynthesis of cofactors", "ATP-dependent chromatin remodeling", 
                      "Mucin type O-glycan biosynthesis")

K_path_compensated <- c("Hippo signaling pathway", "Motor proteins", "Cytokine-cytokine receptor interaction", "Pyruvate metabolism", 
                        "cGMP-PKG signaling pathway", "Focal adhesion", "Efferocytosis")

K_path_insensitive <- c("ECM-receptor interaction", "Complement and coagulation cascades", "Hematopoietic cell lineage",
                        "Arachidonic acid metabolism", "Wnt signaling pathway", "Hippo signaling pathway", "p53 signaling pathway", 
                        "Platelet activation", "Chemokine signaling pathway")

res_K_sensitive <- res_ora_KEGG_sensitive %>% filter(Description %in% K_path_sensitive)
res_K_compensated <- res_ora_KEGG_compensated %>% filter(Description %in% K_path_compensated)
res_K_insensitive <- res_ora_KEGG_insensitive%>% filter(Description %in% K_path_insensitive)

res_K_sensitive <- res_K_sensitive %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_K_compensated <- res_K_compensated %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

res_K_insensitive <- res_K_insensitive %>%
  separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = numerator / denominator)

#res_H_sensitive$log_padjust <- -log10(res_H_sensitive$p.adjust)
#res_H_insensitive$log_padjust <- -log10(res_H_insensitive$p.adjust)
#res_H_compensated$log_padjust <- -log10(res_H_compensated$p.adjust)

p_data <- rbind(res_K_sensitive, res_K_compensated, res_K_insensitive)
p_data$Description <- fct_reorder(p_data$Description, p_data$gene_group)


gene_group_colors = c("Dosage-sensitive" = "#AD002AB2", "Dosage-compensated" = "#00468BB2", "Dosage-insensitive" = "#1B1919B2")

p_gse <- ggplot(p_data, aes(x = GeneRatio_val, y = Description)) +
  geom_point(aes(size = Count, color = pvalue))+
  scale_color_gradient(low = "blue", high = "orange")+
  labs(x = "gene ratio", y = "", title = "")+
  facet_wrap(~gene_group)+
  theme(strip.text.x = element_text(size=10, color="black", face="bold.italic"))+
  theme_bw()+
  #theme(plot.title = element_text(size = 9, face = "bold"))+
  theme(legend.position = "bottom")+
  #scale_x_continuous(limits = c(0, 8))+
  theme(axis.text.y = element_text(color = gene_group_colors[p_data$gene_group]))+
  theme(legend.key.size = unit(0.6, "cm"), legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.spacing.y = unit(0.1, 'cm')) 
p_gse
