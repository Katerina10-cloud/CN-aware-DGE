library(ggplot2)
library(tidyverse)

#Loading the data
statRes_map_noCNV = read.csv('~/model_fit_Python/model_results/results_1/statRes_map_noCNV.csv',header=TRUE)
statRes_map_CNV = read.csv('model_fit_Python/model_results/results_2/statRes_map_CNV.csv',header=TRUE)
metadata = read.csv('model_fit_Python/model_data/metadata.csv',header=TRUE)

save(res_allGenes, file = "~/model_fit_Python/model_results/results_2/res_allGenes.Rdata")

#getting DE genes
sum(statRes_map_noCNV$padj < 0.05 & statRes_map_noCNV$log2FoldChange > 0.5, na.rm=TRUE) #up_regulated
sum(res_allG$padj <= 0.05 & res_allG$B1_1 >= 1.0, na.rm=TRUE) #down-reg
sum(res_allG$B1_1 < 1.0 & res_allG$B1_1 > -1 & res_allG$padj > 0.05, na.rm=TRUE)

deg_up <- subset(statRes_map_noCNV, statRes_map_noCNV$padj < 0.05 & statRes_map_noCNV$log2FoldChange > 0.5)
deg_down <- subset(statRes_map_noCNV, statRes_map_noCNV$padj < 0.05 & statRes_map_noCNV$log2FoldChange < -0.5)



#colnames(deg_up)[2] <- "padj_1"
#deg_up_cnv <- deg_up_cnv %>% mutate(genes = "up")
#deg_up_merged <- deg_up_merged %>% select(1:4)
#deg_down_merged <- cbind(de_down, de_down_cnv)
#deg_up <- deg_up %>% select(1,3,7)
#deg_up <- deg_up %>% remove_rownames %>% column_to_rownames(var="GeneID")
#statRes_map_CNV <- statRes_map_CNV %>% remove_rownames %>% column_to_rownames(var="GeneID")
#deg_down_cnv <- deg_down_cnv[(rownames(deg_down_cnv) %in% rownames(deg_down)),] #delete rows by name
#deg_down_cnv <- deg_down_cnv %>% select(2,6)
deg_up_merged <- deg_up_merged %>% mutate(Difference = B1_2-B1_1)

save(deg_up_cnv, file = "model_fit_Python/model_results/deg_up_cnv.Rdata")
write.csv(rna_cnv_2test, file = "model_data/TCGA/lung_cancer/last_test/rna_cnv_2test.csv")

deg_merged <- deg_merged %>% mutate(difference = B1_2-B1_1)

#Data manipulation and cleaning
resFit_merged <- cbind(statRes_map_noCNV, statRes_map_CNV) %>% 
  select() %>% 
  mutate(Difference = B1_2-B1_1) %>% 
  relocate(B0_2, .before = B1_2) %>% 
  remove_rownames %>% column_to_rownames(var="GeneID")


#Calculate SD across rows
cnv_sd <- transform(cnv_filtered, sd=apply(cnv_filtered, 1, sd, na.rm=TRUE)) %>% 
  subset(cnv_sd, sd < 0.8) %>% 
  select(1,3,5:10)


rna_lusc %>% select(1,3,5:10,11,13,15:20) %>% 
  rna_lusc[(rownames(rna_lusc) %in% rownames(cnv_sd)),] #delete rows by name

cnv_2 <- cnv %>% select(11:20)
cnv_1 <- cnv_1/2
cnv <- cbind(cnv_1, cnv_2)

luad_cnv_tumor <- luad_cnv_tumor %>% remove_rownames %>% column_to_rownames(var="GeneID")

#write.csv(rna_lusc, file="model_fit_Python/model_data/test2/rna_lusc.csv")


deg_down <- deg_down %>% add_column(Group = "B1_1")
resFit_CNV <- resFit_CNV %>% add_column(Group = "B1_2")

resFit_plot <- rbind(resFit_noCNV, resFit_CNV)


#statRes_map_noCNV <- statRes_map_noCNV %>% mutate(cnv_type = ifelse(cnv_mean > 1.7 & cnv_mean < 2.5,"neutral", ifelse(cnv_mean>2.5, "gain", "loss")))
statRes_map_CNV <- statRes_map_CNV %>% mutate(cnv_type = case_when(cnv_mean > 1.7 & cnv_mean <= 2.5 ~ "neutral",
                          cnv_mean <= 1.7 ~ "loss",
                          cnv_mean > 2.5 & cnv_mean <= 4 ~ "gain",
                          cnv_mean > 4 ~ "amplifications"))
#Density plot

ggplot(resFit_merged, aes(Difference)) +
  geom_histogram(bins = 1000) +
  xlim(-4,4) +
  labs(title="B1 difference") +
  theme_classic()

ggplot(resFit_plot,aes(x=B1, fill = Group)) + 
  geom_histogram(data=subset(resFit_plot,Group == 'B1_1'),alpha =0.6, bins = 1000) +
  geom_histogram(data=subset(resFit_plot,Group == 'B1_2'), alpha =0.6, bins = 1000) +
  xlim(-7,12) +
  labs(title="B1 parameters density (MAP) test 2") +
  theme(legend.position="bottom", legend.box = "horizontal") + 
  theme_classic()

#Scatter plot
# scatter NV vs DP
scatter_plot = ggplot(statRes_map_CNV, aes(x = cnv_mean, y = Log2FC)) +
  geom_point(aes(color = cnv_type), size = .8) +
  theme(legend.position = 'bottom') +
  labs(title = "CNV correcrtion") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_hline(yintercept=1, linetype="dashed", color = "blue") +
  geom_hline(yintercept=-1, linetype="dashed", color = "blue") +
  geom_vline(xintercept = 1.8, color = "blue", size = 0.4) +
  geom_vline(xintercept = 2.5, color = "blue", size = 0.4) +
  geom_vline(xintercept = 4, color = "blue", size = 0.4) +
  theme_classic()
  #facet_wrap(~FILTER) +
scatter_plot


# Create new categorical column
resFit_2 <- resFit_2 %>%
  mutate(gene_type = case_when(B1_1 > 0.5 & padj_1 <= 0.05 ~ "up",
                               B1_1 < -0.5 & padj_2 <= 0.05 ~ "down",
                               TRUE ~ "ns")) %>% 
  count(gene_type) %>% 
  distinct(gene_type) %>% pull()

# Add colour, size and alpha (transparency) to volcano plot
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

resFit_2 %>%
  ggplot(aes(x = log2(B1_1),
             y = -log10(padj_1),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 1)),       
                     limits = c(-8, 10))  

#Volcano plot

library(EnhancedVolcano)
p2 <- EnhancedVolcano(resFit_CNV,
                lab = rownames(resFit_CNV),
                x = 'B1',
                y = 'padj',
                title = 'Tumor vs Control (LUSC, CNV)',
                subtitle = "",
                selectLab = c(''),
                legendLabSize = 8.0,
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 1.5,
                labSize = 2.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

library(gridExtra)
library(grid)

grid.arrange(p1, p2, ncol=2)
            
