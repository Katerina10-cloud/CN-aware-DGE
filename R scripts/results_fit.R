library(ggplot2)
library(tidyverse)

# Set input path
path <- "C:/Users/Documents/"
setwd(path)
# Import data
df <- read.csv(paste0(path, ''), row.names = 1)

#Loading the data
statRes_map_noCNV = read.csv('~/model_fit_Python/model_results/results_luad/results_3/statRes_map_noCNV.csv',header=TRUE)
statRes_map_CNV = read.csv('model_fit_Python/model_results/results_2/statRes_map_CNV.csv',header=TRUE)
metadata = read.csv('model_fit_Python/model_data/metadata.csv',header=TRUE)

save(res_noCNV, file = "~/model_fit_Python/model_results/results_luad/results_3/res_noCNV.Rdata")

colnames(statRes_map_noCNV)[3] <- "B1_1"
statRes_map_noCNV <- statRes_map_noCNV %>% select(1,3,7)
res_noCNV <- res_noCNV %>% remove_rownames %>% column_to_rownames(var="Row.names")
res_noCNV <- merge(statRes_map_noCNV, cnv, by = "row.names")

#getting DE genes
sum(statRes_map_noCNV$padj < 0.05 & statRes_map_noCNV$B1_1 >= 0.6, na.rm=TRUE) #up_regulated
sum(statRes_map_noCNV$padj < 0.05 & statRes_map_noCNV$B1_1 <= -0.6, na.rm=TRUE) #down-reg

sum(res_allG$B1_1 < 1.0 & res_allG$B1_1 > -1 & res_allG$padj > 0.05, na.rm=TRUE)
sum(res_noCNV$gene_group == "not significant" & res_noCNV$cn_group == "cn_amplification", na.rm=TRUE)
sum(deg$gene_group == "dosage-insensitive" & deg$cnv == "5", na.rm=TRUE)


# Volcano Plot
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
statRes_map_noCNV$diffexpressed <- "NO"
statRes_map_noCNV$diffexpressed[statRes_map_noCNV$B1_1 >= 0.6 & statRes_map_noCNV$padj < 0.05] <- "UP"
statRes_map_noCNV$diffexpressed[statRes_map_noCNV$B1_1 <= -0.6 & statRes_map_noCNV$padj < 0.05] <- "DOWN"

ggplot(data = statRes_map_noCNV, aes(x = B1_1, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "blue", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), labels = c("Down-regulated", "Not significant", "Up-regulated"))+
  scale_x_continuous(breaks = seq(-8, 12, 2))+
  labs(title="Tumor vs Normal DGE (LUAD, n_genes=23 080, DEG=32%)",x="Effect size (B1_1)")+
  theme_classic()



#deg_up <- subset(statRes_map_NOcnv, statRes_map_NOcnv$padj_1 <= 0.05 & statRes_map_NOcnv$B1_1 > 0.5)
#deg_down <- subset(statRes_map_NOcnv, statRes_map_NOcnv$padj_1 <= 0.05 & statRes_map_NOcnv$B1_1 < -0.5)
#deg <- rbind(deg_up, deg_down)

#colnames(deg_up)[2] <- "padj_1"
#deg_up_cnv <- deg_up_cnv %>% mutate(genes = "up")
#deg_up_merged <- deg_up_merged %>% select(1:4)
#deg_down_merged <- cbind(de_down, de_down_cnv)
#deg_up <- deg_up %>% select(1,3,7)
#deg_merged <- deg_merged %>% remove_rownames %>% column_to_rownames(var="Row.names")
#statRes_map_CNV <- statRes_map_CNV %>% remove_rownames %>% column_to_rownames(var="GeneID")
#deg_down_cnv <- deg_down_cnv[(rownames(deg_down_cnv) %in% rownames(deg_down)),] #delete rows by name
#deg_down_cnv <- deg_down_cnv %>% select(2,6)

deg_up_merged <- deg_up_merged %>% mutate(Difference = B1_2-B1_1)

save(deg_up_cnv, file = "model_fit_Python/model_results/deg_up_cnv.Rdata")
write.csv(rna_cnv_2test, file = "model_data/TCGA/lung_cancer/last_test/rna_cnv_2test.csv")


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


cnv_2 <- cnv %>% select(11:20)
luad_cnv <- luad_cnv/2
cnv <- cbind(cnv_1, cnv_2)


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


#library(gridExtra)
#library(grid)
#grid.arrange(p1, p2, ncol=2)

data_barplot_geneDosage <- data.frame(
  gene_group = rep(c("super-dosage", "d-sensitive", "d-insensitive"), each = 6),
  cn_group = rep(c("0", "1", "2", "3", "4", "5"), 3),
  number_of_genes = c(1, 23, 32, 26, 23, 65, 0, 1, 178, 483, 36, 796, 0, 0, 1, 399, 816, 600)
)

barplot <- ggplot(data_barplot_geneDosage, aes(fill = cn_group, y = number_of_genes, x = gene_group))+
  geom_bar(stat = "identity")+
  labs(x='Gene group', y='Frequency', title='CNV dosage effect on DGE, (LUAD, n_genes=3482)')+
  scale_fill_manual('Position', values=c('red', 'green', 'steelblue', 'violet', 'pink', 'coral2'))+
  #geom_text(aes(gene_group, label = number_of_genes), size = 3, position=position_dodge2(width=0.5))+
  theme_minimal()+
  #facet_wrap("cn_group", ncol=2)+
  guides(fill=guide_legend("CN group"))
barplot


  
