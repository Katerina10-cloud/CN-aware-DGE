### Impact of CNV on DGE analysis ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA/")
pkgs <- c("tidyverse", "biomaRt", "DESeq2", "ggplot2", "smplot2", "pheatmap")
sapply(pkgs, require, character.only = TRUE)


lihc_seg_annot <- readRDS("liver/cnv_seg_annot.RDS")
cnv_lihc <- readRDS("liver/cnv_tumor.RDS")
rna_tum <- readRDS("liver/rna_tumor.RDS")
rna_norm <- readRDS("liver/rna_normal.RDS")

# Plot CN-ratio across all chromosome

# Focus on individual chromosome
#luad_seg_filt <- luad_seg_annot %>% dplyr::filter(chr == "chr7")
#lusc_seg_filt <- lusc_seg_annot %>% dplyr::filter(chr == "chr8")

lihc_seg_filt <- lihc_seg_annot %>% dplyr::filter(chr == "chr8")
cnv_lihc <- cnv_lihc[ rownames(cnv_lihc) %in% lihc_seg_filt$GeneID , ]

cnv_lihc <- log2(cnv_lihc/2)

cnv_mean <- cnv_lihc %>% as.data.frame() %>% 
  dplyr::mutate(cnv_mean = rowMeans(cnv_lihc)) %>% 
  dplyr::select(cnv_mean)

colnames(cnv_mean) <- "cnv_mean"

gene_list <- rownames(cnv_mean)

rna_tum <- rna_tum[gene_list,] %>% dplyr::select(1:50)
rna_norm <- rna_norm[gene_list,] %>% dplyr::select(1:50)
rna <- cbind(rna_norm, rna_tum)


# Z-score trasformation
rna_log_normalized <- rna %>% as.matrix() %>% DESeq2::varianceStabilizingTransformation()
rna_zscore <- t(scale(t(rna_log_normalized)))
rna_zscore_normal <- rna_zscore %>% as.data.frame() %>% dplyr::select(1:50)
rna_zscore_tumor <- rna_zscore %>% as.data.frame() %>% dplyr::select(1:50) %>% as.matrix()

rna_zscore_normal <- rna_zscore_normal %>% 
  as.data.frame() %>%
  dplyr::mutate(zscore_mean = rowMeans(rna_zscore_normal)) %>% 
  dplyr::select(zscore_mean)

rna_zscore_tumor <- rna_zscore_tumor %>% 
  as.data.frame() %>% 
  dplyr::mutate(zscore_mean = rowMeans(rna_zscore_tumor)) %>% 
  dplyr::select(zscore_mean)

cnv_rna <- merge(cnv_mean, rna_zscore_tumor, by = "row.names")
colnames(cnv_rna) <- c("geneID", "cnv_mean", "zscore_mean")

# LUAD - select 20 genes 
gene_list_1 <- c("SNORA31", "SNORA20", "SNORD46", "AGK", "TAS2R38", "TRBV12-3", "TRBJ1-2", "ATG9B", "ATP6V0E2", "CDK5",
                 "CRYGN", "RHEB", "TMUB1", "ZNF398", "BET1", "DLX6-AS1", "ARPC1A", "COPS6", "PSMC2", "DPY19L1",
                 "FKBP14", "POLM", "MALSU1", "INTS1", "SEC61G", "PSPH", "VOPP1", "CCT6A", "NFE2L3", "HNRNPA2B1")

gene_list_2 <- c("EGFR", "PSPHP1", "ZNF713", "SUMF2", "VSTM2A", "CDC42P2", "ELDR", "LANCL2", "TUBBP6", "HAUS6P1",
                 "AGR2", "HDAC9", "AGMO", "C1GALT1", "COL28A1", "KIAA0895", "TBX20", "ZNRF2P1", "CRYZP1", "S100A11P1", 
                 "PMS2P7", "GNB2", "IFT22", "TRIM24", "TMEM60", "RPL41", "TRBV9", "TRBV8-2", "TRPV5", "PRSS2")

# LUSC - select 20 genes 
gene_list_1 <- c("MCPH1-AS1", "DEFB1", "NRG1", "CCAR2", "INTS10", "STC1", "FUT10", "UNC5D", "IKBKB", "IDO1",
                 "PENK", "XKR4", "CA8", "ZNF704", "TRIQK", "COL14A1", "NRBP2", "LINC00861", "GSDMD")

gene_list_2 <- c("DLGAP2", "ARHGEF10", "MCPH1", "PPP1R3B", "SFRP1", "UNC5D", "MTMR7", "NAT1", "TNFRSF10B", "BNIP3L",
                 "UBXN8", "RNF122", "KAT6A", "ATP6V1H", "UBE2V2", "RPS20", "FUT10", "LYPD2", "LACTB2", "UBE2W")

# BRCA - select 20 genes 
gene_list_1 <- c("ENTPD4", "RPL23AP55", "NEFL", "STC1", "DUSP4", "TEX15", "TTI2", "TPT1P8", "TPT1P8", "DKK4",
                 "HOOK3", "ADGRA2", "CEBPD", "SOX17", "CYP7B1", "CLVS1", "ADHFE1", "FABP4", "ANGPT1", "ENPP2", 
                 "PKHD1L1", "KLF10", "NECAB1")

gene_list_2 <- c("CCAR2", "RPL23AP55", "STMN4", "SNORA12", "GOLGA7", "NKX6-3", "PLAT", "RPL5P23", "HOOK3", "ADRB3",
                 "STAR", "PCMTD1", "OPRK1", "RPL37P6", "BHLHE22", "CA8", "SBSPON", "C8orf89", "CA13", "LY6K", "GPR20","NDRG1")

# LIHC - select 20 genes
gene_list_1 <- c("ADAMDEC1", "ESCO2", "ESCO2", "SGK3", "ADAM32", "ZNF703", "IDO1", "ANK1", "ZMAT4", "RNF170",
                 "HGSNAT", "EFCAB1", "ST18", "NKAIN3", "PMP2", "SNTG1", "DECR1", "LY6E", "KLF10", "ANKRD46", "PCMTD1", "LINC00968")

gene_list_2 <- c("C8orf48", "FAM86B1", "CCAR2", "PINX1", "ZNF395", "TDRP", "ADGRA2", "ZMAT4", "SLC20A2", "THAP1",
                 "SNTG1", "SNAI2", "LYN", "GGH", "DNAJC5B", "VCPIP1", "PI15", "TPD52", "RUNX1T1", "ARC", "SPAG1")


data_plot_1 <- cnv_rna[ cnv_rna$geneID %in% gene_list_1 , ]
data_plot_2 <- cnv_rna[ cnv_rna$geneID %in% gene_list_2 , ]

# Spearman correlation 

corr_plot2 <- ggplot2::ggplot(data_plot_2, aes(x = cnv_mean, y = zscore_mean)) +
  geom_point(shape = 21, fill = 'black', size = 1.4) +
  geom_smooth(method="lm",formula=y~x, color="red", fill="black", se=T) +
  smplot2::sm_statCorr(fit.params = list(color = "indianred"), separate_by = "\n", corr_method = 'spearman') +
  labs(x="CN ratio (log2)", y = "mRNA Z-score")+
  #geom_vline(xintercept = c(0.0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(0.0), col = "gray", linetype = 'dashed') +
  scale_y_continuous(limits = c(-0.3, 0.75), 
                     breaks = seq(-0.3, 0.75, by = 0.25),
                     expand = c(0, 0))+
  theme_classic()+
  ggtitle("LIHC - chr 8")+
  theme(plot.title = element_text(hjust = 0.5))

corr_plot2  
corr_plot1 + corr_plot2

# Mapping genes to chromosome locations #

ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
gene_list <- c(gene_list_1, gene_list_2)

gene_locations <- biomaRt::getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
  filters = "hgnc_symbol", 
  values = gene_list, 
  mart = ensembl
)
gene_locations <- gene_locations[!(gene_locations$band == ""), ]

gene_loc <- gene_locations %>% dplyr::select(hgnc_symbol,band) %>% 
  dplyr::rename(geneID = hgnc_symbol)


data_heatmap <- rbind(data_plot_1, data_plot_2)
colnames(data_heatmap) <- c("geneID", "CNV", "RNA")

data_heatmap <- data_heatmap[ data_heatmap$geneID %in% gene_locations$hgnc_symbol , ]
data_heatmap <- dplyr::left_join(data_heatmap, gene_loc, by = "geneID")

# Heatmap
data_long <- t(as.matrix(data_heatmap[, c(-1,-4)]))
colnames(data_long) <- data_heatmap$band 

heatmap <- pheatmap::pheatmap(data_long, 
         cluster_rows = F,   
         cluster_cols = F,   
         color = colorRampPalette(c("#8491B4B2", "white", "#F39B7FB2"))(100),  
         border_color = "white",  
         main = "LIHC - chr 8",  
         fontsize_row = 10,       
         fontsize_col = 8)
