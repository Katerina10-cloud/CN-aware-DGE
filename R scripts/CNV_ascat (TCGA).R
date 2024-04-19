# Exploratory analysis of CNV and RNAseq data (TCGA portal)

library(ggplot2)
library(tidyverse)
library(GenomicRanges)
library(biomaRt)
library(Homo.sapiens)
library("org.Hs.eg.db")
library(AnnotationHub)

#Replace by Checking Condition on Character Column
cnv_tumor$Chromosome[cnv_tumor$Chromosome == "Y"] <- "chrY"

#Mapping genomic regions to gene symbols
ah <- AnnotationHub()
query(ah, c("Gencode", "gff", "human"))
ah["AH49556"]
gc <- ah[["AH49556"]]
gr <- GRanges(cnv_tumor)

#Intersect my regions of interest with the annotation
overlaps <- findOverlaps(gr, gc)
genes <- extractList(gc$gene_name, as(overlaps, "List")) %>% 
  unstrsplit(unique(genes), ";")
cnv_annotation <- paste(as.character(gr), genes)

#Splitting a large data file in multiple columns
library(stringr)
library(tidyverse)
cnv_annotation <- as.data.frame(cnv_annotation)
colnames(cnv_s5_annotation)[3] <- "gene"
cnv_annotation$cnv_annotation <- gsub(' ', ':', cnv_annotation$cnv_annotation)
cnv_annotation <- data.frame(do.call("rbind", strsplit(cnv_annotation$cnv_annotation, split = ":", fixed = TRUE)))
cnv_tumor$combined <- str_c(cnv_tumor$Start, "-", cnv_tumor$End) #join multiple strings into a single string
cnv_anno <- merge(cnv_tumor, cnv_annotation, by = "combined") %>% 
  select(3,4,5,7,9) 

#Split delimited strings in a column and insert as new rows
library(tidyr)
library(tidyverse)
cnv_anno <- cnv_anno %>% 
  mutate(GeneID = strsplit(as.character(gene), ";")) %>% 
  unnest(GeneID)

s10_cnv <- s10_cnv %>% select(2,6) %>% na.omit()
colnames(s10_cnv)[2] <- "s10_tumor"
cnv <- cnv[!duplicated(cnv$GeneID), ] #remove dublicates
cnv <- merge(cnv, cnv_2, by = "GeneID")
