# Exploratory analysis of CNV and RNAseq data (TCGA portal)

library(ggplot2)
library(tidyverse)
library(GenomicRanges)
library(biomaRt)
library(Homo.sapiens)
library(org.Hs.eg.db)
library(AnnotationHub)

#Replace by Checking Condition on Character Column
x <- luad_cnv_seg$chr
luad_cnv_seg$chr <- paste("chr",x, sep = "")

#Mapping genomic regions to gene symbols
ah <- AnnotationHub()
query(ah, c("Gencode", "gff", "human"))
ah["AH49556"]
gc <- ah[["AH49556"]]
gr <- GRanges(s12)

#Intersect my regions of interest with the annotation
overlaps <- findOverlaps(gr, gc)
genes <- extractList(gc$gene_name, as(overlaps, "List")) 
genes <- unstrsplit(unique(genes), sep = ";")
cnv_annotation <- paste(as.character(gr), genes)

#Splitting a large data file in multiple columns
cnv_annotation <- as.data.frame(cnv_annotation)
#colnames(cnv_s5_annotation)[3] <- "gene"
cnv_annotation$cnv_annotation <- gsub(' ', ':', cnv_annotation$cnv_annotation)
cnv_annotation <- data.frame(do.call("rbind", strsplit(cnv_annotation$cnv_annotation, split = ":", fixed = TRUE)))
colnames(cnv_annotation) <- c("chr", "range", "gene")
s12$range <- str_c(s12$start, "-", s12$end) #join multiple strings into a single string
s12_cnv_anno <- merge(s12, cnv_annotation, by = "range") %>% rename(chr.x = "chr") %>% as_tibble() %>% 
  dplyr::select(chr,start,end,range,gene,seg_mean,gene,sample)

#Split delimited strings in a column and insert as new rows
s12_cnv_anno <- s12_cnv_anno %>% 
  mutate(GeneID = strsplit(as.character(gene), ";")) %>% 
  unnest(GeneID)  

s12_cnv_anno <- s12_cnv_anno[!duplicated(s12_cnv_anno$GeneID), ] #remove dublicates
s12_cnv_anno <- s12_cnv_anno %>% dplyr::select(-range,-gene)

save(s12_cnv_anno, file = "s12_cnv_anno.Rdata")                          
