# Exploratory analysis of CNV and RNAseq data (TCGA portal)

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/TCGA")
pkgs <- c("tidyverse", "ggplot2", "GenomicRanges", "biomaRt", "Homo.sapience", "org.Hs.eg.db", "AnnotationHub")
sapply(pkgs, require, character.only = TRUE)

luad_seg <- readRDS("lung/LUAD/data/cnv_seg_annot.RDS")
brca_seg <- readRDS("brca/cnv_seg_annot.RDS")
lihc_seg <- readRDS("liver/cnv_seg.RDS")

#Replace by Checking Condition on Character Column
x <- cnv_seg$chr
cnv_seg$chr <- paste("chr",x, sep = "")

#Mapping genomic regions to gene symbols
ah <- AnnotationHub()
query(ah, c("Gencode", "gff", "human"))
ah["AH49556"]
gc <- ah[["AH49556"]]
gr <- GenomicRanges::GRanges(cnv_seg)

#Intersect my regions of interest with the annotation
overlaps <- IRanges::findOverlaps(gr, gc)
genes <- IRanges::extractList(gc$gene_name, as(overlaps, "List")) 
genes <- unstrsplit(unique(genes), sep = ";")
cnv_annotation <- paste(as.character(gr), genes)

#Splitting a large data file in multiple columns
cnv_annotation <- as.data.frame(cnv_annotation)
cnv_annotation$cnv_annotation <- gsub(' ', ':', cnv_annotation$cnv_annotation)
cnv_annotation <- data.frame(do.call("rbind", strsplit(cnv_annotation$cnv_annotation, split = ":", fixed = TRUE)))
colnames(cnv_annotation) <- c("chr", "range", "gene")
cnv_seg$range <- str_c(cnv_seg$start, "-", cnv_seg$end) #join multiple strings into a single string

cnv_anno <- merge(cnv_seg, cnv_annotation, by = "range") %>% 
  rename(chr.x = "chr") %>% 
  as_tibble() %>% 
  dplyr::select(chr,start,end,range,gene,seg_mean,gene,sample)


#Split delimited strings in a column and insert as new rows
cnv_anno <- cnv_anno %>% 
  dplyr::mutate(GeneID = strsplit(as.character(gene), ";")) %>% 
  tidyr::unnest(GeneID)  

#luad_cnv_anno <- luad_cnv_anno[!duplicated(luad_cnv_anno$GeneID), ] #remove dublicates
cnv_anno <- cnv_anno %>% dplyr::select(-gene)

saveRDS(cnv_anno, file = "lung/LUSC/cnv_seg_annot.RDS")                          
