
# Copy-number-aware Differential Expression Simulation

```{r setup, include=FALSE}
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
```
## Introduction

This document gives some extra details regarding Differential Gene Expression (DGE) simulation informed by Copy Number (CN) signal.
RNA counts are simulated using 'compcodeR' simulator function `generateSyntheticData`, in particular, normal gene expression profiles of mean and dispersion pairs estimated from TCGA-BRCA dataset using DESeq2. 
Copy Number data are simulated using TCGA CN gene specific cancer data (absolute numbers) + sampling method.
We construct a simulated dataset of Negative Binomial data from two condition (tumor vs normal) by inroducing CN multiplicative effect and additive effect coming from other source in order to simulate DGE.

```{r}
setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
library(tidyverse)
library(DESeq2)
library(compcodeR)

# Formatting RNA counts dataset #
load("TCGA/breast_cancer/data/brca_rna_normal.Rdata")
brca_rna_normal <- brca_rna_normal[,1:40]
brca_rna_normal <- brca_rna_normal[1:12000,]
letters <- c(1:12000)
rownames(brca_rna_normal) <- letters

# Generate metadata #
metadata <- data.frame(patID = colnames(brca_rna_normal),
                       condition = rep(c("0", "1"), each = 20))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 

```
```{r}
# Construct DESEQDataSet Object #
dds <- DESeq2::DESeqDataSetFromMatrix(countData = brca_rna_normal, 
                              colData = metadata, 
                              design = ~ condition)
dds <- DESeq2::DESeq(dds)

dispersion <- abs(dds@rowRanges@elementMetadata@listData[["dispMAP"]]) %>% na.omit() %>% as.vector()
intercept <- abs(dds@rowRanges@elementMetadata@listData[["Intercept"]]) %>% na.omit() %>% as.vector()
```

```{r}

# compcodeR simulator function: Simulation of normal RNA read counts#

rna_counts_sim <- compcodeR::generateSyntheticData(dataset = "rna_counts_sim", n.vars = 11958, 
                                                   samples.per.cond = 20, n.diffexp = 0, 
                                                   repl.id = 1, seqdepth = 1e7, 
                                                   fraction.upregulated = 0.0, 
                                                   between.group.diffdisp = F, 
                                                   filter.threshold.total = 1, 
                                                   filter.threshold.mediancpm = 0, 
                                                   fraction.non.overdispersed = 0, 
                                                   relmeans = intercept,
                                                   dispersions = dispersion,
                                                   random.outlier.high.prob = 0,
                                                   random.outlier.low.prob = 0,
                                                   output.file = "rna_counts_sim.rds")

rna_counts_sim <- rna_counts_sim@count.matrix
```
```{r}
#`item{n.vars}{number of genes}
#`item{samples.per.cond}{number of replicates per condition}
#`item{relmeans}{vector of the mean of the intercept betas (log2 scale)}
#`item{dispersions}{vector of the standard deviation of the intercept betas (log2 scale)}
```
```{r}
# Generate metadata #
metadata <- data.frame(patID = colnames(rna_counts_sim),
                       condition = rep(c("0", "1"), each = 20))
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var = "patID") 
```
```{r}
dds <- DESeq2::DESeqDataSetFromMatrix(countData = rna_counts_sim, 
                              colData = metadata, 
                              design = ~ condition)
dds <- DESeq2::DESeq(dds)
res <- results(dds, tidy=TRUE) %>% as.tibble()
```

### Volcano plot of RNA simulated counts

```{r, echo=FALSE}
res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange >= 1.0 & res$padj < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange <= -1.0 & res$padj < 0.05] <- "DOWN"

res <- res %>% dplyr::filter(diffexpressed == "NO")
```
```{r, echo=FALSE}
library(gridExtra)
library(ggplot2)
library(colorspace)
library(ggpubr)

p1 <- ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1.0, 1.0), col = "darkgreen", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "darkgreen", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("gray", "gray", "red"))+
  scale_x_continuous(breaks = seq(-10, 10, 2))+
  labs(title="Simulated RNA counts",x="effect size (log2)")+
  theme_bw()+
  theme(legend.position="none")+
  font("xy.text", size = 10, color = "black")+
  font("xlab", size = 10)+
  font("ylab", size = 10)+
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))
p1
```
