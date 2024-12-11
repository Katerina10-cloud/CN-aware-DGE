# Copy-Number-aware Differential Gene Expression (DGE)

## Overview

The presence of biological signals coming from different sources and affecting gene expression make it
desirable to take them into account in differential expression analysis to provide more comprehensive insight regarding
transcriptional patterns of cancer. In DGE analysis between cancer and normal tissues results can be both driven and masked by CNV, potentially leading to false positives and false negatives and could have significant impact on downstream analysis.



## Metodological design
Our CN-aware approach directly integrates Copy Number (CN) data into *Generalized Linear Model (GLM)* statistical framework to model the *mean* of gene counts using *Negative Binomial* distribution as a *log linear function* of the covariates. The model assumes the linear relationship between the CN and expected mean.
This method aims to separate genes whose expression is primarily driven by CNVs from those regulated by independent biological processes. This distinction is crucial in minimizing the confounding effects of CNVs, which can obscure true regulatory changes and lead to misinterpretations in traditional DGE analyses.


$`K \sim NB(mean = \mu_{ij}, dispersion = \alpha_i)`$ 

$`\mu_{ij} = s_{j} (\frac{CN_{ij}}{2}) q_{ij}`$  

$`\log(q_{ij}) = \beta_{i0} X_{j0} + \beta_{i1} X_{j1}`$ 

$`K_{ij} =`$ mRNA counts

$`CN_{ij} =`$ copy number integers (0,1,2,...,5)

$`X_j = [1, X_{j1}] = `$ covariate vector (design matrix, including intercept)

$`\beta_i = [\beta_{i0}, \beta_{i1}] = `$ regression coefficient vector (base level, effect size)

$`s_j =`$ library size factor





### References

1. Michael I Love, Wolfgang Huber, and Simon Anders. Moderated estimation of fold change and dispersion for rna-seq data with deseq2. Genome biology, 15(12):1–21, 2014. doi:10.1186/s13059-014-0550-8.

2. Boris Muzellec, Maria Telenczuk, Vincent Cabeli, and Mathieu Andreux. Pydeseq2: a python package for bulk rna-seq differential expression analysis. bioRxiv, pages 2022–12, 2022. doi:10.1101/2022.12.14.520412.

3. Anqi Zhu, Joseph G Ibrahim, and Michael I Love. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics, 35(12):2084–2092, 2019. doi:10.1093/bioinformatics/bty895.
