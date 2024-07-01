# Copy-Number-aware Differential Gene Expression (DGE)

## Overview

The presence of biological signals coming from different omics and affecting gene expression make it
desirable to take them into account in DE analysis to provide more comprehensive insight regarding
transcriptional patterns of cancer.


## Metodological design
The focus is to combine RNA-seq and Copy Number cancer data in the DE analysis context using *Generalized Linear Model (GLM)* to model the *mean* and *dispersion* of gene counts using *Negative Binomial* distribution as a *log linear function* of the covariates.

$`K \sim NB(mean = \mu_{ij}, dispersion = \alpha_i)`$ 

$`\mu_{ij} = s_{j} (\frac{CN_{ij}}{2}) q_{ij}`$  

$`\log(q_{ij}) = \beta_{i0} X_{j0} + \beta_{i1} X_{j1}`$ 

$`K_{ij} =`$ mRNA counts

$`CN_{ij} =`$ copy number integers (0,1,2,...,5)

$`X_j = [1, X_{j1}] = `$ covariate vector (design matrix, including intercept)

$`\beta_i = [\beta_{i0}, \beta_{i1}] = `$ regression coefficient vector (base level, effect size)

$`s_j =`$ library size factor





*Statistical testing (Wald test)*:

$`H_0 : \beta_{1i} = 0`$  vs  $`H_1: \beta_{1i} \neq 0`$



## DE simulation

##### CN saturation (sigmoid curve)

$`\sigma(CN_{ij}) | \log(\frac{CN_{ij}}{2}) = \tau`$

$`\sigma(CN_{ij}) = \frac{2\exp^\tau}{1+\exp^\tau}`$



### References

1. Michael I Love, Wolfgang Huber, and Simon Anders. Moderated estimation of fold change and dispersion for rna-seq data with deseq2. Genome biology, 15(12):1–21, 2014. doi:10.1186/s13059-014-0550-8.

2. Boris Muzellec, Maria Telenczuk, Vincent Cabeli, and Mathieu Andreux. Pydeseq2: a python package for bulk rna-seq differential expression analysis. bioRxiv, pages 2022–12, 2022. doi:10.1101/2022.12.14.520412.

3. Anqi Zhu, Joseph G Ibrahim, and Michael I Love. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics, 35(12):2084–2092, 2019. doi:10.1093/bioinformatics/bty895.
