# Copy-Number-aware Differential Gene Expression

## Overview

The presence of biological signals coming from different omics and affecting gene expression make it
desirable to take them into account in DE analysis to provide more comprehensive insight regarding
transcriptional patterns of cancer.


## Metodological design
The focus is to combine RNA-seq and Copy Number cancer data in the DE analysis context using *Generalized Linear Model (GLM)* to model the *mean* and *dispersion* of gene counts using *Negative Binomial* distribution as a *log linear function* of the covariates.

$$K \sim NB(mean = \mu_i_j, dispersion = \alpha_i)$$

$$\mu_i_j = s_j q_i_j$$

$$\log(q_i_j) = \beta_i_0 X_j_0 \displaystyle \frac{CN_i_j}{2} + \beta_i_1 X_j_1$$

$K_i_j =$ mRNA counts
$CN_i_j =$ copy number
$X_j = [1, X_j_1]=$ covariate vector (design matrix)
$\beta_i = [\beta_i_0, \beta_i_1]=$ regression coeffizient vector
$s_j =$ library size factor

*Statistical testing (Wald test)*

$$H_0 : \beta_1_i = 0 vs H_1: \beta_1_i \neq 0$$



## Data simulation

RNA counts are simulated using *OmicsSIMLA* simulator using normal gene expression profiles: <https://omicssimla.sourceforge.io/download.html>

Copy Number data are simulated using TCGA cancer data + custom sampling method.

Differential gene expression (tumor-normal) is simulated by introducing CN multiplicative signal into RNA data.

### Example

``` r
cnv_1 <- sapply(1:36, function(x) sample(x=c(0.5,1,2,3), size = 500, replace=TRUE, prob = c(.20, .60, .10, .10)))
cnv_2 <- sapply(1:36, function(x) sample(x=c(1,2,3,4), size = 11000, replace=TRUE, prob = c(.05, .70, .10, .10)))
cnv_normal <- matrix(2, nrow(rna_normal), 36)
cnv <- cbind(cnv_normal, cnv_tumor)

```

RNA counts matrix of normal sample group is multiplied by CN matrix of tumor samples.

### Example

``` r
cnv <- apply(cnv, 1, function(x) x/2)
cnv <- apply(cnv, 1, function(x) x+10e-9)

rna_cnv <- rna_counts * cnv

```

