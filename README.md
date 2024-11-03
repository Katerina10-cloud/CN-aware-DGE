# Copy-Number-aware Differential Gene Expression (DGE)

## Overview

The presence of biological signals coming from different omics and affecting gene expression make it
desirable to take them into account in DE analysis to provide more comprehensive insight regarding
transcriptional patterns of cancer. In differential gene expression analysis between cancer and normal tissues results can be both driven and masked by CNV, thus leading to false positives and false negatives and could have significant impact on downstream analysis.


## Metodological design
The focus is to directly integrate Copy Number cancer data into  *Generalized Linear Model (GLM)* to model the *mean* and of gene counts using *Negative Binomial* distribution as a *log linear function* of the covariates. The model assumes the linear relationship between the CN and expected mean.


