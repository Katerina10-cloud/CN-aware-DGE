### Bayesian network - Survival analysis ###

setwd("/Users/katsiarynadavydzenka/Documents/PhD_AI/")
pkgs <- c("tidyverse", "DESeq2", "bnlearn", "gRain" )
sapply(pkgs, require, character.only = TRUE)

# Prepare data #
