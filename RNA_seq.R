rm(list = ls())  # This clear the environment of any variables, old plots, that kind of thing.
#fetching files
getwd()
setwd("C:/Users/Acer/Desktop/FILES_LONDON_07_22/2ND_YEAR_BIOCHEMISTRY/Omics_Summer_Term/code")
getwd()

#install.packages("raster","stplanr","tidyverse","sf")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("fgsea")
