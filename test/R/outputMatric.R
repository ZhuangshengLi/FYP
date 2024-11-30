rm(list=ls())
options(stringsAsFactors = F)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
library(tximport) # Import transcript-level abundances and counts for transcript- and gene-level analysis packages
library(dplyr)
library(tidyverse) # ggplot2, stringr, dplyr, tidyr, readr, purrr, tibble, forcats
library(data.table) # Reading files with multi-core support
setwd("~/fyp/test/")


#### Processing raw Salmon files ####
## Load the mapping file between transcript_id and symbol
t2s <- fread("~/fyp/test/name_map/t2s_gencode.txt", data.table = F, header = F); head(t2s)

## Find all quant.sf file paths and import Salmon files for processing and aggregation
files <- list.files(pattern="*quant.sf", recursive=T, full.names = T); files  # Display all files in the directory that meet the criteria
txi <- tximport(files, type = "salmon", tx2gene = t2s)

## Extract sample names from the folders as the row names for counts
cn <- sapply(strsplit(files, '\\/'), function(x) x[length(x)-1]); cn
colnames(txi$counts) <- gsub('_quant', '', cn); colnames(txi$counts)

## Extract counts/tpm expression matrix
counts <- as.data.frame(apply(txi$counts, 2, as.integer)) # Convert counts to integers
rownames(counts) <- rownames(txi$counts) 
tpm <- as.data.frame(txi$abundance)  ### Abundance represents gene TPM values
colnames(tpm) <- colnames(txi$counts)


#### Import or create sample information, rename columns, and group ####
b <- read.csv('./SraRunTable.txt')
b
name_list <- b$source_name[match(colnames(counts), b$Run)]; name_list
nlgl <- data.frame(row.names = colnames(counts),
                   name_list = name_list,
                   group_list = name_list)
fix(nlgl)
name_list <- nlgl$name_list
colnames(counts) <- name_list
colnames(tpm) <- name_list

group_list <- nlgl$group_list
gl <- data.frame(row.names = colnames(counts),
                 group_list = group_list)


#### Initial filtering of low expression genes ####
# Filter genes that have counts greater than 1 in at least the number of replicates
keep_feature <- rowSums(counts > 1) >= 2               # ncol(counts) / length(table(group_list)) 
table(keep_feature)  # Check the filtering status
counts_filt <- counts[keep_feature, ] # Replace counts with the filtered gene matrix (keeping higher expressed genes)
tpm_filt <- tpm[keep_feature, ]


#### Save data ####
counts_raw = counts  
counts = counts_filt
tpm = tpm_filt

save(counts_raw, counts, tpm,
     group_list, gl, txi, 
     file='salmon/1.counts.Rdata')

