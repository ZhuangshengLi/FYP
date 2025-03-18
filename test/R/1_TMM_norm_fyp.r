rm(list = ls())

#title: "Tutorial1 :RNA seq data normalisation by TMM normalisation"
#author: "Tanakamol Mahawan"
#date: "2023-12-09"
#output: html_document


# Introduction 

#RNA-seq normalisation is crucial for adjusting raw data to account for factors 
#that can affect the accuracy of analysis. These factors include sequencing depth, 
#transcript length, and variability between samples or batches. Normalisation methods help to ensure that gene expression values are comparable across different samples or conditions, improving the reliability of downstream analyses. However, the best normalisation method can vary depending on the data and research question. It’s a vital step in RNA-seq data analysis. 

###### 1.Load Packages#########

#loads necessary R packages for the tutorial requires functions from these packages for data manipulation, statistical analysis, and visualization. You may need to install some of them if they are not installed.


library(DESeq2)
library(edgeR)
library(tmaptools)
library(tidyverse)
library(reshape2)

source("/home/sherry/pipline/ML_PDACBiomarker/R/custom_functions.r") # source the built custom functions

setwd("/home/sherry/pipline/ML_PDACBiomarker/R/1.TMM_norm")

data_raw <- read.csv("/home/sherry/pipline/ML_PDACBiomarker/R/1.TMM_norm/data/double.csv", 
                     row.name = 1,
                     check.names = FALSE)

# Subset data for selected patients and remove NAs
data_raw <- na.omit(data_raw)
data_clean <- data_raw[rowSums(data_raw) > 0, ]
dim(data_clean)


# Create DGEList object
dgList <- DGEList(data_clean)
keep <- filterByExpr(dgList)
dgList <- dgList[keep, , keep.lib.sizes = FALSE]

# Log-transform counts per million (CPM)
dgList <- calcNormFactors(dgList, method="TMM")
log_cpm <- cpm(dgList, log = TRUE)


boxplot(log_cpm, main = "Before_normalised_GSE154290")


# PCA plot
data01 <- t(as.matrix(dgList)) + 1
data01 <- log2(data01)
png(filename = "GSE154290_data_before_normalisation.png", 
    width = 2000, height = 1600, res=200)
PCA_plot(data = data01, title = "GSE154290 data before normalisation")
dev.off()


#Normalisation using TMM

# Normalize data using TMM
tmm_factors <- calcNormFactors(dgList, method = "TMM")

# Display normalization factors
tmm_factors$samples$norm.factors

# Without log-transform normalized counts
tmm_normalised_counts <- cpm(tmm_factors, log = FALSE, normalised.lib.sizes = TRUE)

# Log-transform normalized counts
tmm_normalised_counts_log <- cpm(tmm_factors, log = TRUE, normalised.lib.sizes = TRUE)


#Normalisation using DESeq2
counts <- dgList$counts

dds <- DESeqDataSetFromMatrix(countData = counts, colData = DataFrame(condition = factor(rep(1, ncol(counts)))), design = ~ 1)
dds <- estimateSizeFactors(dds)
deseq_sf <- sizeFactors(dds)
head(deseq_sf)
deseq_norm_counts <- counts(dds, normalized = TRUE)


# UQN(Upper Quartile Normalisation)
dgList_uq <- calcNormFactors(dgList, method = "upperquartile")

uq_counts <- cpm(dgList_uq, log = FALSE, normalized.lib.sizes = TRUE)
uq_counts_log <- cpm(dgList_uq, log = TRUE, normalized.lib.sizes = TRUE)



# Box plot of normalized data
png(filename = "BoxPlots/TMM_Normalized_GSE154290.png", width = 2000, height = 2000, res=200)
boxplot(tmm_normalised_counts, main = "TMM Normalized_GSE154290")
dev.off()

png(filename = "BoxPlots/TMM_log_Normalized_GSE154290.png", width = 2000, height = 2000, res=200)
boxplot(tmm_normalised_counts_log, main = "TMM_log Normalized_GSE154290")
dev.off()

png(filename = "BoxPlots/DESeq_Normalized_GSE154290.png", width = 2000, height = 2000, res=200)
boxplot(deseq_norm_counts, main = "DESeq Normalized_GSE154290")
dev.off()

png(filename = "BoxPlots/UQN_GSE154290.png", width = 2000, height = 2000, res=200)
boxplot(uq_counts, main = "UQN_GSE154290")
dev.off()

png(filename = "BoxPlots/UQN_log_GSE154290.png", width = 2000, height = 2000, res=200)
boxplot(deseq_norm_counts, main = "UQN_log_GSE154290")
dev.off()


# PCA plot for TMM normalized data
data02 <- as.matrix(t(tmm_normalised_counts))
png(filename = "PCAPlots/GSE154290_data_normalized_by_TMM.png", width = 2000, height = 2000, res=200)
PCA_plot(data = data02, title = "GSE154290 data normalized by TMM")
dev.off()

# PCA plot for TMM_log normalized data
data02 <- as.matrix(t(tmm_normalised_counts_log))
png(filename = "PCAPlots/GSE154290_data_normalized_by_TMM_log.png", width = 2000, height = 2000, res=200)
PCA_plot(data = data02, title = "GSE154290 data normalized by TMM_log")
dev.off()

# PCA plot for DESeq normalized data
data_deseq <- as.matrix(t(deseq_norm_counts))
png(filename = "PCAPlots/GSE154290_data_normalized_by_DESeq.png", width = 2000, height = 2000, res=200)
PCA_plot(data = data_deseq, title = "GSE154290 data normalized by DESeq")
dev.off()

# PCA plot for UQN normalized data
data_uqn <- as.matrix(t(uq_counts))
png(filename = "PCAPlots/GSE154290_data_normalized_by_UQN.png", width = 2000, height = 2000, res=200)
PCA_plot(data = data_uqn, title = "GSE154290 data normalized by UQN")
dev.off()

# PCA plot for UQN_log normalized data
data_uqn_log <- as.matrix(t(uq_counts_log))
png(filename = "PCAPlots/GSE154290_data_normalized_by_UQN_log.png", width = 2000, height = 2000, res=200)
PCA_plot(data = data_uqn_log, title = "GSE154290 data normalized by UQN_log")
dev.off()



########9.Save Results#########


# Save normalized data
normalized_data <- data.frame(
  Gene = rownames(tmm_normalised_counts),
  TMM = tmm_normalised_counts,
  TMM_log = tmm_normalised_counts_log,
  DESeq = deseq_norm_counts,
  UQN = uq_counts,
  UQN_log = uq_counts_log
)

write.csv(normalized_data, "normalised_GSE154290.csv")


# Save workspace
save.image("1_TMM_norm.rdata")

