rm(list = ls())

library(DESeq2)
library(edgeR)
library(tmaptools)
library(tidyverse)
library(reshape2)

source("/home/sherry/pipline/ML_PDACBiomarker/R/custom_functions.r") # source the built custom functions



data_raw <- read.csv("/home/sherry/pipline/ML_PDACBiomarker/data_R/GSE79668_count.csv", stringsAsFactors = FALSE)
rownames(data_raw) <- data_raw[,1]
data_raw <- data_raw[,-1]
pheno <- read.csv("/home/sherry/pipline/ML_PDACBiomarker/data_R/pheno_GSE79668.csv", stringsAsFactors = FALSE)
head(data_raw)


head(data_raw)

data_clean <- data_raw[, colnames(data_raw) %in% pheno$Sample_ID, ]
data_clean <- na.omit(data_clean)
dim(data_clean)


# Create DGEList object
dgList <- DGEList(data_clean)
# View the first few counts
head(dgList[["counts"]])

# Assign sample groups
samp_groups <- as.factor(pheno$CLASS)
dgList[["samples"]]$group <- samp_groups

# Add additional information to samples table
dgList[["samples"]]$lib.size



######Filtering lowly expressed genes######

#Next,we filter genes based on expression levels and log-transforms counts 
#per million using function filterByExpr. Filtering reduces noise and focuses the analysis on genes with sufficient expression levels.


# Filter genes
keep <- filterByExpr(y = dgList, design = samp_groups)
dgList <- dgList[keep, ]

# Log-transform counts per million (CPM)
log_cpm <- cpm(dgList, log = TRUE)

#####Exploratory Data Analysis####

#Visualisation helps assess the datas distribution and identify potential outliers 
#or patterns.We generate a box plot and PCA plot for visualizing the data before normalisation.


# Box plot of unnormalized data
boxplot(log_cpm, main = "Unnormalized_GSE79668")

# PCA plot
data01 <- t(as.matrix(dgList)) + 1
data01 <- log2(data01)
png(filename = "GSE79668_data_before_normalisation.png", 
    width = 2000, height = 1600, res=200)
PCA_plot(data = data01, group = pheno$CLASS, title = "GSE79668 data before normalisation")
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

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = pheno,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
deseq_sf <- sizeFactors(dds)
head(deseq_sf)
deseq_norm_counts <- counts(dds, normalized = TRUE)


# UQN(Upper Quartile Normalisation)
dgList_uq <- calcNormFactors(dgList, method = "upperquartile")

uq_counts <- cpm(dgList_uq, log = FALSE, normalized.lib.sizes = TRUE)
uq_counts_log <- cpm(dgList_uq, log = TRUE, normalized.lib.sizes = TRUE)



######## Data visualisation ##########

# Box plot of normalized data
png(filename = "BoxPlots/TMM_Normalized_GSE79668.png", width = 2000, height = 1600, res=200)
boxplot(tmm_normalised_counts, main = "TMM Normalized_GSE79668")
dev.off()

png(filename = "BoxPlots/TMM_log_Normalized_GSE79668.png", width = 2000, height = 1600, res=200)
boxplot(tmm_normalised_counts_log, main = "TMM_log Normalized_GSE79668")
dev.off()

png(filename = "BoxPlots/DESeq_Normalized_GSE79668.png", width = 2000, height = 1600, res=200)
boxplot(deseq_norm_counts, main = "DESeq Normalized_GSE79668")
dev.off()

png(filename = "BoxPlots/UQN_GSE79668.png", width = 2000, height = 1600, res=200)
boxplot(uq_counts, main = "UQN_GSE79668")
dev.off()

png(filename = "BoxPlots/UQN_log_GSE79668.png", width = 2000, height = 1600, res=200)
boxplot(deseq_norm_counts, main = "UQN_log_GSE79668")
dev.off()


# PCA plot for TMM normalized data
data02 <- as.matrix(t(tmm_normalised_counts))
png(filename = "PCAPlots/GSE79668_data_normalized_by_TMM.png", width = 2000, height = 1600, res=200)
PCA_plot(data = data02, group = pheno$CLASS, title = "GSE79668 data normalized by TMM")
dev.off()

# PCA plot for TMM_log normalized data
data02 <- as.matrix(t(tmm_normalised_counts_log))
png(filename = "PCAPlots/GSE79668_data_normalized_by_TMM_log.png", width = 2000, height = 1600, res=200)
PCA_plot(data = data02, group = pheno$CLASS, title = "GSE79668 data normalized by TMM_log")
dev.off()

# PCA plot for DESeq normalized data
data_deseq <- as.matrix(t(deseq_norm_counts))
png(filename = "PCAPlots/GSE79668_data_normalized_by_DESeq.png", width = 2000, height = 1600, res=200)
PCA_plot(data = data_deseq, group = pheno$CLASS, title = "GSE79668 data normalized by DESeq")
dev.off()

# PCA plot for UQN normalized data
data_uqn <- as.matrix(t(uq_counts))
png(filename = "PCAPlots/GSE79668_data_normalized_by_UQN.png", width = 2000, height = 1600, res=200)
PCA_plot(data = data_uqn, group = pheno$CLASS, title = "GSE79668 data normalized by UQN")
dev.off()

# PCA plot for UQN_log normalized data
data_uqn_log <- as.matrix(t(uq_counts_log))
png(filename = "PCAPlots/GSE79668_data_normalized_by_UQN_log.png", width = 2000, height = 1600, res=200)
PCA_plot(data = data_uqn_log, group = pheno$CLASS, title = "GSE79668 data normalized by UQN_log")
dev.off()


###### Filtering lowly expressed genes######
keep <- filterByExpr(y = dgList, design = samp_groups)
dgList <- dgList[keep, ]

# Get filtered gene name list
common_genes <- rownames(dgList$counts)
log_cpm_filtered <- log_cpm[common_genes, ]

# TMM
tmm_normalised_counts <- tmm_normalised_counts[common_genes, ]
tmm_normalised_counts_log <- tmm_normalised_counts_log[common_genes, ]

# DESeq
deseq_norm_counts <- deseq_norm_counts[common_genes, ]

# UQN
uq_counts <- uq_counts[common_genes, ]
uq_counts_log <- uq_counts_log[common_genes, ]

rownames(log_cpm_filtered) <- common_genes
rownames(tmm_normalised_counts) <- common_genes
rownames(tmm_normalised_counts_log) <- common_genes
rownames(deseq_norm_counts) <- common_genes
rownames(uq_counts) <- common_genes
rownames(uq_counts_log) <- common_genes

########Save Results#########

save(
  log_cpm_filtered,
  tmm_normalised_counts,        
  tmm_normalised_counts_log,    
  deseq_norm_counts,            
  uq_counts,                    
  uq_counts_log,                
  pheno,                        
  file = "1_TMM_norm.Rdata"
)


# Save normalized data
normalized_data <- data.frame(
  log_cpm = log_cpm_filtered,
  TMM = tmm_normalised_counts,
  TMM_log = tmm_normalised_counts_log,
  DESeq = deseq_norm_counts,
  UQN = uq_counts,
  UQN_log = uq_counts_log
)

write.csv(normalized_data, "normalised_GSE79668.csv")

# Save workspace
save.image("1_TMM_norm.rdata")




