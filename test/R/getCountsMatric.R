## Extract read count matrix from output of featureCounts
## Read count matrix: raw counts for all genes and samples

### Environment Setup
rm(list=ls())
options(stringsAsFactors = F) 
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr tibble forcats
if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("edgeR")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  BiocManager::install("preprocessCore")
}

if (!requireNamespace("scran", quietly = TRUE)) {
  BiocManager::install("scran")
}

library(preprocessCore)
library(edgeR)
library(DESeq2)
library(scran)
library(limma)
library(data.table) # Multicore file reading
setwd("~/fyp/test")

#### Process counts to obtain expression matrix ####
a1 <- fread('./counts/counts.txt',
            header = T,data.table = F) # Load counts, set the first column as column names
colnames(a1)
counts <- a1[,7:ncol(a1)] # Extract the counts part for gene expression as counts
rownames(counts) <- a1$Geneid # Set gene names as row names
# Rename sample names
colnames(counts)
colnames(counts) <- gsub('/home/test/align/bam/','', # Remove sample name prefix
                         gsub('_sorted.bam','',  colnames(counts))) # Remove sample name suffix

#### Import or construct sample information, rename columns and group samples ####
b <- read.csv('./SraRunTable.txt')
b
name_list <- b$source_name[match(colnames(counts),b$Run)]; name_list  # Select the required sample information column
nlgl <- data.frame(row.names=colnames(counts),
                   name_list=name_list,
                   group_list=name_list)
fix(nlgl)  # Manually edit to construct sample names and grouping information
name_list <- nlgl$name_list
colnames(counts) <- name_list # Update sample names
group_list <- nlgl$group_list
gl <- data.frame(row.names=colnames(counts), # Construct a data frame of sample names and group information
                 group_list=group_list)

#### counts to TPM conversion ####
# Note that the original unscreened counts matrix should be converted
### Extract Geneid and Length (transcript length) from the original counts file counts.txt and calculate TPM
geneid_efflen <- subset(a1,select = c("Geneid","Length"))
colnames(geneid_efflen) <- c("geneid","efflen")  

### Extract efflen corresponding to geneid in counts
efflen <- geneid_efflen[match(rownames(counts),
                              geneid_efflen$geneid),
                        "efflen"]

# log
log_counts <- log2(counts+1)


### Calculate TPM formula
# TPM (Transcripts Per Kilobase Million)
counts2TPM <- function(count=count, efflength=efflen){
  RPK <- count/(efflength/1000)   # Reads Per Kilobase, length normalization
  PMSC_rpk <- sum(RPK)/1e6        # Per million scaling factor, depth normalization
  RPK/PMSC_rpk              
}  

tpm <- as.data.frame(apply(counts,2,counts2TPM))
colSums(tpm)

### Calculate TPM then log-normalize
# TPM followed by log normalization
tpm_log <- log2(tpm + 1)

# Calculate CPM using edgeR
dge <- DGEList(counts = counts)
cpm <- edgeR::cpm(dge)

### Calculate CPM then log-normalize
# CPM followed by log normalization
cpm_log <- log2(cpm + 1)

# Calculate RPKM or FPKM
efflen_kb <- efflen/1000
total_counts <- colSums(counts)/1e6
counts2RPKM <- function(count, efflen_kb, total_counts){
  RPKM <- count/(efflen_kb*total_counts)
  return(RPKM)
}
rpkm <- as.data.frame(mapply(function(col, tc) counts2RPKM(counts[, col], efflen_kb, tc),
                             1:ncol(counts), total_counts))
colnames(rpkm) <- colnames(counts)
rownames(rpkm) <- rownames(counts)

### Calculate RPKM then log-normalize
# RPKM followed by log normalization
rpkm_log <- log2(rpkm + 1)

# Calculate TMM
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
tmm_normalized_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

### Calculate TMM then log-normalize
# TMM followed by log normalization
tmm_log <- log2(tmm_normalized_counts + 1)

# Calculate DESeq2 Size Factor normalization
dds <- DESeqDataSetFromMatrix(countData = counts, colData = DataFrame(condition = factor(rep(1, ncol(counts)))), design = ~ 1)
dds <- estimateSizeFactors(dds)
deseq2_normalized_counts <- counts(dds, normalized = TRUE)

### Calculate DESeq2 Size Factor then log-normalize
# DESeq2 followed by log normalization
deseq2_log <- log2(deseq2_normalized_counts + 1)

# Calculate Upper Quartile Normalization
dge <- calcNormFactors(dge, method = "upperquartile")
upper_quartile_normalized_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

### Calculate Upper Quartile then log-normalize
# Upper Quartile followed by log normalization
uqn_log <- log2(upper_quartile_normalized_counts + 1)

# Merge all duplicate symbols
g2s <- fread('./name_map/g2s_vm25_gencode.txt',header = F,data.table = F) # Load information file extracted from gencode gtf file
colnames(g2s) <- c("geneid","symbol")


symbol <- g2s[match(rownames(counts),g2s$geneid),"symbol"] # Match symbols corresponding to row names in counts
table(duplicated(symbol))  # Count duplicate gene names

### Use aggregate to merge genes with the same symbol in the symbol column
counts <- aggregate(counts, by=list(symbol), FUN=sum)
counts <- column_to_rownames(counts,'Group.1')

tpm <- aggregate(tpm, by=list(symbol), FUN=sum) ### Use aggregate to merge genes with the same symbol in the symbol column
tpm <- column_to_rownames(tpm,'Group.1')

tpm_log <- aggregate(tpm_log, by=list(symbol), FUN=sum)
tpm_log <- column_to_rownames(tpm_log, 'Group.1')

cpm <- aggregate(cpm, by=list(symbol), FUN=sum)
cpm <- column_to_rownames(cpm, 'Group.1')

cpm_log <- aggregate(cpm_log, by=list(symbol), FUN=sum)
cpm_log <- column_to_rownames(cpm_log, 'Group.1')

rpkm <- aggregate(rpkm, by=list(symbol), FUN=sum)
rpkm <- column_to_rownames(rpkm, 'Group.1')

rpkm_log <- aggregate(rpkm_log, by=list(symbol), FUN=sum)
rpkm_log <- column_to_rownames(rpkm_log, 'Group.1')

tmm <- aggregate(tmm_normalized_counts, by=list(symbol), FUN=sum)
tmm <- column_to_rownames(tmm, 'Group.1')

tmm_log <- aggregate(tmm_log, by=list(symbol), FUN=sum)
tmm_log <- column_to_rownames(tmm_log, 'Group.1')

DESeqSF <- aggregate(deseq2_normalized_counts, by=list(symbol), FUN=median)
DESeqSF <- column_to_rownames(DESeqSF, 'Group.1')

DESeq2_log <- aggregate(deseq2_log, by=list(symbol), FUN=median)
DESeq2_log <- column_to_rownames(DESeq2_log, 'Group.1')

uqn <- aggregate(upper_quartile_normalized_counts, by=list(symbol), FUN=median)
uqn <- column_to_rownames(uqn, 'Group.1')

uqn_log <- aggregate(uqn_log, by=list(symbol), FUN=median)
uqn_log <- column_to_rownames(uqn_log, 'Group.1')


log <- aggregate(log_counts, by=list(symbol), FUN=median)
log <- column_to_rownames(log, 'Group.1')

#### Initial filtering of low-expression genes #### (Filtering criteria are not unique and depend on the situation)
# Filter genes (rows) with counts greater than 1 in at least the number of replicates
keep_feature <- rowSums(counts>1) >= 2
table(keep_feature)  # Check filtering results, FALSE for the number of low-expression genes (rows), TRUE for the number of genes to retain

counts_filt <- counts[keep_feature, ] # Replace counts with the filtered gene matrix (retain higher expression genes)
tpm_filt <- tpm[keep_feature, ]
tpm_log_filt <- tpm_log[keep_feature, ]
cpm_filt <- cpm[keep_feature, ]
cpm_log_filt <- cpm_log[keep_feature, ]
rpkm_filt <- rpkm[keep_feature, ]
rpkm_log_filt <- rpkm_log[keep_feature, ]
tmm_filt <- tmm[keep_feature, ]
tmm_log_filt <- tmm_log[keep_feature, ]
DESeqSF_filt <- DESeqSF[keep_feature, ]
DESeq2_log_filt <- DESeq2_log[keep_feature, ]
uqn_filt <- uqn[keep_feature, ]
uqn_log_filt <- uqn_log[keep_feature, ]
log_filt <- log[keep_feature, ]

#### Save data ####
counts_raw=counts # Rename for convenience in subsequent analysis
counts=counts_filt
tpm=tpm_filt
tpm_log=tpm_log_filt
cpm=cpm_filt
cpm_log=cpm_log_filt
rpkm=rpkm_filt
rpkm_log=rpkm_log_filt
tmm=tmm_filt
tmm_log=tmm_log_filt
DESeqSF=DESeqSF_filt
DESeq2_log=DESeq2_log_filt
uqn=uqn_filt
uqn_log=uqn_log_filt
log=log_filt

save(counts_raw, counts,
     tpm, 
     tpm_log,
     cpm,
     cpm_log,
     rpkm,
     rpkm_log,
     tmm,
     tmm_log,
     DESeqSF,
     DESeq2_log,
     uqn,
     uqn_log,
     log,
     group_list, gl, 
     file= "~/fyp/test/R/normalized.Rdata")




