# Load required libraries
library(ggplot2)
library(reshape2)
library(GGally)

# Load data
load("/home/sherry/test/R/normalized.Rdata")

# Extract normalized data
cpm_data <- as.matrix(cpm)
rpkm_data <- as.matrix(rpkm)
tpm_data <- as.matrix(tpm)
tmm_data <- as.matrix(tmm)
DESeqSF_data <- as.matrix(DESeqSF)
uqn_data <- as.matrix(uqn)
log_data <- as.matrix(log)

# Combine data into a data frame
combined_data <- data.frame(
  CPM = as.vector(cpm_data),
  RPKM = as.vector(rpkm_data),
  TPM = as.vector(tpm_data),
  TMM = as.vector(tmm_data),
  DESeq_SF = as.vector(DESeqSF_data),
  UQN = as.vector(uqn_data),
  Log = as.vector(log_data)
)

# Calculate correlation matrix
cor_matrix <- cor(combined_data, method = "pearson")

# Print correlation matrix
print(cor_matrix)

# Convert correlation matrix to long format for heatmap plotting
cor_data <- melt(cor_matrix)

# Plot correlation heatmap
ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#B43970", high = "#282A62", mid = "#FFFFFF", midpoint = 0.5) +
  theme_minimal() +
  labs(title = "Correlation Heatmap of Normalization Methods", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


