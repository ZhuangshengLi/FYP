rm(list=ls())
# Load required libraries
library(ggplot2)
library(reshape2)
library(GGally)

# Load data
load("~/fyp/test/R/normalized.Rdata")

# Extract normalized data
cpm_data <- as.matrix(cpm)
rpkm_data <- as.matrix(rpkm)
tpm_data <- as.matrix(tpm)
tmm_data <- as.matrix(tmm)
DESeqSF_data <- as.matrix(DESeqSF)
uqn_data <- as.matrix(uqn)
counts_raw_data <- as.matrix(counts_raw)
counts_data <- as.matrix(counts)
tpm_log_data <- as.matrix(tpm_log)
cpm_log_data <- as.matrix(cpm_log)
rpkm_log_data <- as.matrix(rpkm_log)
tmm_log_data <- as.matrix(tmm_log)
DESeq2_log_data <- as.matrix(DESeq2_log)
uqn_log_data <- as.matrix(uqn_log)
log_data <- as.matrix(log)

# Combine data into a data frame
combined_data <- data.frame(
  CPM = as.vector(cpm_data),
  RPKM = as.vector(rpkm_data),
  TPM = as.vector(tpm_data),
  TMM = as.vector(tmm_data),
  DESeq_SF = as.vector(DESeqSF_data),
  UQN = as.vector(uqn_data),
  Log = as.vector(log_data),
  TPM_Log = as.vector(tpm_log_data),
  CPM_Log = as.vector(cpm_log_data),
  RPKM_Log = as.vector(rpkm_log_data),
  TMM_Log = as.vector(tmm_log_data),
  DESeq2_Log = as.vector(DESeq2_log_data),
  UQN_Log = as.vector(uqn_log_data)
)

# Calculate correlation matrix
cor_matrix <- cor(combined_data, method = "pearson")

# Print correlation matrix
print(cor_matrix)

# Convert correlation matrix to long format for heatmap plotting
cor_data <- melt(cor_matrix)

# Plot correlation heatmap
heatmap_plot <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#B43970", high = "#282A62", mid = "#FFFFFF", midpoint = 0.5) +
  theme_minimal() +
  labs(title = "Correlation Heatmap of Normalization Methods", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),  # 设置图表背景为白色
    plot.background = element_rect(fill = "white", color = NA)    # 设置绘图背景为白色
  )

print(heatmap_plot)

ggsave("./R/correlation_heatmap.png", plot = heatmap_plot, width = 8, height = 6, dpi = 300)



